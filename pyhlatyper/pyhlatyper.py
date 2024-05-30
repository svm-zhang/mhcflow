#!/usr/bin/env python

from __future__ import annotations

import argparse
import itertools
import math
import numpy as np
import os
import polars as pl
import pysam
import re
import sys


def parse_cmd() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        metavar="FILE",
        type=str,
        required=True,
        help="specify path to BAM file",
    )
    parser.add_argument(
        "--freq",
        metavar="FILE",
        type=str,
        required=True,
        help="specify path to HLA frequency file",
    )
    parser.add_argument(
        "--race",
        metavar="STR",
        type=str,
        default="Unknown",
        help="specify path to HLA frequency file",
    )
    parser.add_argument(
        "--outdir",
        metavar="DIR",
        type=str,
        required=True,
        help="specify path to output directory",
    )

    return parser


def get_alleles(bam: str) -> list(str):
    bamf = pysam.AlignmentFile(bam, "rb")
    return bamf.references


def get_rg(bam: str) -> dict(str, str):
    bamh = pysam.AlignmentFile(bam, "rb")
    rg = bamh.header.get("RG", None)
    if rg is None:
        raise ValueError("[ERROR] Found no RG information in BAM")

    return rg[0]


def parse_cigar(cigar_str: str):
    cigar_iter = itertools.groupby(cigar_str, lambda k: k.isdigit())
    cigar_list = []
    for _, n in cigar_iter:
        op = int("".join(n)), "".join(next(cigar_iter)[1])
        cigar_list.append(op)
    return cigar_list


def parse_md(md_str: str):
    md_iter = itertools.groupby(
        md_str, lambda k: k.isalpha() or not k.isalnum()
    )
    return ["".join(group) for c, group in md_iter if not c or group]


def aln_has_indel(cigar_str: str):
    cigar = set(k[1] for k in parse_cigar(cigar_str=cigar_str))
    return "I" in cigar or "D" in cigar or "S" in cigar


def score_log_liklihood(base_qs, md, scale=math.exp(23)):
    score: float = 0.0
    start, end = 0, 0
    for i in range(len(md)):
        # print(start, end, md[i])
        if md[i].startswith("^"):
            continue
        block = np.array([])
        if md[i].isalpha():
            # mismatches
            end = start + len(md[i])
            block = np.array([10 ** (-k / 10) / 3 for k in base_qs[start:end]])
        else:
            # matches
            end = start + int(md[i])
            block = np.array([1 - 10 ** (-k / 10) for k in base_qs[start:end]])
        score += np.sum(np.log(np.multiply(block, scale)))
        start = end
        # print(start, end, score)
    return score


def extract_gene_from_allele(allele: str) -> str:
    allele = allele.lower()
    # A*01:01:01, hla_a_01_01_01
    return re.sub("(\\*|_)*([0-9]+[a-z]*(:|_)*)", "", allele)


def extract_supertype_from_allele(allele: str) -> str:
    allele = allele.lower()
    a_iter = itertools.groupby(allele, lambda k: k.isalnum())
    return "_".join(
        ["".join(group) for i, (c, group) in enumerate(a_iter) if c and i < 8]
    )


def extract_alignments(
    bam: str, allele: str, freq_df: pl.DataFrame
) -> pl.DataFrame | None:
    hla_gene = extract_gene_from_allele(allele)
    bamf = pysam.AlignmentFile(bam, "rb")
    n_reads = bamf.count(contig=allele)
    if n_reads < 50:
        return None

    scores: list[float] = []
    ids: list[str] = []
    for aln in bamf.fetch(contig=allele):
        if aln.is_qcfail or aln.is_supplementary:
            continue

        if not aln.is_proper_pair:
            continue

        if aln_has_indel(cigar_str=aln.cigarstring):
            continue

        if not aln.has_tag("MD"):
            continue
        md_str = aln.get_tag("MD")
        # md_str = "12G6C7^A5"
        md = parse_md(md_str=md_str)

        score = score_log_liklihood(aln.query_qualities, md)
        scores += [np.float32(score)]
        ids += [aln.qname]

    res_df = pl.DataFrame(
        {
            "ids": ids,
            "scores": scores,
        },
    )
    # this makes sure only paired reads counted towards final score
    res_df = res_df.with_columns(pl.col("ids").count().over("ids").alias("n"))
    res_df = res_df.filter(pl.col("n") == 2)
    res_df = res_df.drop("n")
    # FIXME: need to add frequency prior to the log-liklihood score
    res_df = res_df.group_by("ids").agg(pl.col("scores").sum())
    res_df = res_df.with_columns(allele=pl.lit(allele), gene=pl.lit(hla_gene))
    print(allele)
    return res_df


def get_winners(allele_scores: pl.DataFrame) -> pl.DataFrame:
    tot_scores = allele_scores.group_by(["allele", "gene"]).agg(
        pl.col("scores").sum()
    )
    winners = tot_scores.filter(
        pl.col("scores") == pl.col("scores").max().over("gene")
    )
    winners = winners.rename({"scores": "tot_scores"})
    return winners


def score_second(score_table: pl.DataFrame, a1_winners: pl.DataFrame):

    score_table = score_table.join(a1_winners, on=["ids", "gene"], how="left")
    score_table = score_table.with_columns(
        pl.col("scores_right").fill_null(0.0)
    )
    score_table = score_table.with_columns(
        factor=pl.col("scores") / (pl.col("scores") + pl.col("scores_right"))
    )
    score_table = score_table.with_columns(
        scores=pl.col("scores") * pl.col("factor")
    )
    return score_table


def main():

    parser = parse_cmd()
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    freq_df = pl.read_csv(args.freq, separator="\t")
    # remove supertype has all zero pop frequency
    freq_df = freq_df.filter(
        pl.fold(0, lambda acc, s: acc + s, pl.all().exclude(pl.String)) > 0.0
    )

    rg = get_rg(bam=args.bam)
    sid = rg.get("SM")
    alleles = get_alleles(bam=args.bam)
    alleles_df = pl.DataFrame({"Allele": alleles})
    alleles_df = alleles_df.with_columns(
        pl.col("Allele")
        .map_elements(extract_supertype_from_allele, return_dtype=pl.String)
        .alias("SuperType")
    )
    alleles_df = alleles_df.join(
        freq_df, left_on="SuperType", right_on="Allele", how="inner"
    )
    alleles = alleles_df["Allele"].to_list()

    out_a1 = f"{args.outdir}/{sid}.a1.tsv"
    out_a2 = f"{args.outdir}/{sid}.a2.tsv"

    a1_scores = pl.DataFrame()
    if not os.path.exists(out_a1):
        print("score first allele")
        score_tables: list[pl.DataFrame] = []
        for allele in alleles:
            score_tables.append(
                extract_alignments(
                    bam=args.bam, allele=allele, freq_df=freq_df
                )
            )
        a1_scores = pl.concat([s for s in score_tables if s is not None])
        a1_scores.write_csv(out_a1, separator="\t")

    a1_scores = pl.read_csv(out_a1, separator="\t")

    a1_winners = get_winners(allele_scores=a1_scores)
    winner_scores = a1_scores.join(
        a1_winners, on=["gene", "allele"], how="inner"
    )
    if not os.path.exists(out_a2):
        print("score second allele")
        a2_score_tables: list[pl.DataFrame] = []
        for gene in a1_winners["gene"].unique().to_list():
            a2_score_tables.append(
                score_second(
                    a1_scores.filter(pl.col("gene") == gene),
                    winner_scores.filter(pl.col("gene") == gene),
                )
            )
        a2_scores = pl.concat(a2_score_tables)
        a2_scores.write_csv(out_a2, separator="\t")

    a2_scores = pl.read_csv(out_a2, separator="\t")
    a2_winners = get_winners(allele_scores=a2_scores)

    hla_res = f"{args.outdir}/{sid}.hlatyping.res.tsv"
    hla_res_df = pl.concat([a1_winners, a2_winners])
    hla_res_df = hla_res_df.with_columns(sample=pl.lit(sid)).sort(by="allele")
    print(hla_res_df)
    hla_res_df.write_csv(hla_res, separator="\t")


if __name__ == "__main__":
    main()
