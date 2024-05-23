from __future__ import annotations

import argparse
from dataclasses import dataclass
import itertools
import math
import numpy as np
import polars as pl
import pysam
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

    return parser


def get_alleles(bam: str) -> list(str):
    bamf = pysam.AlignmentFile(bam, "rb")
    return bamf.references


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


def extract_alignments(bam: str, allele: str, freq_df: pl.DataFrame):
    bamf = pysam.AlignmentFile(bam, "rb")
    n_reads = bamf.count(contig=allele)
    if n_reads < 50:
        return

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
    res_df = res_df.groupby("ids").agg(pl.col("scores").sum())
    res_df = res_df.with_columns(allele=pl.lit(allele))
    # FIXME: need to add frequency prior to the log-liklihood score
    tot_df = pl.DataFrame(
        {
            "ids": ["total"],
            "scores": [res_df["scores"].sum()],
            "allele": [allele],
        }
    )
    res_df = pl.concat([res_df, tot_df])
    res_df = res_df.with_columns(pl.col("scores").cast(pl.Float32))
    print(allele)
    # print(res_df.sort(by="ids").head(10))
    # print(res_df.filter(pl.col("ids") == "total"))
    # if allele == "hla_a_01_01_02":
    #    sys.exit()


def main():

    parser = parse_cmd()
    args = parser.parse_args()

    freq_df = pl.read_csv(args.freq, separator="\t")
    print(freq_df.head())

    alleles = get_alleles(bam=args.bam)
    print(alleles[0:5])

    for allele in alleles:
        extract_alignments(bam=args.bam, allele=allele, freq_df=freq_df)


if __name__ == "__main__":
    main()
