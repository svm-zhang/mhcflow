import multiprocessing as mp
import time
from functools import partial

import numpy as np
import polars as pl
import pysam
from tinyscibio import (
    BAMetadata,
    _PathLike,
    bed_to_df,
    make_dir,
    parse_path,
    walk_bam,
)

from .dumper import _bam2fq_from_idx
from .helper import (
    FileManifest,
    _check_rg_exists,
    _check_single_rg,
    _get_sm,
    _verify_prev_run,
)
from .logger import logger
from .tag_builder import PREBUILD_TYPE, build

CHR6 = ["6", "chr6", "NC00006", "CM000668"]


def _fish_unplaced(
    bam_fspath: _PathLike, prebuilt_tag: PREBUILD_TYPE, out: _PathLike
) -> tuple[pl.DataFrame, _PathLike]:
    logger.info("Fish unplaced sequence with tag pattern")
    out = parse_path(out)
    logdir = out.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    done = logdir / f"{out.stem}.done"
    if done.exists():
        logger.info(
            "Found result of fished sequences with tag pattern "
            f"from previous run: {done}. Skip."
        )
        df = pl.read_csv(out, separator="\t")
        return (df, done)
    qnames = []
    start_t = time.time()
    with pysam.AlignmentFile(str(bam_fspath), "rb") as bamf:
        for aln in bamf.fetch(until_eof=True):
            if not aln.is_unmapped:
                continue
            assert aln.query_sequence is not None
            match = list(prebuilt_tag.iter(aln.query_sequence))
            if match:
                qnames += [aln.query_name]
    logger.info(
        f"Fish unplaced sequence with tag pattern: {time.time() - start_t} sec"
    )
    logger.info(f"Fished {len(qnames)} unplaced sequence with tag pattern")
    merged_qnames = (
        pl.DataFrame({"qnames": qnames}) if qnames else pl.DataFrame()
    )
    merged_qnames.write_csv(out, separator="\t")
    done.touch()
    return merged_qnames, done


def _fish_one_hla(region: str, bam_fspath: str) -> pl.DataFrame:
    return walk_bam(bam_fspath, region, exclude=0, return_qname=True)


def _fish_multi_hla(
    bed_fsapth: _PathLike,
    bam_fspath: _PathLike,
    out: _PathLike,
    nproc: int = 4,
) -> tuple[pl.DataFrame, _PathLike]:
    logger.info(f"Fish sequence mapped to regions defined in {bed_fsapth}.")
    out = parse_path(out)
    logdir = out.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    done = logdir / f"{out.stem}.done"
    if done.exists():
        logger.info(
            "Found result of fished sequences with tag pattern "
            f"from previous run: {done}. Skip."
        )
        df = pl.read_csv(out, separator="\t")
        return (df, done)
    df = bed_to_df(bed_fsapth)
    regions = [f"{row[0]}:{row[1]}-{row[2]}" for row in df.rows()]
    nproc = min(len(regions), nproc)  # nproc set to minimum of these 2 values
    qnames: list[pl.DataFrame] = []
    start_t = time.time()
    with mp.get_context("spawn").Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_fish_one_hla, bam_fspath=str(bam_fspath)),
            regions,
        ):
            if res is not None:
                qnames += [res.select("qnames")]
    logger.info(
        "Fished sequence mapped to regions defined in BED file: "
        f"{time.time() - start_t} sec."
    )
    logger.info(
        f"Fished {len(qnames)} mapped to HLA regions defined in BED file"
    )
    merged_qnames = pl.concat(qnames) if qnames else pl.DataFrame()
    merged_qnames.write_csv(out, separator="\t")
    done.touch()
    return (merged_qnames, done)


def _fish_one_region(
    region: tuple[str, int, int],
    bam_fspath: _PathLike,
    prebuilt_tag: PREBUILD_TYPE,
) -> pl.DataFrame | None:
    qnames = []
    sn, start, end = region
    with pysam.AlignmentFile(str(bam_fspath), "rb") as bamf:
        for aln in bamf.fetch(contig=sn, start=start, stop=end):
            assert aln.query_sequence is not None
            match = list(prebuilt_tag.iter(aln.query_sequence))
            if match:
                qnames += [aln.query_name]
    return pl.DataFrame({"qnames": qnames}) if qnames else None


def _fish_multi_regions(
    split_regions: list[tuple[str, int, int]],
    bam_fspath: _PathLike,
    prebuilt_tag: PREBUILD_TYPE,
    out: _PathLike,
    nproc: int = 4,
) -> tuple[pl.DataFrame, _PathLike]:
    logger.info("Fish sequence with tag pattern.")
    out = parse_path(out)
    logdir = out.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    done = logdir / f"{out.stem}.done"
    if done.exists():
        logger.info(
            "Found result of fished sequences with tag pattern "
            f"from previous run: {done}. Skip."
        )
        df = pl.read_csv(out, separator="\t")
        return (df, done)
    start_t = time.time()
    qnames: list[pl.DataFrame] = []
    with mp.Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(
                _fish_one_region,
                bam_fspath=bam_fspath,
                prebuilt_tag=prebuilt_tag,
            ),
            split_regions,
        ):
            if res is not None:
                qnames += [res]
    logger.info(f"Fish sequence with tag pattern: {time.time() - start_t} sec")
    logger.info(
        f"Fished {sum([df.shape[0] for df in qnames])} readss with tag pattern."
    )
    # merged_qnames will be empty dataframe if there were no df returned
    merged_qnames = pl.concat(qnames) if qnames else pl.DataFrame()
    # empty file will be generated when merged_qnames is empty df
    merged_qnames.write_csv(out, separator="\t")
    done.touch()
    return (merged_qnames, done)


def _split_regions(
    regions: dict[str, list[int]],
    by: str,
    n_splits: int = 4,
    split_size: int = 500_000,
) -> list[tuple[str, int, int]]:
    if by not in ["len", "num"]:
        raise ValueError(f"Unsupported split by method {by=}")
    logger.info("Split regions into smaller intevals.")
    splits = []
    for sn, region in regions.items():
        start, end = region
        n_split = (
            ((end - start + 1) // split_size) + 1 if by == "len" else n_splits
        )
        s = np.array_split(np.arange(end), n_split)
        splits += [(sn, int(split[0]), int(split[-1])) for split in s]

    logger.info(f"Split regions into {len(splits)} intevals.")
    return splits


def _run_fisher(
    bam_fspath: _PathLike,
    tag_fspath: _PathLike,
    bed_fspath: _PathLike,
    outdir: _PathLike,
    prebuild_method: str = "ahocorasick",
    nproc: int = 4,
    overwrite: bool = False,
) -> FileManifest:
    logger.info("Start to fish HLA-relevant reads.")
    outdir = parse_path(outdir)
    make_dir(outdir, parents=True, exist_ok=True)

    bametadata = BAMetadata(str(bam_fspath))
    _check_rg_exists(bametadata)
    _check_single_rg(bametadata)
    rg = bametadata.read_groups[0]
    sm = _get_sm(rg)

    fm_json = outdir / f"{sm}.fisher.file_manifest.json"
    if fm_json.exists():
        logger.info(
            f"Detected file manifest from previous fisher run: {str(fm_json)}"
        )
        fisher_fm = FileManifest._from_json(fm_json)
        if _verify_prev_run(fisher_fm, overwrite):
            return fisher_fm

    fisher_fm = FileManifest()
    logdir = outdir / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    fisher_done = logdir / f"{sm}.fisher.done"

    prebuilt_tag = build(tag_fspath, method=prebuild_method)

    bam_fspath = bametadata.fspath
    hla_qname_out = outdir / f"{sm}.fisher.hla_bed.idx"
    hla_bed_qnames, hla_bed_done = _fish_multi_hla(
        bed_fspath, bam_fspath, hla_qname_out
    )

    regions = {
        sn: [1, ln] for sn, ln in bametadata.seqmap().items() if sn in CHR6
    }
    splits = _split_regions(regions, n_splits=nproc, by="num")
    chr6_qname_out = outdir / f"{sm}.fisher.chr6.idx"
    chr6_qnames, chr6_done = _fish_multi_regions(
        splits, bam_fspath, prebuilt_tag, chr6_qname_out, nproc
    )
    unplaced_qname_out = outdir / f"{sm}.fisher.unplaced.idx"
    unplaced_qnames, unplaced_done = _fish_unplaced(
        bam_fspath, prebuilt_tag, unplaced_qname_out
    )

    merged_qnames = pl.concat(
        [hla_bed_qnames, chr6_qnames, unplaced_qnames]
    ).unique()
    fisher_idx_out = outdir / f"{sm}.fisher.idx.final.tsv"
    merged_qnames.write_csv(fisher_idx_out, separator="\t")

    dumper_fm = _bam2fq_from_idx(fisher_idx_out, bametadata, outdir, nproc)

    logger.info(f"Fished {merged_qnames.shape[0]} in total.")
    fisher_done.touch()

    logger.info("Register all relevant files to manifest.")
    # register all relevant files to manifest
    fisher_fm._register_inputs(
        bam=bametadata.fspath, tag=tag_fspath, bed=bed_fspath
    )
    fisher_fm._register_aux(done=fisher_done, myself=fm_json)
    fisher_fm._register_outputs(
        **dumper_fm.outputs,
        hla_fished_idx=hla_qname_out,
        chr6_fished_idx=chr6_qname_out,
        unplaced_fished_idx=unplaced_qname_out,
        fished_all_idx=fisher_idx_out,
    )
    fisher_fm._register_intermediate(**dumper_fm.outputs)
    fisher_fm._register_intermediate_aux(
        **dumper_fm.intermediate_aux,
        **dumper_fm.aux,
        hla_bed_done=hla_bed_done,
        chr6_done=chr6_done,
        unplaced_done=unplaced_done,
    )
    # dump manifest to disk in json format.
    logger.info(f"Dump manifest to {fm_json}")
    fisher_fm._to_json(json_out=fm_json)

    return fisher_fm
