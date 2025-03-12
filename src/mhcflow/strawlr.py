import multiprocessing as mp
import time
from functools import partial

import numpy as np
import polars as pl
import pysam
from tinyscibio import BAMetadata, _PathLike, bed_to_df, walk_bam

from .logger import logger
from .tag_builder import build

CHR6 = ["6", "chr6", "NC00006", "CM000668"]


def _fish_unplaced(bam_fspath, prebuilt_tag) -> pl.Series:
    logger.info("Fish unplaced sequence with tag pattern")
    qnames = []
    start_t = time.time()
    with pysam.AlignmentFile(bam_fspath, "rb") as bamf:
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
    return pl.Series("qnames", qnames)


def _fish_one_hla(region: str, bam_fspath: str) -> pl.DataFrame:
    return walk_bam(bam_fspath, region, exclude=0, return_qname=True)


def _fish_multi_hla(
    bed_fsapth: _PathLike, bam_fspath: str, nproc: int = 4
) -> list[pl.Series]:
    logger.info(f"Fish sequence mapped to regions defined in {bed_fsapth}.")
    df = bed_to_df(bed_fsapth)
    regions = [f"{row[0]}:{row[1]}-{row[2]}" for row in df.rows()]
    print(regions)
    nproc = min(len(regions), nproc)  # nproc set to minimum of these 2 values
    qnames: list[pl.Series] = []
    start_t = time.time()
    with mp.get_context("spawn").Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_fish_one_hla, bam_fspath=bam_fspath),
            regions,
        ):
            if res is not None:
                qnames += [res["qnames"]]
    logger.info(
        "Fished sequence mapped to regions defined in BED file: "
        f"{time.time() - start_t} sec."
    )
    logger.info(
        f"Fished {len(qnames)} mapped to HLA regions defined in BED file"
    )
    return qnames


def _fish_one_region(region, bam_fspath, prebuilt_tag) -> pl.Series | None:
    qnames = []
    sn, start, end = region
    with pysam.AlignmentFile(str(bam_fspath), "rb") as bamf:
        for aln in bamf.fetch(contig=sn, start=start, stop=end):
            assert aln.query_sequence is not None
            match = list(prebuilt_tag.iter(aln.query_sequence))
            if match:
                qnames += [aln.query_name]
    return pl.Series("qnames", qnames) if qnames else None


def _fish_multi_regions(
    split_regions: list[tuple[str | int]],
    bam_fspath,
    prebuilt_tag,
    nproc: int = 4,
) -> list[pl.Series]:
    logger.info("Fish sequence with tag pattern.")
    start_t = time.time()
    qnames: list[pl.Series] = []
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
    logger.info(f"Fished {len(qnames)} sequences with tag pattern.")
    return qnames


def _split_regions(
    regions: dict[str, list[int]],
    by: str,
    n_splits: int = 4,
    split_size: int = 500_000,
) -> list[tuple[str | int]]:
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


def run_strawlr(
    bam_fspath: str,
    tag_seq_fspath: _PathLike,
    hla_bed_fspath: _PathLike,
    prebuild_method: str = "ahocorasick",
    nproc: int = 4,
):
    prebuilt_tag = build(tag_seq_fspath, method=prebuild_method)

    bametadata = BAMetadata(bam_fspath)

    fished_qnames = _fish_multi_hla(hla_bed_fspath, bam_fspath)

    regions = {
        sn: [1, ln] for sn, ln in bametadata.seqmap().items() if sn in CHR6
    }
    splits = _split_regions(regions, n_splits=nproc, by="num")
    fished_qnames += _fish_multi_regions(
        splits, bam_fspath, prebuilt_tag, nproc
    )
    fished_qnames += [_fish_unplaced(bam_fspath, prebuilt_tag)]
    merged_qnames = pl.concat(fished_qnames).unique()

    logger.info(f"Fished {merged_qnames.shape[0]} in total.")
