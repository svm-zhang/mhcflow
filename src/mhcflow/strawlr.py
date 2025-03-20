import multiprocessing as mp
import shutil
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

from .helper import (
    FileManifest,
    _check_rg_exists,
    _check_single_rg,
    _get_sm,
)
from .logger import logger
from .runnable import _extract_from_bam
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
    bed_fsapth: _PathLike, bam_fspath: _PathLike, nproc: int = 4
) -> list[pl.Series]:
    logger.info(f"Fish sequence mapped to regions defined in {bed_fsapth}.")
    df = bed_to_df(bed_fsapth)
    regions = [f"{row[0]}:{row[1]}-{row[2]}" for row in df.rows()]
    nproc = min(len(regions), nproc)  # nproc set to minimum of these 2 values
    qnames: list[pl.Series] = []
    start_t = time.time()
    with mp.get_context("spawn").Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_fish_one_hla, bam_fspath=str(bam_fspath)),
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

    fisher_fm = FileManifest()
    fm_json = outdir / f"{sm}.fisher.file_manifest.json"
    if fm_json.exists():
        # load json and check existence of done
        if not overwrite:
            fisher_fm = fisher_fm._from_json(fm_json)
            fisher_done = parse_path(fisher_fm.aux.get("done", ""))
            if fisher_done.exists():
                logger.info(
                    "Found done file for fisher from previous run: "
                    f"{fisher_done}. Skip."
                )
                return fisher_fm
        logger.info(
            f"Overwrite specified. Remove results from previous run: {outdir}"
        )
        # TODO: implement a _clean method to FileManifest to avoid below
        shutil.rmtree(outdir)
        make_dir(outdir, parents=True, exist_ok=True)
    fisher_done = outdir / f"{sm}.fisher.done"

    prebuilt_tag = build(tag_fspath, method=prebuild_method)

    bam_fspath = bametadata.fspath
    fished_qnames = _fish_multi_hla(bed_fspath, bam_fspath)

    regions = {
        sn: [1, ln] for sn, ln in bametadata.seqmap().items() if sn in CHR6
    }
    splits = _split_regions(regions, n_splits=nproc, by="num")
    fished_qnames += _fish_multi_regions(
        splits, bam_fspath, prebuilt_tag, nproc
    )
    fished_qnames += [_fish_unplaced(bam_fspath, prebuilt_tag)]
    merged_qnames = pl.concat(fished_qnames).unique()

    # split all fished read names into batches
    qnames_batches = np.array_split(merged_qnames.to_numpy(), nproc)
    idxs = []
    for i in range(len(qnames_batches)):
        qname_batch_fspath = outdir / f"{sm}.fisher.{i}.idxs"
        pl.DataFrame({"qnames": qnames_batches[i]}).write_csv(
            qname_batch_fspath, include_header=False
        )
        idxs.append(qname_batch_fspath)

    # extract reads to fastq given read names
    r1s, r2s = [], []
    with mp.Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_extract_from_bam, bam_fspath=bametadata.fspath), idxs
        ):
            r1, r2 = res
            r1s.append(r1)
            r2s.append(r2)
    logger.info(f"Fished {merged_qnames.shape[0]} in total.")
    fisher_done.touch()

    logger.info("Register all relevant files to manifest.")
    # register all relevant files to manifest
    fisher_fm._register_inputs(
        bam=bametadata.fspath, tag=tag_fspath, bed=bed_fspath
    )
    fisher_fm._register_aux(done=fisher_done)
    fisher_fm._register_outputs(idxs=idxs, r1s=r1s, r2s=r2s)
    fisher_fm._register_intermediate(r1s=r1s, r2s=r2s)
    # dump manifest to disk in json format.
    logger.info(f"Dump manifest to {fm_json}")
    fisher_fm._register_aux(myself=fm_json)
    fisher_fm._to_json(json_out=fm_json)

    return fisher_fm
