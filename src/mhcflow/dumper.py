import multiprocessing as mp
from functools import partial

import numpy as np
import polars as pl
from tinyscibio import BAMetadata, _PathLike, make_dir, parse_path

from .helper import (
    FileManifest,
    _check_rg_exists,
    _check_single_rg,
    _get_sm,
    _verify_prev_run,
)
from .logger import logger
from .runnable import _samtools_fastq


def _bam2fq_from_idx(
    idx_fspath: _PathLike,
    bametadata: BAMetadata,
    outdir: _PathLike,
    nproc: int = 1,
    overwrite: bool = False,
) -> FileManifest:
    logger.info(f"Start to dump reads to fastq given ids: {idx_fspath}")
    # set up output dir
    outdir = parse_path(outdir)
    make_dir(outdir, parents=True, exist_ok=True)

    # get sm field from rg
    _check_rg_exists(bametadata)
    _check_single_rg(bametadata)
    rg = bametadata.read_groups[0]
    sm = _get_sm(rg)

    fm_json = outdir / f"{sm}.dumper.file_manifest.json"
    if fm_json.exists():
        logger.info(
            f"Detected file manifest from previous dumper run: {str(fm_json)}"
        )
        fm = FileManifest._from_json(fm_json)
        if _verify_prev_run(fm, overwrite):
            return fm

    fm = FileManifest()
    logdir = outdir / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    dumper_done = logdir / f"{sm}.dumper.done"
    idx_df = pl.read_csv(idx_fspath, separator="\t")
    if "qnames" not in idx_df.columns:
        raise pl.exceptions.ColumnNotFoundError()

    # split read ids into batches and output each batch to file
    # each idx file will be input to samtools view -N
    qnames = idx_df["qnames"].unique().to_numpy()
    qnames_batches = np.array_split(qnames, nproc)
    idxs = []
    for i in range(len(qnames_batches)):
        qname_batch_fspath = outdir / f"{sm}.fisher.idxs.{i}.txt"
        pl.DataFrame({"qnames": qnames_batches[i]}).write_csv(
            qname_batch_fspath, include_header=False
        )
        idxs.append(qname_batch_fspath)
    # run samtools fastq to extract read and dump to fq
    r1s, r2s = [], []
    dones, logs = [], []
    with mp.Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_samtools_fastq, bam_fspath=bametadata.fspath), idxs
        ):
            r1, r2, done, log = res
            r1s.append(r1)
            r2s.append(r2)
            dones.append(done)
            logs.append(log)
    dumper_done.touch()
    fm._register_inputs(idx=idx_fspath, bam=bametadata.fspath)
    fm._register_outputs(r1s=r1s, r2s=r2s, idxs=idxs)
    fm._register_aux(done=dumper_done)
    fm._register_intermediate_aux(dones=dones, logs=logs)
    fm._to_json(fm_json)
    return fm
