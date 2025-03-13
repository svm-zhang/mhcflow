import multiprocessing as mp
import shlex
import subprocess as sp
from functools import partial
from pathlib import Path

import numpy as np
import polars as pl
from tinyscibio import BAMetadata, _PathLike, parse_path

from .logger import logger


def _extract_from_bam(
    batch_fspath: tuple[int, Path, Path, Path, Path], bam_fspath: str
):
    logger.initialize()
    _, qname_batch_fspath, done_fspath, out_r1, out_r2 = batch_fspath
    if done_fspath.exists():
        logger.info(
            "Extraction of reads into fastq has been done "
            f"for batch {qname_batch_fspath}."
        )
        return
    qname_batch_fspath = str(qname_batch_fspath)
    try:
        cmd_1 = f"samtools view -h -N {qname_batch_fspath} {bam_fspath}"
        p1 = sp.Popen(shlex.split(cmd_1), stdout=sp.PIPE)
        cmd_2 = "samtools sort -n"
        p2 = sp.Popen(shlex.split(cmd_2), stdin=p1.stdout, stdout=sp.PIPE)
        cmd_3 = f"samtools fastq -n -1 {out_r1} -2 {out_r2} -0 /dev/null -s /dev/null"
        merged_cmd = " | ".join([cmd_1, cmd_2, cmd_3])
        logger.info(f"Extract reads into fastq using cmd: {merged_cmd}")
        p3 = sp.Popen(
            shlex.split(cmd_3),
            stdin=p2.stdout,
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL,
        )
        p3.communicate()
        p1.wait()
        p2.wait()
        done_fspath.touch()
    except Exception as e:
        print(e)


def _dispatch_qnames(
    qnames: pl.Series, outdir: _PathLike, nproc: int = 1
) -> list[tuple]:
    outdir = parse_path(outdir)
    qname_batch_fspaths = list(outdir.glob("*split*ids.txt"))
    if len(qname_batch_fspaths) != 0 and len(qname_batch_fspaths) != nproc:
        logger.info(
            "Detected inconsistent qname batching "
            f"between previous ({len(qname_batch_fspaths)}) and "
            f"current ({nproc}) run."
        )
        logger.info(f"Respect the current batching setting: {nproc}.")
        logger.info("Remove batching files from previous run.")
        prev_batch_fq_fspaths = list(outdir.glob("*split*fastq"))
        prev_batch_done_fspath = list(outdir.glob("*done$"))
        _ = [
            p.unlink(missing_ok=True)
            for p in qname_batch_fspaths
            + prev_batch_fq_fspaths
            + prev_batch_done_fspath
        ]
    qname_batches = np.array_split(qnames.to_numpy(), nproc)
    qname_batch_fspaths = []
    for i in range(len(qname_batches)):
        qname_batch_fspath = outdir / f"qnames.split.{i}.ids.txt"
        qname_batch_done_fspath = outdir / f"qnames.split.{i}.ids.done"
        batch_r1 = outdir / f"qnames.split.{i}.R1.fastq"
        batch_r2 = outdir / f"qnames.split.{i}.R2.fastq"
        pl.DataFrame({"qnames": qname_batches[i]}).write_csv(
            qname_batch_fspath, include_header=False
        )
        qname_batch_fspaths.append(
            (
                i,
                qname_batch_fspath,
                qname_batch_done_fspath,
                batch_r1,
                batch_r2,
            )
        )
    return qname_batch_fspaths


def dump_to_fastq(
    bametadata: BAMetadata,
    qname_batches: list[tuple],
    nproc: int = 1,
) -> None:
    with mp.Pool(processes=nproc) as pool:
        for _ in pool.imap_unordered(
            partial(_extract_from_bam, bam_fspath=bametadata.fspath),
            qname_batches,
        ):
            pass
    for qname_batch in qname_batches:
        batch_id, qname_batch_fspath, done_fspath, r1, r2 = qname_batch
        if (not r1.exists() or not r2.exists()) and not done_fspath.exists():
            raise FileNotFoundError(
                f"Failed to find either {r1} or {r2} fastq file. "
                f"Extraction failed for qnames in batch {batch_id}: "
                f"{str(qname_batch_fspath)}."
            )
