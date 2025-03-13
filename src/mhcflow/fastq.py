import multiprocessing as mp
import shlex
import subprocess as sp
from functools import partial
from pathlib import Path

import numpy as np
import polars as pl
from tinyscibio import BAMetadata, _PathLike, parse_path


def _extract_from_bam(batch_fspath: tuple[Path, Path, Path], bam_fspath: str):
    qname_batch_fspath, out_r1, out_r2 = batch_fspath
    qname_batch_fspath = str(qname_batch_fspath)
    try:
        cmd_1 = f"samtools view -h -N {qname_batch_fspath} {bam_fspath}"
        p1 = sp.Popen(shlex.split(cmd_1), stdout=sp.PIPE)
        cmd_2 = "samtools sort -n"
        p2 = sp.Popen(shlex.split(cmd_2), stdin=p1.stdout, stdout=sp.PIPE)
        cmd_3 = f"samtools fastq -n -1 {out_r1} -2 {out_r2} -0 /dev/null -s /dev/null"
        print(" | ".join([cmd_1, cmd_2, cmd_3]))
        p3 = sp.Popen(shlex.split(cmd_3), stdin=p2.stdout)
        p3.communicate()
        p1.wait()
        p2.wait()
    except Exception as e:
        print(e)


def dump_to_fastq(
    bametadata: BAMetadata,
    qnames: pl.Series,
    outdir: _PathLike,
    nproc: int = 1,
) -> list[tuple]:
    # FIXME: move this to main
    outdir = parse_path(outdir)
    qname_batches = np.array_split(qnames.to_numpy(), nproc)
    qname_batch_fspaths = []
    for i in range(len(qname_batches)):
        qname_batch_fspath = outdir / f"qnames.split.{i}.ids.txt"
        batch_r1 = outdir / f"qnames.split.{i}.R1.fastq"
        batch_r2 = outdir / f"qnames.split.{i}.R2.fastq"
        pl.DataFrame({"qnames": qname_batches[i]}).write_csv(
            qname_batch_fspath, include_header=False
        )
        qname_batch_fspaths.append((qname_batch_fspath, batch_r1, batch_r2))

    print(len(qname_batches))
    with mp.Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_extract_from_bam, bam_fspath=bametadata.fspath),
            qname_batch_fspaths,
        ):
            pass
    return qname_batch_fspaths
