import subprocess as sp
from pathlib import Path

from .logger import logger


def _novoalign(manifest, ref, rg_sm):
    _, qname_batch_fspath, done_fspath, r1, r2 = manifest
    bam = r1
    try:
        cmd_1 = [
            "novoalign",
            "-d",
            ref,
            "-F",
            "STDFQ",
            "-R",
            "0",
            "-r",
            "All",
            "-o",
            "FullNW",
            "-o",
            "SAM",
            rg_sm,
            "-f",
            r1,
            r2,
        ]
        p1 = sp.Popen(cmd_1, stdout=sp.PIPE)
        cmd_2 = ["samtools", "view", "-bh", "-o", bam]
        p2 = sp.Popen(cmd_2, stdin=p1.stdout, stdout=sp.PIPE)
        cmd_str = " | ".join([" ".join(cmd_1), " ".join(cmd_2)])
        logger.info(f"Extract reads into fastq using cmd: {cmd_str}")
        p2.communicate()
        p1.wait()
    except Exception as e:
        print(e)


def _concat(bam_list_fspath, out_bam):
    try:
        cmd = ["samtools", "cat", "-o", out_bam, "-b", bam_list_fspath]
        p = sp.Popen(cmd, stdout=sp.PIPE)
        p.communicate()
    except Exception as e:
        print(e)


def _run_realigner(
    fqs: list[tuple], ref: Path, out_bam: Path, nproc: int = 4
) -> None:
    pass
