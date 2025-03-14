import multiprocessing as mp
import subprocess as sp
from functools import partial
from pathlib import Path

from .logger import logger


def _novoalign(
    task: tuple[Path, Path, Path, Path, Path], nix: Path, rg: list
) -> Path:
    r1, r2, out_realn_bam, out_realn_log, out_realn_done = task
    logger.initialize()
    if out_realn_done.exists():
        logger.info(f"Realignment for {r1.name}, {r2.name} has been done.")
        return out_realn_bam
    rg_str = "@RG\t" + "\t".join(rg)
    try:
        cmd_1 = [
            "novoalign",
            "-d",
            str(nix),
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
            rg_str,
            "-f",
            str(r1),
            str(r2),
        ]
        cmd_2 = ["samtools", "view", "-bh", "-o", str(out_realn_bam)]
        cmd_str = " | ".join([" ".join(cmd_1), " ".join(cmd_2)])
        logger.info(
            f"Realign fished reads to HLA reference using cmd: {cmd_str}."
        )
        with open(out_realn_log, "w") as f:
            p1 = sp.Popen(cmd_1, stdout=sp.PIPE, stderr=f)
            p2 = sp.Popen(cmd_2, stdin=p1.stdout, stdout=sp.PIPE)
            p2.communicate()
            p1.wait()
        if not out_realn_bam.exists():
            raise FileNotFoundError(
                f"Failed to find realigned BAM: {out_realn_bam}."
                f"Realignment failed for read pair {r1}, {r2}."
            )
        out_realn_done.touch()
    except Exception as e:
        print(e)
        raise SystemExit
    return out_realn_bam


def _concat(bam_list_fspath: Path, out_bam: Path) -> None:
    logger.info("Concatenate individual bam files.")
    out_cat_done = out_bam.with_suffix(".done")
    if out_cat_done.exists():
        logger.info(
            "Found concatenated realigned BAM file from "
            f"previous run: {out_bam}."
        )
        return
    try:
        cmd = [
            "samtools",
            "cat",
            "-o",
            str(out_bam),
            "-b",
            str(bam_list_fspath),
        ]
        logger.info(" ".join(cmd))
        p = sp.Popen(cmd, stdout=sp.PIPE)
        p.communicate()
        out_cat_done.touch()
    except Exception as e:
        print(e)


def _run_realigner(
    fished_fqs: list[tuple[Path, Path]],
    ref: Path,
    outdir: Path,
    rg,
    nproc: int = 1,
) -> Path:
    # this should recover the read group info in the original BAM
    sm = rg.get("SM")
    out_bam = outdir / f"{sm}.hla.realn.bam"
    if out_bam.exists():
        logger.info(f"Found previously realigned BAM file: {out_bam}. Skip.")
        return out_bam
    rg = [f"{k}:{v}" for k, v in rg.items()]
    realn_tasks = []
    for i in range(len(fished_fqs)):
        out_realn_bam = outdir / f"hla.realn.split.{i}.bam"
        out_realn_log = outdir / f"hla.realn.split.{i}.log"
        out_realn_done = outdir / f"hla.realn.split.{i}.done"
        r1, r2 = fished_fqs[i]
        realn_tasks.append(
            (r1, r2, out_realn_bam, out_realn_log, out_realn_done)
        )
    bams = []
    with mp.Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_novoalign, nix=ref, rg=rg), realn_tasks
        ):
            bams.append(res)
    bam_list_fspath = outdir / "bams.list.txt"
    with open(bam_list_fspath, "w") as f:
        f.write("\n".join([str(bam) for bam in bams]))
    _concat(bam_list_fspath, out_bam)
    return out_bam
