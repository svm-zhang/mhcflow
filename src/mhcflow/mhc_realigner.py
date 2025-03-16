import multiprocessing as mp
import subprocess as sp
from functools import partial
from pathlib import Path

from .logger import logger


def _novoalign(
    task: tuple[Path, Path, Path, Path, Path], nix: Path, rg: list
) -> Path:
    r1, r2, bam_out, realn_log, realn_done = task
    logger.initialize()
    if realn_done.exists():
        logger.info(f"Realignment for {r1.name}, {r2.name} has been done.")
        return bam_out
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
        cmd_2 = ["samtools", "view", "-bh", "-o", str(bam_out)]
        cmd_str = " | ".join([" ".join(cmd_1), " ".join(cmd_2)])
        logger.info(
            "Realign to HLA reference with fished reads: "
            f"{r1.name}, {r2.name}."
        )
        with open(realn_log, "w") as f:
            f.write(f"{cmd_str}\n")
            p1 = sp.Popen(cmd_1, stdout=sp.PIPE, stderr=f)
            p2 = sp.Popen(cmd_2, stdin=p1.stdout, stdout=sp.PIPE)
            p2.communicate()
            p1.wait()
        if not bam_out.exists():
            raise FileNotFoundError(
                f"Failed to find realigned BAM: {bam_out}."
                f"Realignment failed for read pair {r1}, {r2}."
            )
        realn_done.touch()
    except Exception as e:
        print(e)
        raise SystemExit
    return bam_out


def _concat(bam_list_fspath: Path, bam_out: Path) -> None:
    logger.info("Concatenate individual bam files.")
    cat_log = bam_out.with_suffix(".log")
    cat_done = bam_out.with_suffix(".done")
    if cat_done.exists():
        logger.info(
            "Found concatenated realigned BAM file from "
            f"previous run: {bam_out}."
        )
        return
    try:
        cmd = [
            "samtools",
            "cat",
            "-o",
            str(bam_out),
            "-b",
            str(bam_list_fspath),
        ]
        with open(cat_log, "w") as f:
            f.write(" ".join(cmd) + "\n")
            p = sp.Popen(cmd, stdout=f, stderr=sp.STDOUT)
            p.communicate()
        cat_done.touch()
    except Exception as e:
        print(e)
        raise SystemExit


def _sort(bam_in: Path, bam_out: Path, nproc: int = 1) -> None:
    logger.info(f"Sort concatenated BAM file: {bam_in.name}.")
    bai = bam_out.with_suffix(".bam.bai")
    sort_log = bam_out.with_suffix(".sort.log")
    sort_done = bam_out.with_suffix(".sort.done")
    if sort_done.exists():
        logger.info(
            f"Found sorted BAM result from previous run: {str(bam_out)}. Skip."
        )
        return
    try:
        cmd = [
            "samtools",
            "sort",
            "-@",
            f"{nproc}",
            "--write-index",
            "-o",
            f"{str(bam_out)}##idx##{str(bai)}",
            str(bam_in),
        ]
        with open(str(sort_log), "w") as f:
            f.write(" ".join(cmd) + "\n")
            p = sp.Popen(cmd, stdout=f, stderr=sp.STDOUT)
            p.communicate()
        sort_done.touch()
    except Exception as e:
        print(e)
        raise SystemExit


# TODO: add clean param
def _run_realigner(
    fished_fqs: list[tuple[Path, Path]],
    ref: Path,
    outdir: Path,
    rg,
    nproc: int = 1,
) -> Path:
    logger.info("Realign fished reads to HLA reference.")
    # this should recover the read group info in the original BAM
    sm = rg.get("SM")
    realn_bam = outdir / f"{sm}.hla.realn.bam"
    if realn_bam.exists():
        logger.info(f"Found previously realigned BAM file: {realn_bam}. Skip.")
        return realn_bam
    rg = [f"{k}:{v}" for k, v in rg.items()]
    realn_tasks = []
    for i in range(len(fished_fqs)):
        split_bam_out = outdir / f"hla.realn.split.{i}.bam"
        split_log = outdir / f"hla.realn.split.{i}.log"
        split_done = outdir / f"hla.realn.split.{i}.done"
        r1, r2 = fished_fqs[i]
        realn_tasks.append((r1, r2, split_bam_out, split_log, split_done))
    bams = []
    with mp.Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_novoalign, nix=ref, rg=rg), realn_tasks
        ):
            bams.append(res)
    concat_bam = outdir / "hla.realn.merged.bam"
    bam_list_fspath = outdir / "bams.list.txt"
    with open(bam_list_fspath, "w") as f:
        f.write("\n".join([str(bam) for bam in bams]))
    _concat(bam_list_fspath, concat_bam)

    _sort(bam_in=concat_bam, bam_out=realn_bam, nproc=nproc)
    logger.info(f"Realignment result in {str(realn_bam)}")
    return realn_bam
