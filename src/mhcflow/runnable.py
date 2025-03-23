import shlex
import subprocess as sp
import sys

from tinyscibio import _PathLike, make_dir, parse_path

from .logger import logger


def _samtools_fastq(
    idx_fspath: _PathLike, bam_fspath: _PathLike
) -> tuple[_PathLike, _PathLike, _PathLike, _PathLike]:
    idx_fspath = parse_path(idx_fspath)
    logdir = idx_fspath.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    idx_done = logdir / f"{idx_fspath.stem}.done"
    idx_log = logdir / f"{idx_fspath.stem}.log"
    r1 = idx_fspath.with_suffix(".R1.fastq")
    r2 = idx_fspath.with_suffix(".R2.fastq")
    logger.initialize()
    if idx_done.exists():
        logger.info(
            "Found extraction done file from previous run: "
            f"{idx_fspath.name}. Skip."
        )
        return (r1, r2, idx_done, idx_log)

    try:
        logger.info(
            f"Extract reads given {idx_fspath.name} to {r1.name} {r2.name}"
        )
        with open(idx_log, "a") as f:
            cmd_1 = f"samtools view -h -N {str(idx_fspath)} {str(bam_fspath)}"
            cmd_2 = "samtools sort -n"
            cmd_3 = (
                f"samtools fastq -n -1 {r1} -2 {r2} -0 /dev/null -s /dev/null"
            )
            cmd_str = " | ".join([cmd_1, cmd_2, cmd_3])
            f.write(f"{cmd_str}\n")
            p1 = sp.Popen(shlex.split(cmd_1), stdout=sp.PIPE, stderr=f)
            p2 = sp.Popen(
                shlex.split(cmd_2), stdin=p1.stdout, stdout=sp.PIPE, stderr=f
            )
            p3 = sp.Popen(
                shlex.split(cmd_3),
                stdin=p2.stdout,
                stdout=f,
                stderr=sp.STDOUT,
            )
            p3.communicate()
            p1.wait()
            p2.wait()

        idx_done.touch()
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    return (r1, r2, idx_done, idx_log)


def _novoalign(
    task: tuple[_PathLike, _PathLike, _PathLike],
    fa: _PathLike,
    rg: dict[str, str],
) -> tuple[_PathLike, _PathLike, _PathLike]:
    r1, r2, bam_out = task
    r1 = parse_path(r1)
    r2 = parse_path(r2)
    bam_out = parse_path(bam_out)
    logdir = bam_out.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    realn_done = logdir / f"{bam_out.stem}.done"
    realn_log = logdir / f"{bam_out.stem}.log"
    logger.initialize()
    if realn_done.exists():
        logger.info(f"Realignment for {r1.name}, {r2.name} has been done.")
        return (bam_out, realn_log, realn_done)
    rg_lst = [f"{k}:{v}" for k, v in rg.items()]
    rg_str = "@RG\t" + "\t".join(rg_lst)
    nix = parse_path(fa).with_suffix(".nix")
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
        return (bam_out, realn_log, realn_done)
    except Exception as e:
        logger.error(e)
        sys.exit(1)


def _concat(
    bam_list_fspath: _PathLike, bam_out: _PathLike
) -> tuple[_PathLike, _PathLike, _PathLike]:
    logger.info("Concatenate individual bam files.")
    bam_out = parse_path(bam_out)
    logdir = bam_out.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    cat_done = logdir / f"{bam_out.stem}.done"
    cat_log = logdir / f"{bam_out.stem}.log"
    if cat_done.exists():
        logger.info(
            "Found concatenated realigned BAM file from "
            f"previous run: {bam_out}."
        )
        return bam_out, cat_log, cat_done
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
        return (bam_out, cat_log, cat_done)
    except Exception as e:
        logger.error(e)
        sys.exit(1)


def _sort(
    bam_in: _PathLike, bam_out: _PathLike, nproc: int = 1
) -> tuple[_PathLike, _PathLike, _PathLike]:
    bam_in = parse_path(bam_in)
    bam_out = parse_path(bam_out)
    bai = bam_out.with_suffix(".bam.bai")
    logdir = bam_out.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    sort_done = logdir / f"{bam_out.stem}.done"
    sort_log = logdir / f"{bam_out.stem}.log"
    logger.info(f"Sort concatenated BAM file: {bam_in.name}.")
    if sort_done.exists():
        logger.info(
            f"Found sorted BAM result from previous run: {str(bam_out)}. Skip."
        )
        return (bam_out, sort_log, sort_done)
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
        return (bam_out, sort_log, sort_done)
    except Exception as e:
        logger.error(e)
        sys.exit(1)


def _novoindex(fa: _PathLike) -> tuple[_PathLike, _PathLike, _PathLike]:
    logger.info(f"Index Fasta using novoindex: {fa}")
    nix = parse_path(fa).with_suffix(".nix")
    logdir = nix.parent / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    index_done = logdir / f"{nix.stem}.novoindex.done"
    index_log = logdir / f"{nix.stem}.novoindex.log"
    if index_done.exists():
        logger.info(f"Found .nix index for Fasta file: {fa}. Skip.")
        return (nix, index_log, index_done)
    try:
        cmd = ["novoindex", str(nix), str(fa)]
        with open(index_log, "w") as f:
            f.write(f"{' '.join(cmd)}\n")
            p = sp.Popen(cmd, stdout=f, stderr=sp.STDOUT)
            p.communicate()
        index_done.touch()
        return (nix, index_log, index_done)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
