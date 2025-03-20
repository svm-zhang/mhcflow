import multiprocessing as mp
import subprocess as sp
import sys
from functools import partial
from pathlib import Path

from tinyscibio import BAMetadata, _PathLike, make_dir, parse_path

from .helper import FileManifest, _check_rg_exists, _check_single_rg, _get_sm
from .logger import logger


def _novoalign(
    task: tuple[_PathLike, _PathLike, _PathLike],
    fa: _PathLike,
    rg: dict[str, str],
) -> tuple[Path, Path, Path]:
    r1, r2, bam_out = task
    r1 = parse_path(r1)
    r2 = parse_path(r2)
    bam_out = parse_path(bam_out)
    realn_log = bam_out.with_suffix(".log")
    realn_done = bam_out.with_suffix(".done")
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
    cat_log = bam_out.with_suffix(".log")
    cat_done = bam_out.with_suffix(".done")
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
    sort_log = bam_out.with_suffix(".sort.log")
    sort_done = bam_out.with_suffix(".sort.done")
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


def _run_realigner(
    bam_fspath,
    ref: _PathLike,
    fisher_fm_json: _PathLike,
    outdir: _PathLike,
    nproc: int = 1,
) -> FileManifest:
    logger.info("Realign fished reads to HLA reference.")
    bametadata = BAMetadata(bam_fspath)
    _check_rg_exists(bametadata)
    _check_single_rg(bametadata)
    rg = bametadata.read_groups[0]
    sm = _get_sm(rg)

    outdir = parse_path(outdir)
    make_dir(outdir, parents=True, exist_ok=True)

    realigner_fm = FileManifest()
    realigner_fm_json = outdir / f"{sm}.realinger.file_manifest.json"
    if realigner_fm_json.exists():
        pass

    realigner_done = outdir / f"{sm}.realigner.done"

    fisher_fm_json = parse_path(fisher_fm_json)
    if not fisher_fm_json.exists():
        raise FileNotFoundError()

    # getting relevant outputs from fisher step.
    fisher_fm = FileManifest._from_json(fisher_fm_json)
    r1s = fisher_fm.outputs.get("r1s", [])
    r2s = fisher_fm.outputs.get("r2s", [])
    if not r1s or not r2s:
        raise FileNotFoundError()
    assert isinstance(r1s, list) and isinstance(r2s, list)
    if len(r1s) != len(r2s):
        # r1 and r2 must come in pairs
        raise ValueError()

    realn_tasks = []
    for i in range(len(r1s)):
        split_bam_out = outdir / f"{sm}.hla.realn.{i}.bam"
        r1, r2 = r1s[i], r2s[i]
        realn_tasks.append((r1, r2, split_bam_out))
    bams = []
    logs = []
    dones = []
    with mp.Pool(processes=nproc) as pool:
        for res in pool.imap_unordered(
            partial(_novoalign, fa=ref, rg=rg), realn_tasks
        ):
            bam_out, realn_log, realn_done = res
            logs.append(realn_log)
            dones.append(realn_done)
            bams.append(bam_out)

    concat_bam = outdir / f"{sm}.hla.realn.merged.bam"
    bam_list_fspath = outdir / "bams.list.txt"
    with open(bam_list_fspath, "w") as f:
        f.write("\n".join([str(bam) for bam in bams]))
    _, concat_log, concat_done = _concat(bam_list_fspath, concat_bam)

    realn_bam = outdir / f"{sm}.hla.realn.bam"
    _, sort_log, sort_done = _sort(
        bam_in=concat_bam, bam_out=realn_bam, nproc=nproc
    )

    logger.info(f"Realignment result in {str(realn_bam)}")
    realigner_done.touch()

    realigner_fm._register_inputs(fisher_json=fisher_fm_json, r1s=r1s, r2s=r2s)
    realigner_fm._register_outputs(realn_bam=realn_bam)
    realigner_fm._register_aux(done=realigner_done)
    realigner_fm._register_intermediate(
        bams=bams,
        concat_bam=concat_bam,
        concat_bam_list=bam_list_fspath,
    )
    realigner_fm._register_intermediate_aux(
        realn_dones=dones,
        realn_logs=logs,
        concat_log=concat_log,
        concat_done=concat_done,
        sort_log=sort_log,
        sort_done=sort_done,
    )
    realigner_fm._to_json(realigner_fm_json)
    return realigner_fm
