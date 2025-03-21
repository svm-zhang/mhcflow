import multiprocessing as mp
import shutil
from functools import partial

from tinyscibio import BAMetadata, _PathLike, make_dir, parse_path

from .helper import FileManifest, _check_rg_exists, _check_single_rg, _get_sm
from .logger import logger
from .runnable import _concat, _novoalign, _sort


def _run_realigner(
    bam_fspath,
    ref: _PathLike,
    fisher_fm_json: _PathLike,
    outdir: _PathLike,
    nproc: int = 1,
    overwrite: bool = False,
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
    realigner_fm_json = outdir / f"{sm}.realigner.file_manifest.json"
    # check json and skip entirely if all recorded done files exists
    # otherwise; keep going
    if realigner_fm_json.exists():
        if not overwrite:
            realigner_fm = FileManifest._from_json(realigner_fm_json)
            realigner_done = parse_path(realigner_fm.aux.get("done", ""))
            intermediate_dones = []
            for k, v in realigner_fm.intermediate_aux.items():
                if not k.endswith("done") and not k.endswith("dones"):
                    continue
                intermediate_dones += v if isinstance(v, list) else [v]
            n_not_exists = [
                f for f in intermediate_dones if not parse_path(f).exists()
            ]
            if realigner_done.exists() and not n_not_exists:
                logger.info(
                    "Found all done files for realigner from previous run. "
                    "Skip."
                )
                return realigner_fm
            logger.info(
                "Failed to skip realigner entirely."
                "Missing either realigner done or intermediate done files."
                f"Please check: {realigner_done}, {n_not_exists}"
            )
        else:
            logger.info(
                "Overwrite specified. "
                f"Remove realigner results from previous run: {outdir}"
            )
            # TODO: need a _clean method.
            shutil.rmtree(outdir)
            make_dir(outdir, parents=True, exist_ok=True)

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
    realigner_fm._register_aux(done=realigner_done, myself=realigner_fm_json)
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
