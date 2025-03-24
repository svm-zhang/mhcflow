from pathlib import Path

import polars as pl
from pyfaidx import Faidx
from tinyscibio import BAMetadata, _PathLike, make_dir, parse_path

from .helper import (
    FileManifest,
    _check_rg_exists,
    _check_single_rg,
    _get_sm,
    _verify_prev_run,
)
from .logger import logger
from .realigner import _run_realigner
from .runnable import _novoindex


def _dump_seq(allele: str, fa: Faidx, out: Path) -> None:
    with open(out, "a") as f:
        sequence = fa.fetch(allele, 1, fa.index[allele].rlen)
        f.write(f">{allele}\n{str(sequence)}\n")


def _run_finalizer(
    bam_fspath: _PathLike,
    ref: _PathLike,
    fisher_fm_json: _PathLike,
    typer_res_fspath: _PathLike,
    outdir: _PathLike,
    nproc: int = 1,
    overwrite: bool = False,
) -> FileManifest:
    logger.info("Finalize sample HLA reference and realignemnt.")
    bametadata = BAMetadata(str(bam_fspath))
    _check_rg_exists(bametadata)
    _check_single_rg(bametadata)
    rg = bametadata.read_groups[0]
    sm = _get_sm(rg)

    outdir = parse_path(outdir)
    make_dir(outdir, parents=True, exist_ok=True)

    fm_json = outdir / f"{sm}.finalizer.file_manifest.json"
    if fm_json.exists():
        logger.info(
            f"Detected file manifest from previous run: {str(fm_json)}"
        )
        finalizer_fm = FileManifest._from_json(fm_json)
        if _verify_prev_run(finalizer_fm, overwrite):
            return finalizer_fm

    finalizer_fm = FileManifest()
    logdir = outdir / "log"
    make_dir(logdir, parents=True, exist_ok=True)
    finalizer_done = logdir / f"{sm}.finalizer.done"

    ref = parse_path(ref)
    fai = ref.parent / parse_path(f"{ref.name}.fai")
    if not fai.exists():
        raise FileNotFoundError

    typer_res = pl.read_csv(typer_res_fspath, separator="\t")
    # 3 locus * 2 alleles = 6
    if typer_res.shape[0] != 6:
        raise ValueError

    fa_out = outdir / f"{sm}.hla.fasta"
    fa = Faidx(ref)
    # do not forget to take unique for homozygous genotype of HLA gene.
    # I dont want duplicated sequences in the fasta.
    _ = list(
        map(
            lambda x: _dump_seq(x, fa, fa_out),
            typer_res["allele"].unique().to_list(),
        )
    )
    nix, index_log, index_done = _novoindex(fa_out)

    realigner_fm = _run_realigner(
        bam_fspath, nix, fisher_fm_json, outdir, nproc, overwrite
    )
    finalizer_done.touch()

    # record file manifest for finalizer
    finalizer_fm._register_inputs(hlaref=str(ref), typer_res=typer_res_fspath)
    finalizer_fm._register_aux(done=finalizer_done, myself=fm_json)
    finalizer_fm._register_outputs(
        sample_hlaref=str(fa_out),
        **realigner_fm.outputs,
    )
    finalizer_fm._register_intermediate(**realigner_fm.intermediates)
    finalizer_fm._register_intermediate_aux(
        index_log=index_log,
        index_done=index_done,
        **realigner_fm.intermediate_aux,
    )
    finalizer_fm._to_json(fm_json)

    return finalizer_fm
