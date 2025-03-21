import subprocess as sp
from pathlib import Path

import polars as pl
from pyfaidx import Faidx
from tinyscibio import BAMetadata, _PathLike, make_dir, parse_path

from .helper import FileManifest, _check_rg_exists, _check_single_rg, _get_sm
from .logger import logger
from .realigner import _run_realigner


def dump_seq(allele: str, fa: Faidx, out: Path):
    with open(out, "a") as f:
        sequence = fa.fetch(allele, 1, fa.index[allele].rlen)
        f.write(f">{allele}\n{str(sequence)}\n")


def _run_finalizer(
    bam_fspath: _PathLike,
    ref: _PathLike,
    fisher_fm: FileManifest,
    realigner_fm: FileManifest,
    typer_res: pl.DataFrame,
    outdir: _PathLike,
    nproc: int = 1,
    overwrite: bool = False,
) -> FileManifest:
    bametadata = BAMetadata(str(bam_fspath))
    _check_rg_exists(bametadata)
    _check_single_rg(bametadata)
    rg = bametadata.read_groups[0]
    sm = _get_sm(rg)

    finalizer_fm = FileManifest()

    ref = parse_path(ref)
    fai = ref.parent / parse_path(f"{ref.name}.fai")
    if not fai.exists():
        raise FileNotFoundError

    outdir = parse_path(outdir)
    make_dir(outdir, parents=True, exist_ok=True)

    # 3 locus * 2 alleles = 6
    if typer_res.shape[0] != 6:
        raise ValueError

    fa_out = outdir / f"{sm}.hla.fasta"
    if fa_out.exists():
        fa_out.unlink()
    fa = Faidx(ref)
    # do not forget to take unique for homozygous genotype of HLA gene.
    # I dont want duplicated sequences in the fasta.
    _ = list(
        map(
            lambda x: dump_seq(x, fa, fa_out),
            typer_res["allele"].unique().to_list(),
        )
    )
    _novoindex(fa_out)

    # _run_realigner(fished_fqs, fa_out, outdir, rg, nproc)
    return finalizer_fm
