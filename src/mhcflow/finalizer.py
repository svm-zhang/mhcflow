import subprocess as sp
from pathlib import Path

import polars as pl
from pyfaidx import Faidx
from tinyscibio import make_dir, parse_path

from .logger import logger
from .realigner import _run_realigner


def _novoindex(fa: Path) -> None:
    nix = fa.with_suffix(".nix")
    index_log = nix.with_suffix(".novoindex.log")
    index_done = nix.with_suffix(".novoindex.done")
    if index_done.exists():
        logger.info(f"Found .nix index for Fasta file: {fa}. Skip.")
        return
    try:
        cmd = ["novoindex", str(nix), str(fa)]
        with open(index_log, "w") as f:
            f.write(f"{' '.join(cmd)}\n")
            p = sp.Popen(cmd, stdout=f, stderr=sp.STDOUT)
            p.communicate()
        index_done.touch()
    except Exception as e:
        print(e)
        raise SystemExit
    return


def dump_seq(allele: str, fa: Faidx, out: Path):
    with open(out, "a") as f:
        sequence = fa.fetch(allele, 1, fa.index[allele].rlen)
        f.write(f">{allele}\n{str(sequence)}\n")


def finalize(
    typer_res: pl.DataFrame,
    ref: Path,
    fished_fqs: list[tuple[Path, Path]],
    rg: dict[str, str],
    outdir: Path,
    nproc: int = 1,
):
    fai = ref.parent / parse_path(f"{ref.name}.fai")
    if not fai.exists():
        raise FileNotFoundError

    make_dir(outdir, parents=True, exist_ok=True)

    # 3 locus * 2 alleles = 6
    if typer_res.shape[0] != 6:
        raise ValueError

    sample = rg.get("SM")
    fa_out = outdir / f"{sample}.hla.fasta"
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

    _run_realigner(fished_fqs, fa_out, outdir, rg, nproc)
