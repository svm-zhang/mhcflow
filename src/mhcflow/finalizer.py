from pathlib import Path

import polars as pl
from pyfaidx import Faidx
from tinyscibio import make_dir, parse_path

from .mhc_realigner import _run_realigner


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
    _ = list(
        map(lambda x: dump_seq(x, fa, fa_out), typer_res["allele"].to_list())
    )

    _run_realigner(fished_fqs, ref, outdir, rg, nproc)
