from mhctyper import run_mhctyper
from tinyscibio import make_dir, parse_path

from .cli import parse_cmd
from .finalizer import _run_finalizer
from .fisher import _run_fisher
from .logger import logger
from .realigner import _run_realigner


def run_mhcflow() -> int:
    parser = parse_cmd()
    args = parser.parse_args()
    logger.initialize()
    logger.info("Start running mhcflow.")

    make_dir(args.outdir, parents=True, exist_ok=True)

    ref = args.ref.with_suffix(".nix")
    if not ref.exists():
        raise FileNotFoundError(
            f"Failed to find HLA reference novoalign index file: {ref}"
        )

    out_fisher_dir = args.outdir / "fisher"
    fisher_fm = _run_fisher(
        args.bam,
        args.tag,
        args.bed,
        out_fisher_dir,
        nproc=args.nproc,
        overwrite=args.overwrite,
    )
    fisher_fm_json = fisher_fm.aux.get("myself", "")

    out_realn_dir = args.outdir / "realigner"
    realigner_fm = _run_realigner(
        args.bam, args.ref, fisher_fm_json, out_realn_dir, args.nproc
    )

    realn_bam = realigner_fm.outputs.get("realn_bam", "")
    assert isinstance(realn_bam, str)
    typer_dir = args.outdir / "typer"
    _, typer_res_fspath = run_mhctyper(
        bam=parse_path(realn_bam),
        freq=args.freq,
        outdir=typer_dir,
        min_ecnt=args.min_ecnt,
        nproc=args.nproc,
        debug=args.debug,
        overwrite=args.overwrite,
    )

    finalizer_dir = args.outdir / "finalizer"
    finalizer_fm = _run_finalizer(
        args.bam,
        args.ref,
        fisher_fm_json,
        typer_res_fspath,
        finalizer_dir,
        nproc=args.nproc,
        overwrite=args.overwrite,
    )

    if not args.no_clean:
        logger.info("Clean intermediate files.")
        finalizer_fm._clean_attr("intermediates")
        realigner_fm._clean_attr("intermediates")
        fisher_fm._clean_attr("intermediates")
    logger.info("Finished running mhcflow.")
    return 0
