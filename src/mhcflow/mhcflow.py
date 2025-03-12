import polars as pl
from mhctyper import run_mhctyper
from tinyscibio import BAMetadata, make_dir

from .cli import parse_cmd
from .fastq import dump_to_fastq
from .logger import logger
from .strawlr import run_strawlr


def run_mhcflow():
    parser = parse_cmd()
    args = parser.parse_args()

    logger.initialize()
    make_dir(args.outdir, parents=True, exist_ok=True)

    bametadata = BAMetadata(args.bam)
    rg = bametadata.read_groups
    if not rg:
        raise ValueError("Failed to find any read group information in BAM.")
    if len(rg) > 1:
        raise ValueError("Found more than one read group information in BAM.")
    sm = bametadata.read_groups[0].get("SM", None)
    if sm is None:
        raise ValueError(f"Failed to find SM field in read group {rg}.")
    logger.info(sm)

    out_fisher_dir = args.outdir / "fisher"
    make_dir(out_fisher_dir, parents=True, exist_ok=True)
    out_fished_qnames = out_fisher_dir / f"{sm}.fished.ids.txt"
    if not out_fished_qnames.exists():
        # run strawlr to get read ids
        fished_qnames = run_strawlr(bametadata, args.tag, args.bed)
        fished_qnames.to_frame().write_csv(out_fished_qnames, separator="\t")
    else:
        logger.info(
            "Found results of fished reads from previous run: "
            f"{out_fished_qnames}"
        )
        fished_qnames = pl.read_csv(
            out_fished_qnames, separator="\t"
        ).to_series()

    print(fished_qnames.head())
    print(fished_qnames.shape)

    # get read sequence

    # realign them using novoalign

    # call mhctyper
    # FIXME: update mhctyper to also make it library
    # run_mhctyper()

    # get result and finalize
    pass
