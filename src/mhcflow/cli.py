import argparse

from tinyscibio import parse_path


def parse_cmd() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to BAM file",
    )
    parser.add_argument(
        "--ref",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to HLA reference in Fasta format.",
    )
    parser.add_argument(
        "--tag",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to HLA tag sequence file.",
    )
    parser.add_argument(
        "--bed",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to HLA region in BED format.",
    )
    parser.add_argument(
        "--freq",
        metavar="FILE",
        type=parse_path,
        required=True,
        help="specify path to HLA frequency file.",
    )
    parser.add_argument(
        "--outdir",
        metavar="DIR",
        type=parse_path,
        required=True,
        help="specify path to output folder.",
    )
    parser.add_argument(
        "--min-ecnt",
        metavar="INT",
        type=int,
        default=999,
        help="specify minimum # of mm events (999).",
    )
    parser.add_argument(
        "--nproc",
        metavar="INT",
        type=int,
        default=8,
        help="specify # processes to use (8).",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="specify to overwrite scores."
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="specify to NOT clean intermediate files.",
    )
    parser.add_argument(
        "--debug", action="store_true", help="specify to enter debug mode."
    )
    return parser
