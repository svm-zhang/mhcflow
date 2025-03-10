import time
from functools import partial
from multiprocessing import get_context

import polars as pl
import pysam
from tinyscibio import _PathLike, parse_path, walk_bam

from .logger import logger
from .tag_builder import build


def _get_qnames_from_seq_with_tag(bam_fspath, prebuilt_tag) -> set[str]:
    start_t = time.time()
    qnames = set()
    with pysam.AlignmentFile(str(bam_fspath), "rb") as bamf:
        for aln in bamf.fetch("6"):
            assert aln.query_sequence is not None
            match = list(prebuilt_tag.iter(aln.query_sequence))
            if match:
                qnames.add(aln.query_name)
    print(f"{time.time() - start_t} sec")
    return qnames


def _get_qnames_from_unplace_seq_with_tag(
    bam_fspath, prebuilt_tag
) -> set[str]:
    qnames = set()
    start_t = time.time()
    with pysam.AlignmentFile(bam_fspath, "rb") as bamf:
        for aln in bamf.fetch(until_eof=True):
            if not aln.is_unmapped:
                continue
            assert aln.query_sequence is not None
            match = list(prebuilt_tag.iter(aln.query_sequence))
            if match:
                qnames.add(aln.query_name)
    print(f"{time.time() - start_t} sec")
    print(len(qnames))
    return qnames


def _get_qnames_from_seq_mapped_to_class1(bam_fspath) -> set[str]:
    regions = [
        ("6", 29909037, 29913661),
        ("6", 31236526, 31239869),
        ("6", 31321649, 31324964),
    ]
    qnames = set()
    start_t = time.time()
    with pysam.AlignmentFile(bam_fspath, "rb") as bamf:
        for region in regions:
            contig, start, end = region
            for aln in bamf.fetch(contig, start, end):
                qnames.add(aln.query_name)
    print(f"{time.time() - start_t} sec")
    print(len(qnames))
    return qnames


def run_strawlr(
    bam_fspath: _PathLike,
    tag_seq_fspath: _PathLike,
    prebuild_method: str = "ahocorasick",
):
    bam_fspath = parse_path(bam_fspath)

    prebuilt_tag = build(tag_seq_fspath)
    # print(type(prebuilt_tag))
    # print("TGCTGGAGTGTCCCAAGAGAGATGCAAAGTGTCTGAAT" in prebuilt_tag)

    qnames_set_1 = _get_qnames_from_seq_with_tag(bam_fspath, prebuilt_tag)
    print(len(qnames_set_1))
    qnames_set_2 = _get_qnames_from_unplace_seq_with_tag(
        bam_fspath, prebuilt_tag
    )
    print(len(qnames_set_2))
    qnames_set_3 = _get_qnames_from_seq_mapped_to_class1(bam_fspath)
    print(len(qnames_set_3))
