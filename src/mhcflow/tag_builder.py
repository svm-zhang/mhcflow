import pickle
import time
from typing import TypeAlias

import ahocorasick
from tinyscibio import _PathLike, parse_path

from .logger import logger

PREBUILD_METHOD = ["ahocorasick"]

PREBUILD_TYPE: TypeAlias = ahocorasick.Automaton


def _build_ahocorasick_automaton(
    tag_fspath: _PathLike, out: _PathLike
) -> PREBUILD_TYPE:
    logger.info("Start to build automaton using ahocorasick algorithm.")
    A = ahocorasick.Automaton(ahocorasick.STORE_LENGTH)
    start_t = time.time()
    with open(tag_fspath, "r") as fIN:
        for kmer in fIN.readlines():
            A.add_word(kmer.strip())
    A.make_automaton()
    logger.info(f"Automaton built in: {time.time() - start_t} sec")
    logger.info(f"Dump automaton to file: {out}")
    with open(out, "wb") as fOUT:
        pickle.dump(A, fOUT)
    return A


def _build_bloom_filter(tag_fspath: _PathLike, out: _PathLike) -> None:
    raise NotImplementedError(
        "Building bloom filter method has not been implemented."
    )


def _load_ahocorasick_automaton(fspath: _PathLike) -> PREBUILD_TYPE:
    with open(fspath, "rb") as fIN:
        return pickle.load(fIN)


def _load_bloom_filter() -> None:
    raise NotImplementedError(
        "Loading prebuilt bloom filter has not been implemented."
    )


def build(
    tag_seq_fspath: _PathLike, method: str = "ahocorasick"
) -> PREBUILD_TYPE:
    tag_seq_fspath = parse_path(tag_seq_fspath)
    if not tag_seq_fspath.exists():
        raise FileNotFoundError(
            f"Failed to find given path to tag sequence: {tag_seq_fspath}."
        )

    if method not in PREBUILD_METHOD:
        raise ValueError(
            f"Unsupported prebuild method specified: {method}. "
            f"Supported ones are: {PREBUILD_METHOD}."
        )
    prebuild_fspath = tag_seq_fspath.with_suffix(".prebuild")
    if prebuild_fspath.exists():
        return load(prebuild_fspath, method="ahocorasick")

    match method:
        case "ahocorasick":
            return _build_ahocorasick_automaton(
                tag_seq_fspath, prebuild_fspath
            )
        case "bloom":
            _build_bloom_filter(tag_seq_fspath, prebuild_fspath)
        case _:  # no pragma
            pass


def load(
    prebuild_fspath: _PathLike, method: str = "ahocorasick"
) -> PREBUILD_TYPE:
    match method:
        case "ahocorasick":
            return _load_ahocorasick_automaton(prebuild_fspath)
        case "bloom":
            return _load_bloom_filter()
        case _:  # no pragma
            pass
