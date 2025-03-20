import inspect
import json
import shlex
import subprocess as sp
import sys
from dataclasses import dataclass, field
from pathlib import Path

from tinyscibio import BAMetadata, _PathLike, parse_path

from .logger import logger


@dataclass
class FileManifest:
    _inputs: dict[str, _PathLike] = field(default_factory=dict, init=False)
    _outputs: dict[str, _PathLike | list[_PathLike]] = field(
        default_factory=dict, init=False
    )
    _aux: dict[str, _PathLike] = field(default_factory=dict, init=False)
    _intermediates: dict[str, _PathLike | list[_PathLike]] = field(
        default_factory=dict, init=False
    )

    @property
    def inputs(self) -> dict[str, _PathLike]:
        return self._inputs

    @property
    def outputs(self) -> dict[str, _PathLike | list[_PathLike]]:
        return self._outputs

    @property
    def aux(self) -> dict[str, _PathLike]:
        return self._aux

    @property
    def intermediates(self) -> dict[str, _PathLike | list[_PathLike]]:
        return self._intermediates

    def _register_inputs(self, **kwargs: _PathLike) -> None:
        self._inputs.update(**kwargs)

    def _register_outputs(self, **kwargs: _PathLike | list[_PathLike]) -> None:
        self._outputs.update(**kwargs)

    def _register_aux(self, **kwargs: _PathLike) -> None:
        self._aux.update(**kwargs)

    def _register_intermediate(
        self, **kwargs: _PathLike | list[_PathLike]
    ) -> None:
        self._intermediates.update(**kwargs)

    @classmethod
    def _from_json(cls, json_fspath: _PathLike) -> "FileManifest":
        with open(json_fspath, "r") as f:
            test = json.load(f)
        fm = cls()
        for k, v in test.items():
            setattr(fm, f"_{k}", v)
        return fm

    def _to_json(self, json_out: _PathLike) -> None:
        attrs = {}
        for attr_k, attr_v in inspect.getmembers(FileManifest):
            if not isinstance(attr_v, property):
                continue
            attr_v = getattr(self, attr_k)
            for item in attr_v.keys():
                if isinstance(attr_v[item], list):
                    attr_v[item] = [str(i) for i in attr_v[item]]
                else:
                    attr_v[item] = str(attr_v[item])
            attrs[attr_k] = attr_v
        with open(json_out, "w") as f:
            json.dump(attrs, f)

    # TODO: implement
    def _clean(self) -> None:
        raise NotImplementedError


def _extract_from_bam(
    idx_fspath: _PathLike, bam_fspath: _PathLike
) -> tuple[Path, Path]:
    logger.initialize()
    idx_fspath = parse_path(idx_fspath)
    r1 = idx_fspath.with_suffix(".R1.fastq")
    r2 = idx_fspath.with_suffix(".R2.fastq")

    if r1.exists() and r2.exists():
        logger.info(
            f"Found {r1} and {r2} extracted for read ids in {idx_fspath} file."
        )
        return (r1, r2)
    try:
        cmd_1 = f"samtools view -h -N {str(idx_fspath)} {str(bam_fspath)}"
        p1 = sp.Popen(shlex.split(cmd_1), stdout=sp.PIPE)
        cmd_2 = "samtools sort -n"
        p2 = sp.Popen(shlex.split(cmd_2), stdin=p1.stdout, stdout=sp.PIPE)
        cmd_3 = f"samtools fastq -n -1 {r1} -2 {r2} -0 /dev/null -s /dev/null"
        merged_cmd = " | ".join([cmd_1, cmd_2, cmd_3])
        logger.info(f"Extract reads into fastq using cmd: {merged_cmd}")
        p3 = sp.Popen(
            shlex.split(cmd_3),
            stdin=p2.stdout,
            stdout=sp.DEVNULL,
            stderr=sp.DEVNULL,
        )
        p3.communicate()
        p1.wait()
        p2.wait()
    except Exception as e:
        print(e)
    return (r1, r2)


def _check_rg_exists(bametadata: BAMetadata) -> None:
    try:
        rg = bametadata.read_groups
        if not rg:
            raise ValueError(
                "Failed to find any read group information in BAM: "
                f"{bametadata.fspath}"
            )
    except ValueError as e:
        logger.error(e)
        sys.exit(1)


def _check_single_rg(bametadata: BAMetadata) -> None:
    try:
        rg = bametadata.read_groups
        if len(rg) > 1:
            raise ValueError(
                f"Found more than one read group information in BAM: {rg}"
            )
    except ValueError as e:
        logger.error(e)
        sys.exit(1)


def _get_sm(rg: dict[str, str]) -> str:
    try:
        sm = rg.get("SM", None)
        if sm is None:
            raise ValueError(f"Failed to find SM field in read group: {rg}.")
        return sm
    except ValueError as e:
        logger.error(e)
        sys.exit(1)
