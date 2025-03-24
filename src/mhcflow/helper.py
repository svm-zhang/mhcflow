import inspect
import json
import sys
from collections.abc import Sequence
from dataclasses import dataclass, field

from tinyscibio import BAMetadata, _PathLike, parse_path

from .logger import logger


@dataclass
class FileManifest:
    _inputs: dict[str, _PathLike | Sequence[_PathLike]] = field(
        default_factory=dict, init=False
    )
    _outputs: dict[str, _PathLike | Sequence[_PathLike]] = field(
        default_factory=dict, init=False
    )
    _aux: dict[str, _PathLike] = field(default_factory=dict, init=False)
    _intermediates: dict[str, _PathLike | Sequence[_PathLike]] = field(
        default_factory=dict, init=False
    )
    _intermediate_aux: dict[str, _PathLike | Sequence[_PathLike]] = field(
        default_factory=dict, init=False
    )

    @property
    def inputs(self) -> dict[str, _PathLike | Sequence[_PathLike]]:
        return self._inputs

    @property
    def outputs(self) -> dict[str, _PathLike | Sequence[_PathLike]]:
        return self._outputs

    @property
    def aux(self) -> dict[str, _PathLike]:
        return self._aux

    @property
    def intermediate_aux(self) -> dict[str, _PathLike | Sequence[_PathLike]]:
        return self._intermediate_aux

    @property
    def intermediates(self) -> dict[str, _PathLike | Sequence[_PathLike]]:
        return self._intermediates

    def _register_inputs(
        self, **kwargs: _PathLike | Sequence[_PathLike]
    ) -> None:
        self._inputs.update(**kwargs)

    def _register_outputs(
        self, **kwargs: _PathLike | Sequence[_PathLike]
    ) -> None:
        self._outputs.update(**kwargs)

    def _register_aux(self, **kwargs: _PathLike) -> None:
        self._aux.update(**kwargs)

    def _register_intermediate(
        self, **kwargs: _PathLike | Sequence[_PathLike]
    ) -> None:
        self._intermediates.update(**kwargs)

    def _register_intermediate_aux(
        self, **kwargs: _PathLike | Sequence[_PathLike]
    ) -> None:
        self._intermediate_aux.update(**kwargs)

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

    def _reset(self) -> None:
        for attr, attr_v in inspect.getmembers(FileManifest):
            if not isinstance(attr_v, property):
                continue
            if attr == "inputs":
                continue
            self._clean_attr(attr)

    def _clean_attrs(self, attrs: list[str]) -> None:
        for attr in attrs:
            self._clean_attr(attr)

    def _clean_attr(self, attr: str) -> None:
        to_be_cleaned = getattr(self, attr)
        if not to_be_cleaned:
            return
        for k, v in to_be_cleaned.items():
            if isinstance(v, list):
                for f in v:
                    parse_path(f).unlink(missing_ok=True)
                continue
            parse_path(v).unlink(missing_ok=True)


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


def _verify_prev_run(fm: FileManifest, overwrite: bool = False) -> bool:
    # if overwrite specify, destroy all files generated previously
    # except everything in fm.inputs.
    if overwrite:
        logger.info("Overwrite specified. Remove results from previous run.")
        fm._reset()
        return False

    logger.info("Verify manifest from previous run.")
    done = parse_path(fm.aux.get("done", ""))
    intermediate_dones = []
    for k, v in fm.intermediate_aux.items():
        if not k.endswith("done") and not k.endswith("dones"):
            continue
        intermediate_dones += v if isinstance(v, list) else [v]
    not_exists = [
        str(f) for f in intermediate_dones if not parse_path(f).exists()
    ]
    not_exists += [str(done)] if not done.exists() else []
    if not not_exists:
        logger.info("Verified all done files from previous run. Skip.")
        return True
    logger.info(
        "Failed to verify previous run due to missing some done files: "
        f"{not_exists}"
    )
    logger.info(
        "Clean intermediate outputs and done files for a fresh new run."
    )
    fm._reset()
    return False
