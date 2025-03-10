from __future__ import annotations

import logging
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from typing import Optional

    from tinyscibio import _PathLike


class Logger(logging.Logger):
    def __init__(self, name: str):
        super().__init__(name)
        self._configured = False

    def initialize(
        self, debug: bool = False, f: Optional[_PathLike] = None
    ) -> None:
        if self._configured:
            return

        self.setLevel(logging.DEBUG if debug else logging.INFO)
        if self.hasHandlers():
            self.handlers.clear()

        console_handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "[%(module)s.py] - [%(funcName)s:%(lineno)d] - [%(levelname)s]: %(message)s"
        )
        console_handler.setFormatter(formatter)
        console_handler.setLevel(logging.INFO)
        self.addHandler(console_handler)

        if debug:
            if f is None:
                raise ValueError(
                    f"Parameter {f=} is required to add a FileHandler to "
                    "logger when debug is true."
                )
            file_handler = logging.FileHandler(f, mode="a")
            file_handler.setFormatter(formatter)
            file_handler.setLevel(logging.DEBUG)
            self.addHandler(file_handler)

        self._configured = True


logging.setLoggerClass(Logger)
# logging.getLogger returns Logger type
# Although setLoggerClass(Logger) sets the logger object as Logger class,
# the getLogger method, according to typeshed, still returns Logger type
# Thats why LSP (e.g. pyright) complains about initialize is not an attribute
# https://github.com/microsoft/pyright/discussions/5418
logger = cast(Logger, logging.getLogger(__name__))
