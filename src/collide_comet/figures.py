"""Saving figures to the project's ``Figures`` folder."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, List, Optional, Union

__all__ = ["figure_dir", "set_figure_dir", "save_figure"]

_FIGURE_DIR: Optional[Path] = None
_ENV_VAR = "COLLIDE_COMET_FIGURES"
_KNOWN_FORMATS = {"png", "pdf", "svg", "jpg", "jpeg", "eps", "tif", "tiff", "ps", "webp"}


def _project_root(start: Optional[Path] = None) -> Optional[Path]:
    """Nearest parent of ``start`` that contains ``pyproject.toml`` or ``.git``."""
    start = Path.cwd() if start is None else Path(start)
    for folder in (start, *start.parents):
        if (folder / "pyproject.toml").exists() or (folder / ".git").exists():
            return folder
    return None


def figure_dir(subfolder: Optional[Union[str, Path]] = None, create: bool = True) -> Path:
    """Folder that figures are saved to.

    Resolution order: :func:`set_figure_dir`, the ``COLLIDE_COMET_FIGURES``
    environment variable, ``<project root>/Figures`` (the project root being
    the nearest parent of the working directory with a ``pyproject.toml`` or
    ``.git``), and finally ``./Figures``.  So notebooks kept in
    ``notebooks/`` write to the repository's top-level ``Figures`` folder.
    """
    if _FIGURE_DIR is not None:
        base = _FIGURE_DIR
    elif os.environ.get(_ENV_VAR):
        base = Path(os.environ[_ENV_VAR])
    else:
        root = _project_root()
        base = (root if root is not None else Path.cwd()) / "Figures"
    if subfolder:
        base = base / subfolder
    if create:
        base.mkdir(parents=True, exist_ok=True)
    return base


def set_figure_dir(path: Optional[Union[str, Path]]) -> Optional[Path]:
    """Set (or with ``None`` reset) the folder used by :func:`save_figure`."""
    global _FIGURE_DIR
    _FIGURE_DIR = None if path is None else Path(path).expanduser()
    return _FIGURE_DIR


def save_figure(
    fig,
    name: str,
    subfolder: Optional[Union[str, Path]] = None,
    formats: Union[str, Iterable[str]] = "png",
    dpi: int = 300,
    **savefig_kwargs,
) -> List[Path]:
    """Save ``fig`` as ``Figures/[subfolder/]name.<format>``.

    Parameters
    ----------
    name : str
        File name without extension (an extension, if given, selects the
        format).
    formats : str or iterable of str
        One or more of ``"png"``, ``"pdf"``, ``"svg"``, ``"jpg"`` ...
    dpi : int
        Resolution for raster formats.

    Returns the list of files written.
    """
    name = str(name)
    stem, ext = os.path.splitext(name)
    if ext.lower().lstrip(".") in _KNOWN_FORMATS:
        name, formats = stem, ext.lstrip(".")
    if isinstance(formats, str):
        formats = [formats]
    folder = figure_dir(subfolder)
    savefig_kwargs.setdefault("bbox_inches", "tight")
    written = []
    for fmt in formats:
        path = folder / f"{name}.{fmt.lstrip('.')}"
        fig.savefig(path, dpi=dpi, **savefig_kwargs)
        written.append(path)
    return written
