"""Colour maps used by the MATLAB figures."""

from __future__ import annotations

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

__all__ = ["redblue", "custom_rainbow", "get_cmap"]


def redblue(n: int = 256) -> ListedColormap:
    """Blue - white - red diverging map (MATLAB File Exchange ``redblue``).

    This is the colour map ``asymPlots.makePlot`` applies to every slice plot.
    """
    half = n // 2
    ramp = np.arange(half) / max(half - 1, 1)
    r = np.r_[ramp, np.ones(n - half)]
    g = np.r_[ramp, ramp[::-1]] if n % 2 == 0 else np.r_[ramp, 1.0, ramp[::-1]]
    b = np.r_[np.ones(n - half), ramp[::-1]]
    r, b = r[:n], b[-n:]
    return ListedColormap(np.column_stack([r, g[:n], b]), name="redblue")


def custom_rainbow(n: int = 256) -> LinearSegmentedColormap:
    """White-navy-cyan-green-yellow-orange-maroon-black map.

    The ``customMap`` defined (but not applied) in ``asymPlots.makePlot``.
    """
    anchors = np.array(
        [
            [255, 255, 255], [2, 5, 95], [2, 241, 229], [23, 142, 5], [251, 246, 13],
            [197, 62, 6], [105, 11, 13], [212, 72, 73], [0, 0, 0],
        ]
    ) / 255.0
    return LinearSegmentedColormap.from_list("collide_custom", anchors, N=n)


def get_cmap(cmap=None):
    """Resolve a colour map: ``None``/"redblue", "custom", a matplotlib name or object."""
    import matplotlib

    if cmap is None or (isinstance(cmap, str) and cmap.lower() == "redblue"):
        return redblue()
    if isinstance(cmap, str) and cmap.lower() in ("custom", "collide_custom"):
        return custom_rainbow()
    if isinstance(cmap, str):
        return matplotlib.colormaps[cmap]
    return cmap
