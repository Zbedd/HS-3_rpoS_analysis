"""House figure system — one style and one layout engine for every figure.

    import viz as fs

    fs.use()                                  # install the rcParams, once
    fig, ax = plt.subplots(figsize=(fs.COL_DOUBLE, 3.0))
    fs.finish_axis(ax)
    fs.save(fig, path)

`style` owns how a figure looks — colours, type, marks, chrome, export.
`layout` owns how panels are arranged. Both are re-exported here, so a
plotting module imports this package and nothing else from `viz`.
"""

from __future__ import annotations

from .style import *          # noqa: F401,F403
from .layout import *         # noqa: F401,F403
from . import layout, style   # noqa: F401

__all__ = list(style.__all__) + list(layout.__all__)
