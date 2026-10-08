"""Render alignment strips and domain-architecture figures."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .interpro import DomainArchitecture
from .residues import RegionDiagnostic
from . import panels


def render_domain_cartoon(
    arches: list[tuple[str, DomainArchitecture, bool]],
    out_path: Path,
    title: str | None = None,
    right_label_family: str | None = None,
    use_full_pfams: bool = False,
) -> None:
    """Render one domain-architecture cartoon as its own figure.

    The drawing itself is `panels.draw_domain_cartoon`, which the manuscript
    figure calls with the same arguments; this wrapper owns only the page.
    See that function for what each argument selects.
    """
    fig, ax = plt.subplots(figsize=(11.5, 0.9 + 0.65 * len(arches)))
    panels.draw_domain_cartoon(
        ax, arches, title=title, right_label_family=right_label_family,
        use_full_pfams=use_full_pfams)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)


def render_msa_panel(
    aligned_fa_path: Path,
    diagnostic_columns: list[int],
    out_path: Path,
    title: str,
    start_col: int,
    end_col: int,
) -> None:
    """Render an MSA slice with diagnostic columns marked.

    Title is rendered as a free-standing figure-level text strip with
    extra top margin so it does not collide with the MsaViz tracks.
    """
    from pymsaviz import MsaViz
    # wrap_length larger than any diagnostic window: these regions do not
    # exceed ~80 columns, so the alignment fits one row.
    span = end_col - start_col + 1
    mv = MsaViz(
        str(aligned_fa_path),
        wrap_length=max(span, 90),
        show_consensus=True,
        show_count=False,
        color_scheme="Clustal",
        start=start_col + 1,                # pyMSAviz is 1-based
        end=end_col + 1,
    )
    in_window = [c for c in diagnostic_columns if start_col <= c < end_col]
    if in_window:
        try:
            mv.add_markers(positions=[c + 1 for c in in_window], color="red", marker="*")
        except Exception:
            pass
    fig = mv.plotfig()
    # ~1 inch of headroom above the MsaViz tracks, regardless of the
    # alignment's natural height, to hold the title strip clear of the boxes.
    w, h = fig.get_size_inches()
    headroom_in = 0.9
    new_h = h + headroom_in
    fig.set_size_inches(w, new_h)
    # Shift the existing axes downward so the top margin becomes the title strip.
    headroom_frac = headroom_in / new_h
    for ax in fig.axes:
        bbox = ax.get_position()
        ax.set_position([bbox.x0,
                         bbox.y0 * (1.0 - headroom_frac),
                         bbox.width,
                         bbox.height * (1.0 - headroom_frac)])
    fig.text(0.02, 1.0 - headroom_frac / 2, title,
             fontsize=12, ha="left", va="center", weight="bold")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
