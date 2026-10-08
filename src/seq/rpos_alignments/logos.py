"""Exploratory RpoS/RpoD sequence logos and shared logo helpers."""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import logomaker
import pandas as pd
import numpy as np

from .residues import RegionDiagnostic
from .labels import Label


AA20 = "ACDEFGHIKLMNPQRSTVWY"
GAP = "-"


# Chemistry colour scheme used by logomaker — legend reference for both
CHEMISTRY_LEGEND = [
    ("polar (G S T Y C Q N)", "#109648"),
    ("basic (K R H)",          "#3060A8"),
    ("acidic (D E)",           "#CC3232"),
    ("hydrophobic (A V L I M P W F)", "#000000"),
]


# --------------------------------------------------------------- frequencies

def _freq_matrix(
    aligned: dict[str, str], member_ids: list[str], cols: list[int],
) -> pd.DataFrame:
    """Per-column amino-acid frequency matrix for a panel subset.

    Index = column number (0-based alignment col). Columns = each of the 20 AAs.
    """
    rows = []
    for col in cols:
        residues = [aligned[i][col] for i in member_ids if aligned[i][col] != GAP]
        c = Counter(residues)
        total = sum(c.values()) or 1
        rows.append({aa: c.get(aa, 0) / total for aa in AA20})
    df = pd.DataFrame(rows, index=cols, columns=list(AA20))
    return df


def _bits_matrix(freq_df: pd.DataFrame) -> pd.DataFrame:
    """Information-content (bits) matrix for the logo height per column.

    Uses Schneider's small-sample correction approximated as 0 (the panel is
    too small for an exact correction; we use plain information content).
    """
    bits_per_col: list[pd.Series] = []
    for _, row in freq_df.iterrows():
        h = -sum(p * np.log2(p) for p in row if p > 0)
        info = max(0.0, np.log2(20) - h)   # log2(20) ≈ 4.32
        bits_per_col.append(row * info)
    out = pd.DataFrame(bits_per_col)
    out.index = freq_df.index
    return out


# --------------------------------------------------------------- sigma logos

def render_sigma_diagnostic_logos(
    aligned: dict[str, str],
    region: RegionDiagnostic,
    rpos_ref_ids: list[str],
    rpod_ref_ids: list[str],
    hs3_id: str,
    out_path: Path,
    region_title: str,
    hs3_track_label: str = "HS-3 rpoS",
) -> None:
    """One supplementary-quality figure with stacked RpoS-vs-RpoD logos.

    Layout:
        Row 1 — RpoS reference consensus logo (information content, bits)
        Row 2 — RpoD reference consensus logo
        Row 3 — HS-3 candidate residues laid out as a per-column track,
                shaded by which panel the residue agrees with (RpoS-like
                green, RpoD-like red, neither grey). Diagnostic columns are
                marked with a red asterisk above the RpoS logo.

    The figure speaks directly to the question the literature considers
    diagnostic: at the columns where the two reference panels disagree, does
    the candidate carry the RpoS or the RpoD consensus residue?
    """
    cols = sorted({c.column_index for c in region.diag_columns})
    if not cols:
        return

    # Extend slightly so the logo has context flanking the diagnostic columns
    pad = 2
    col_min, col_max = max(0, cols[0] - pad), cols[-1] + pad
    window = list(range(col_min, col_max + 1))

    rpos_freq = _freq_matrix(aligned, rpos_ref_ids, window)
    rpod_freq = _freq_matrix(aligned, rpod_ref_ids, window)
    rpos_bits = _bits_matrix(rpos_freq)
    rpod_bits = _bits_matrix(rpod_freq)
    # Use the E. coli RpoS residue number for x-axis labels (anchor of the diag).
    # Build mapping: col → ecoli_rpoS_pos using the RegionDiagnostic columns.
    col_to_pos = {c.column_index: c.ecoli_rpos_pos for c in region.diag_columns}

    # Dedicated marker strip at the top so asterisks never overlap residue
    # letters in the RpoS information-content logo.
    fig = plt.figure(figsize=(max(8, len(window) * 0.32), 6.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[0.30, 3, 3, 1.4], hspace=0.30,
                          left=0.08, right=0.97, top=0.86, bottom=0.18)
    ax_markers = fig.add_subplot(gs[0])
    ax_rpos = fig.add_subplot(gs[1])
    ax_rpod = fig.add_subplot(gs[2])
    ax_hs3 = fig.add_subplot(gs[3])

    # marker strip styling
    ax_markers.set_xlim(-0.5, len(window) - 0.5)
    ax_markers.set_ylim(0, 1)
    ax_markers.set_xticks([])
    ax_markers.set_yticks([])
    for spine in ax_markers.spines.values():
        spine.set_visible(False)

    # --- RpoS logo
    logo_rpos = logomaker.Logo(
        rpos_bits.reset_index(drop=True),
        color_scheme="chemistry", ax=ax_rpos,
        baseline_width=0.3,
    )
    logo_rpos.style_spines(visible=False)
    logo_rpos.style_spines(spines=["left", "bottom"], visible=True)
    ax_rpos.set_ylabel("bits (RpoS refs)", fontsize=10)
    ax_rpos.set_ylim(0, 4.4)
    ax_rpos.set_xticks([])
    ax_rpos.tick_params(axis="y", labelsize=8)
    ax_rpos.set_title("RpoS reference consensus (Group 2 σ⁷⁰)",
                      fontsize=10, loc="left", y=0.94)

    # --- RpoD logo
    logo_rpod = logomaker.Logo(
        rpod_bits.reset_index(drop=True),
        color_scheme="chemistry", ax=ax_rpod,
        baseline_width=0.3,
    )
    logo_rpod.style_spines(visible=False)
    logo_rpod.style_spines(spines=["left", "bottom"], visible=True)
    ax_rpod.set_ylabel("bits (RpoD refs)", fontsize=10)
    ax_rpod.set_ylim(0, 4.4)
    ax_rpod.set_xticks([])
    ax_rpod.tick_params(axis="y", labelsize=8)
    ax_rpod.set_title("RpoD reference consensus (Group 1 σ⁷⁰, housekeeping)",
                      fontsize=10, loc="left", y=0.94)

    # --- HS-3 track
    hs3_seq = aligned[hs3_id]
    verdict_by_col = {c.column_index: c.verdict for c in region.diag_columns}
    palette = {
        "RpoS-like": "#2ca02c",
        "RpoD-like": "#d62728",
        "neither":   "#8c8c8c",
    }
    ax_hs3.set_xlim(-0.5, len(window) - 0.5)
    ax_hs3.set_ylim(0, 1)
    ax_hs3.set_yticks([])
    ax_hs3.tick_params(axis="x", labelsize=8)
    diagnostic_cols_set = {c.column_index for c in region.diag_columns}
    x_tick_positions, x_tick_labels = [], []
    for i, col in enumerate(window):
        is_diag = col in diagnostic_cols_set
        verdict = verdict_by_col.get(col, "neither" if is_diag else None)
        bg = palette.get(verdict, "#f0f0f0")
        edge = "black" if is_diag else "#bbbbbb"
        ax_hs3.add_patch(patches.Rectangle((i - 0.45, 0.05), 0.9, 0.9,
                                            facecolor=bg, edgecolor=edge,
                                            linewidth=0.6))
        res = hs3_seq[col]
        ax_hs3.text(i, 0.50, res if res != GAP else "—",
                    ha="center", va="center",
                    fontsize=10, weight="bold",
                    color="white" if is_diag else "black")
        if col in col_to_pos:
            x_tick_positions.append(i)
            x_tick_labels.append(str(col_to_pos[col]))
    ax_hs3.set_xticks(x_tick_positions)
    ax_hs3.set_xticklabels(x_tick_labels, fontsize=7)
    ax_hs3.set_xlabel("E. coli RpoS residue position",
                     fontsize=9, labelpad=2)
    ax_hs3.set_ylabel(hs3_track_label, fontsize=9, weight="bold")
    for spine in ("top", "right", "left"):
        ax_hs3.spines[spine].set_visible(False)
    ax_hs3.spines["bottom"].set_visible(False)

    # Asterisks go in the marker strip (axes-coords) so they clear the letters.
    for i, col in enumerate(window):
        if col in diagnostic_cols_set:
            ax_markers.plot(i, 0.4, marker="*", markersize=11,
                            color="#d62728", clip_on=False, mec="black",
                            mew=0.5)
    ax_markers.text(-0.5, 0.4, "diagnostic →", ha="right", va="center",
                    fontsize=8, color="#444")

    # Two-row legend strip below the figure: (i) verdict colour key + asterisk,
    # (ii) chemistry colour scheme of the logo letters.
    verdict_handles = [
        Line2D([0], [0], marker="*", color="#d62728", linestyle="None",
               markersize=10, mec="black", mew=0.4,
               label="diagnostic column (RpoS vs RpoD consensus differ)"),
        patches.Patch(facecolor=palette["RpoS-like"], edgecolor="black",
                      label="HS-3 residue matches RpoS consensus"),
        patches.Patch(facecolor=palette["RpoD-like"], edgecolor="black",
                      label="HS-3 residue matches RpoD consensus"),
        patches.Patch(facecolor=palette["neither"], edgecolor="black",
                      label="HS-3 residue matches neither"),
    ]
    leg1 = fig.legend(handles=verdict_handles, loc="lower center",
                      fontsize=8.0, frameon=False, ncol=4,
                      bbox_to_anchor=(0.5, 0.085),
                      columnspacing=1.5, handletextpad=0.4)
    fig.add_artist(leg1)
    chem_handles = [
        patches.Patch(facecolor=col, edgecolor="black", label=lbl)
        for lbl, col in CHEMISTRY_LEGEND
    ]
    fig.legend(handles=chem_handles, loc="lower center",
               fontsize=7.5, frameon=False, ncol=4,
               bbox_to_anchor=(0.5, 0.025), title="Logo-letter colours (chemistry scheme)",
               title_fontsize=7.5,
               columnspacing=1.5, handletextpad=0.4)

    fig.suptitle(region_title, fontsize=12, weight="bold",
                 x=0.08, ha="left", y=0.97)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
