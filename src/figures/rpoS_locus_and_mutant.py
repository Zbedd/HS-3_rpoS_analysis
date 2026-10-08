"""The manuscript's four-panel HS-3 RpoS figure.

    A  chromosome and sigma-factor loci
    B  rpoS neighbourhood in HS-3 and three reference organisms
    C  RpoS domain architecture
    D  pRE118 integration
"""
from __future__ import annotations

from pathlib import Path

import viz as fs  # selects the Agg backend on import

from seq import domain_panels
from seq.rpos_alignments import interpro
from seq.rpos_alignments.paths import Paths as Rq2Paths, load_config as rq2_config
from seq.rpoS_locus import io, junction, panels

ROWS = [{"locus": 2.25}, {"domains": 2.00}, {"construct": 2.30}]
RING_SPLIT = [0.30, 0.70]


def plot(out_path: Path, workdir: Path = None) -> list:
    """Draw panels A–D. Returns the paths written."""
    cfg, paths = io.load(workdir)

    fig, page = fs.stack(fs.COL_DOUBLE, ROWS, pad=0.20)
    axes = []

    left, right = fs.split(page["locus"], RING_SPLIT)
    ax_chrom = fig.add_subplot(left)
    panels.draw_chromosome(
        ax_chrom, io.replicon_length(cfg, paths, "hs3"),
        io.sigma_factors(cfg, paths, "hs3"),
        focus=cfg["sigma_factors"]["focus"],
        housekeeping=cfg["sigma_factors"]["housekeeping"])

    ax_syn = fig.add_subplot(right)
    panels.draw_synteny(ax_syn, io.synteny_tracks(cfg, paths),
                        flanks=(cfg["neighbourhood"]["flank_upstream_bp"],
                                cfg["neighbourhood"]["flank_downstream_bp"]))
    axes += [ax_chrom, ax_syn]

    rq2_paths = Rq2Paths(paths.workdir)
    rq2_cfg = rq2_config(paths.workdir)
    arches = interpro.architecture_rows(rq2_paths, rq2_cfg, ("rpoS",))

    ax_dom = fig.add_subplot(page["domains"])
    domain_panels.draw_domain_cartoon(ax_dom, arches,
                                     type_scale=domain_panels.HOUSE_TYPE,
                                     emphasis="label", label_mode="legend",
                                     align_on="PF00140", shape="arrow",
                                     show_axis=False, scale_bar=(100, "100 aa"))
    axes.append(ax_dom)

    integration = junction.resolve(cfg, paths)
    _, vector_features = io.vector(cfg, paths)
    palette = panels.vector_palette(f.label for f in vector_features)

    ax_int = fig.add_subplot(page["construct"])
    panels.draw_integration(
        ax_int, junction.locus_layout(cfg, paths, integration),
        palette=palette, vector_name=cfg["vector"]["name"])
    axes.append(ax_int)

    fs.tag_panels(axes, labels="ABCD", align=True)
    return fs.save(fig, out_path)
