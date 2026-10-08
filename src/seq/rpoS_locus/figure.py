"""Standalone figures, one panel each.

These exist to check a panel in isolation at the size it will be printed.
The manuscript figure composes the same panel functions — see
`figures.rpoS_locus_and_mutant`.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

import viz as fs

from .. import domain_panels
from . import domains, io, junction, panels
from .paths import Paths


def chromosome(cfg: dict, paths: Paths, out: Path) -> list:
    fig, ax = plt.subplots(figsize=(fs.COL_SINGLE, 3.0))
    panels.draw_chromosome(
        ax, io.replicon_length(cfg, paths, "hs3"),
        io.sigma_factors(cfg, paths, "hs3"),
        focus=cfg["sigma_factors"]["focus"],
        housekeeping=cfg["sigma_factors"]["housekeeping"],
        note="plus pJHS3, 2,078 bp")
    return fs.save(fig, out)


def synteny(cfg: dict, paths: Paths, out: Path) -> list:
    fig, ax = plt.subplots(figsize=(fs.COL_DOUBLE, 2.3))
    panels.draw_synteny(ax, io.synteny_tracks(cfg, paths),
                        flanks=(cfg['neighbourhood']['flank_upstream_bp'],
                                cfg['neighbourhood']['flank_downstream_bp']))
    return fs.save(fig, out)


def integration(cfg: dict, paths: Paths, out: Path) -> list:
    integ = junction.resolve(cfg, paths)
    _, features = io.vector(cfg, paths)
    fig, ax = plt.subplots(figsize=(fs.COL_DOUBLE, 2.6))
    panels.draw_integration(
        ax, junction.locus_layout(cfg, paths, integ),
        palette=panels.vector_palette(f.label for f in features),
        vector_name=cfg["vector"]["name"])
    return fs.save(fig, out)


def plasmid(cfg: dict, paths: Paths, out: Path) -> list:
    integ = junction.resolve(cfg, paths)
    _, vector_features = io.vector(cfg, paths)
    size, features = io.construct(cfg, paths, integ)
    fig, ax = plt.subplots(figsize=(fs.COL_SINGLE, 3.0))
    panels.draw_plasmid(
        ax, size, features,
        name=panels.construct_name(cfg["vector"]["name"]),
        palette=panels.construct_palette(vector_features))
    return fs.save(fig, out)


# ── The supplement: the four other σ factors ──────────────────
# The two functions below draw into an axes rather than making a figure, so
# the standalone panels here and the composed page in
# `figures.alt_factors_locus` read the same config and cannot disagree about
# a flank, a palette or which domain a panel is aligned on.


def draw_alt_synteny(ax, cfg: dict, paths: Paths, key: str) -> None:
    """One locus of the supplement: its neighbourhood across the four genomes."""
    spec = io.locus(cfg, key)
    order = [group["key"] for group in spec["ortholog_groups"]]
    panels.draw_synteny(
        ax, io.synteny_tracks(cfg, paths, key),
        flanks=(spec["flank_upstream_bp"], spec["flank_downstream_bp"]),
        anchor=spec["anchor"],
        fills=panels.locus_palette(order, spec["anchor"]),
        name_locals=True)


def draw_alt_domains(ax, cfg: dict, paths: Paths, key: str,
                     arches: dict) -> None:
    """The same locus as a protein: the σ factor's domains in all four.

    Every Pfam-A signature is drawn, not the σ70 regions alone — RpoN carries
    none of those, and the non-essential region that marks RpoD as a Group 1
    factor is not one of them either.
    """
    domain_panels.draw_domain_cartoon(
        ax, domains.architecture_rows(cfg, paths, key, arches),
        type_scale=domain_panels.HOUSE_TYPE, emphasis="label",
        label_mode="legend", align_on=cfg["alt_factors"]["align_on"][key],
        shape="arrow", show_axis=False, scale_bar=(100, "100 aa"),
        use_full_pfams=True, legend_ncol=3)


def alt_synteny(cfg: dict, paths: Paths, key: str, out: Path) -> list:
    fig, ax = plt.subplots(figsize=(fs.COL_ONEHALF, 2.2))
    draw_alt_synteny(ax, cfg, paths, key)
    return fs.save(fig, out)


def alt_domains(cfg: dict, paths: Paths, key: str, arches: dict,
                out: Path) -> list:
    fig, ax = plt.subplots(figsize=(fs.COL_SINGLE, 2.2))
    draw_alt_domains(ax, cfg, paths, key, arches)
    return fs.save(fig, out)


def render_alt(workdir: Path = None, arches: dict = None) -> None:
    """Each supplement panel on its own, at the size it will be printed."""
    cfg, paths = io.load(workdir)
    paths.alt_output_root.mkdir(parents=True, exist_ok=True)
    if arches is None:
        arches = domains.ensure_scan(cfg, paths)
    for key in cfg["alt_factors"]["order"]:
        for name, render in (("synteny", lambda o: alt_synteny(cfg, paths, key, o)),
                             ("domains", lambda o: alt_domains(cfg, paths, key,
                                                               arches, o))):
            written = render(paths.alt_output_root / f"{key}_{name}")
            print("wrote " + ", ".join(str(p) for p in written))


def render_all(workdir: Path = None) -> None:
    cfg, paths = io.load(workdir)
    paths.ensure()
    out = paths.output_root
    for name, render in (("chromosome", chromosome), ("synteny", synteny),
                         ("integration", integration), ("plasmid", plasmid)):
        written = render(cfg, paths, out / name)
        print("wrote " + ", ".join(str(p) for p in written))
