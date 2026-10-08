"""RpoS/RpoD diagnostics and the four manuscript alignment supplements."""
from __future__ import annotations

from pathlib import Path

from . import align, fetch, figure, interpro, msa_view, report, residues
from .paths import Paths, load_config


def verify_software(cfg: dict) -> None:
    """Hard-fail if any installed tool drifts from ``cfg['software_versions']``.

    Currently pins pyMSAviz (the only pure-Python dependency whose version is
    load-bearing for figure output). InterProScan / Pfam / NCBIfam pins are
    enforced downstream in ``interpro.load_cached_iprscan`` and
    ``interpro.scan_proteins`` against the cached JSON metadata.
    """
    pinned = cfg.get("software_versions") or {}
    expected = pinned.get("pymsaviz")
    if expected:
        import pymsaviz
        got = getattr(pymsaviz, "__version__", None)
        if got and got != expected:
            raise SystemExit(
                f"[pipeline] pyMSAviz version drift: installed {got!r}, "
                f"pinned {expected!r}. Install pinned via "
                f"`pip install pymsaviz=={expected}`, or update "
                f"config.yaml → software_versions.pymsaviz."
            )


def sigma_config(cfg: dict) -> dict:
    """Select the two sigma-factor panels without changing the source config."""
    return {**cfg,
            "references": {p: cfg["references"][p] for p in ("rpoS", "rpoD")},
            "hs3_candidates": {p: cfg["hs3_candidates"][p]
                               for p in ("rpoS", "rpoD")}}


def diagnose(paths: Paths, cfg: dict, hs3_seqs: dict, refs: dict) -> dict:
    """The existing residue diagnostic over all four configured sigma regions."""
    sigma_aln = align.parse_aligned_fasta(
        align.aligned_fasta(paths.alignments_root, "sigma70")
    )
    ecoli_rpos_id = cfg["diagnostic_regions"]["ecoli_rpoS_uniprot"] + "_RpoS"
    hs3_rpos_aln_id = hs3_seqs["rpoS"][0] + "_HS3rpoS"
    rpos_ref_ids = [rid + "_RpoS" for rid, _, _ in refs["rpoS"]]
    rpod_ref_ids = [rid + "_RpoD" for rid, _, _ in refs["rpoD"]]
    region_keys = [k for k in cfg["diagnostic_regions"] if k.startswith("region_")]
    region_diagnostics = {
        key: residues.build_diagnostic(
            sigma_aln, ecoli_rpos_id, hs3_rpos_aln_id,
            rpos_ref_ids, rpod_ref_ids,
            cfg["diagnostic_regions"][key]["name"],
            tuple(cfg["diagnostic_regions"][key]["ecoli_rpoS_range"]),
            cfg["diagnostic"]["consensus_min"],
        )
        for key in region_keys
    }
    return region_diagnostics


def render(paths: Paths, cfg: dict, hs3_seqs: dict, refs: dict,
           arches: dict, region_diagnostics: dict, *, exploratory: bool = False) -> None:
    """Write the four published MSA panels; additional views are opt-in."""
    cfg = sigma_config(cfg)
    ref_labels, hs3_labels, _ = msa_view.build_labels(paths, cfg)
    sig_pres, _ = msa_view.build_sigma70_presentation(paths, cfg, ref_labels, hs3_labels)
    region_panel_titles = {
        "region_1_2": "Sigma70_r1.2 (PF00140)",
        "region_2":   "Sigma70_r2 (PF04542)",
        "region_3":   "Sigma70_r3 (PF04539)",
        "region_4":   "Sigma70_r4 (PF04545)",
    }
    for region_key, diag in region_diagnostics.items():
        if not diag.diag_columns:
            continue
        cols = sorted({c.column_index for c in diag.diag_columns})
        figure.render_msa_panel(
            sig_pres, [c.column_index for c in diag.diag_columns],
            paths.output_rpoS / f"panelB_{region_key}_msa.png",
            title=region_panel_titles.get(region_key, diag.region_name),
            start_col=max(0, cols[0] - 5), end_col=cols[-1] + 5,
        )
    if not exploratory:
        return
    from . import logos
    sigma_rows = interpro.architecture_rows(paths, cfg, ("rpoS", "rpoD"),
                                            arches=arches)
    figure.render_domain_cartoon(
        sigma_rows,
        paths.output_rpoS / "panelA_domain_architecture.png",
    )
    figure.render_domain_cartoon(
        sigma_rows,
        paths.output_rpoS / "panelA_domain_architecture_unfiltered.png",
        use_full_pfams=True,
    )

    sigma_aln = align.parse_aligned_fasta(
        align.aligned_fasta(paths.alignments_root, "sigma70"))
    rpos_ref_ids = [rid + "_RpoS" for rid, _, _ in refs["rpoS"]]
    rpod_ref_ids = [rid + "_RpoD" for rid, _, _ in refs["rpoD"]]
    hs3_rpos_aln_id = hs3_seqs["rpoS"][0] + "_HS3rpoS"
    region_logo_titles = {
        "region_1_2": "Panel C (i) · Information-content logos at σ region 1.2 — "
                      "RpoS vs RpoD consensus with HS-3 rpoS candidate residues",
        "region_2":   "Panel C (ii) · Information-content logos at σ region 2 "
                      "(incl. subregion 2.4, −10 recognition) — RpoS vs RpoD "
                      "consensus with HS-3 rpoS candidate residues",
        "region_3":   "Panel C (iii) · Information-content logos at σ region 3 — "
                      "RpoS vs RpoD consensus with HS-3 rpoS candidate residues",
        "region_4":   "Panel C (iv) · Information-content logos at σ region 4 "
                      "(incl. subregion 4.2, −35 recognition) — RpoS vs RpoD "
                      "consensus with HS-3 rpoS candidate residues",
    }
    for region_key, diag in region_diagnostics.items():
        if not diag.diag_columns:
            continue
        logos.render_sigma_diagnostic_logos(
            sigma_aln, diag, rpos_ref_ids, rpod_ref_ids, hs3_rpos_aln_id,
            paths.output_rpoS / f"panelC_{region_key}_sequence_logos.png",
            region_logo_titles.get(region_key, diag.region_name),
            hs3_track_label="HS-3 rpoS",
        )


def run(workdir: Path, *, exploratory: bool = False) -> None:
    """Run the manuscript RpoS alignment analysis."""
    paths = Paths(workdir=workdir)
    cfg = sigma_config(load_config(workdir))
    verify_software(cfg)
    for directory in (paths.output_rpoS, paths.queries_root,
                      paths.shared_queries_root, paths.alignments_root,
                      paths.raw_hmmer):
        directory.mkdir(parents=True, exist_ok=True)

    hs3_seqs = {}
    for role, acc in cfg["hs3_candidates"].items():
        path = fetch.extract_hs3_protein(paths, acc)
        _, sequence = fetch.read_fasta_record(path)
        hs3_seqs[role] = (acc, sequence)
    refs = {}
    for panel, entries in cfg["references"].items():
        refs[panel] = []
        for entry in entries:
            _, sequence = fetch.read_fasta_record(fetch.get_reference(paths, entry))
            refs[panel].append((entry["id"], entry["name"], sequence))

    sigma_inputs = [(hs3_seqs[role][0] + "_HS3" + role, hs3_seqs[role][1])
                    for role in ("rpoS", "rpoD")]
    for panel, suffix in (("rpoS", "_RpoS"), ("rpoD", "_RpoD")):
        sigma_inputs.extend((rid + suffix, sequence) for rid, _, sequence in refs[panel])
    align.align(align.build_fasta_input(sigma_inputs), paths.alignments_root, "sigma70")

    records = interpro.load_cached_iprscan(paths, cfg)
    records.update(interpro.load_rq2_iprscan(paths, cfg))
    proteins = {acc: sequence for acc, sequence in hs3_seqs.values()}
    proteins.update({rid: sequence for entries in refs.values()
                     for rid, _, sequence in entries})
    if set(proteins) - records.keys():
        # Keep independent submissions separate from the combined workflow's cache.
        records.update(interpro.scan_proteins(paths, proteins, "rpos_refs", cfg=cfg))
    arches = {pid: interpro.architecture_from_iprscan_record(record)
              for pid, record in records.items()}
    diagnostics = diagnose(paths, cfg, hs3_seqs, refs)
    render(paths, cfg, hs3_seqs, refs, arches, diagnostics, exploratory=exploratory)
    report.emit_rpoS(
        paths=paths, cfg=cfg,
        hs3_rpos_acc=hs3_seqs["rpoS"][0], hs3_rpod_acc=hs3_seqs["rpoD"][0],
        hs3_rpos_arch=arches[hs3_seqs["rpoS"][0]],
        hs3_rpod_arch=arches[hs3_seqs["rpoD"][0]],
        ref_arches={p: [arches[rid] for rid, _, _ in refs[p]]
                    for p in ("rpoS", "rpoD")},
        region_diagnostics=diagnostics,
    )
    print(f"RpoS alignments and diagnostics written to {paths.output_rpoS}")
