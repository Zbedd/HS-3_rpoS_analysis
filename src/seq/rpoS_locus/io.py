"""Genome, annotation and vector readers. The only place a coordinate,
a strand or an ortholog assignment is defined.

Coordinates are 1-based inclusive throughout, matching GFF3 and GenBank, so
a number printed on a panel is the number in the annotation file. Sequences
are returned as plain upper-case strings on the plus strand of the replicon.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass, replace
from pathlib import Path
from urllib.parse import unquote

from .paths import Paths, load_config

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


# ── Sequence ──────────────────────────────────────────────────

def read_fasta(path: Path) -> dict:
    """Records keyed on the first token of the header."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records, name, chunks = {}, None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return {k: v.upper() for k, v in records.items()}


def replicon(cfg: dict, paths: Paths, tag: str) -> str:
    """The chromosome sequence for one genome tag."""
    name = cfg["assemblies"][tag]["replicon"]
    records = read_fasta(paths.genome_fna(tag))
    if name not in records:
        raise KeyError(f"{name} not in {paths.genome_fna(tag)}; "
                       f"have {sorted(records)}")
    return records[name]


def mutant_contig(cfg: dict, paths: Paths) -> str:
    records = read_fasta(paths.data(cfg["mutant"]["contigs"]))
    return records[cfg["mutant"]["contig"]]


# ── Annotation ────────────────────────────────────────────────

@dataclass(frozen=True)
class Feature:
    seqid: str
    start: int          # 1-based inclusive
    end: int
    strand: str
    gene: str = ""
    product: str = ""
    locus_tag: str = ""
    protein_id: str = ""
    pseudo: bool = False    # RefSeq calls the reading frame broken; there is
                            # no protein, so nothing can be aligned to it

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2

    @property
    def name(self) -> str:
        return self.gene or self.locus_tag or self.protein_id


def _attributes(field: str) -> dict:
    out = {}
    for item in field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = unquote(value)
    return out


def read_gff(path: Path, seqid: str = None, kinds=("CDS",)) -> list:
    """Features of the given kinds, in coordinate order."""
    features = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 9 or row[2] not in kinds:
                continue
            if seqid is not None and row[0] != seqid:
                continue
            attrs = _attributes(row[8])
            features.append(Feature(
                seqid=row[0], start=int(row[3]), end=int(row[4]), strand=row[6],
                gene=attrs.get("gene", ""), product=attrs.get("product", ""),
                locus_tag=attrs.get("locus_tag", ""),
                protein_id=attrs.get("protein_id", ""),
                pseudo=attrs.get("pseudo", "") == "true",
            ))
    features.sort(key=lambda f: f.start)
    return features


def genome_features(cfg: dict, paths: Paths, tag: str) -> list:
    return read_gff(paths.genome_gff(tag), cfg["assemblies"][tag]["replicon"])


def replicon_length(cfg: dict, paths: Paths, tag: str) -> int:
    """Length from the GFF ##sequence-region line for the named replicon."""
    name = cfg["assemblies"][tag]["replicon"]
    with open(paths.genome_gff(tag), encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4 and parts[1] == name:
                    return int(parts[3])
            elif not line.startswith("#"):
                break
    raise ValueError(f"no ##sequence-region for {name}")


# ── Panel A: the σ factor complement ──────────────────────────

def sigma_factors(cfg: dict, paths: Paths, tag: str = "hs3") -> list:
    """(display name, Feature) for each σ factor, in config order.

    Raises if a configured factor is missing or ambiguous — a σ factor that
    silently vanishes from the panel would read as a genome that lacks it.
    """
    patterns = cfg["sigma_factors"]["patterns"]
    features = genome_features(cfg, paths, tag)
    found = []
    for name, product in patterns.items():
        hits = [f for f in features if f.product == product]
        if len(hits) != 1:
            raise ValueError(
                f"{tag}: expected exactly one {name} ({product!r}), got {len(hits)}")
        found.append((name, hits[0]))
    return found


def rpoS_cds(cfg: dict, paths: Paths) -> Feature:
    """The rpoS CDS from its initiation codon, not from RefSeq's.

    RefSeq starts the CDS at an internal ATG; `config.rpoS_start` carries the
    real one and the reason. Anything that needs the reading frame — the
    integration layout, where the arm falls inside the gene — resolves here,
    so the correction is made once.
    """
    feature = next(f for f in genome_features(cfg, paths, "hs3")
                   if f.gene == "rpoS")
    start = cfg.get("rpoS_start")
    if start is None or start == feature.start:
        return feature
    if start >= feature.start or (feature.start - start) % 3:
        raise ValueError(
            f"rpoS_start {start:,} must be upstream of the annotated "
            f"{feature.start:,} and in frame with it")
    return replace(feature, start=start)


# ── Panel B: a σ factor neighbourhood ─────────────────────────

def locus(cfg: dict, key: str = None) -> dict:
    """The spec for one σ factor locus.

    `None` is rpoS's own, which the manuscript figure draws; a key names one
    of the four the supplement draws. A spec carries the anchor, the flanks
    and the ortholog groups, so every locus reaches the same readers below.
    """
    if key is None or key == cfg["neighbourhood"]["anchor"]:
        return cfg["neighbourhood"]
    return cfg["alt_factors"]["loci"][key]


def _group_key(groups: list, feature: Feature) -> str:
    for group in groups:
        if feature.gene and feature.gene in group["symbols"]:
            return group["key"]
        if any(p.lower() in feature.product.lower() for p in group["products"]):
            return group["key"]
    return ""


def neighbourhood(cfg: dict, paths: Paths, tag: str, key: str = None) -> tuple:
    """(anchor, features) — every CDS within the configured flank of the anchor.

    Features carry an `ortholog` attribute via `ortholog_of`. The anchor is
    matched on its ortholog group rather than on a gene symbol: only K-12
    names σ54 `rpoN`, and the other three annotate it by product alone.
    """
    spec = locus(cfg, key)
    anchor_key = spec["anchor"]
    groups = spec["ortholog_groups"]
    # Pulled symmetrically at the wider of the two reaches; `synteny_tracks`
    # trims to the real flanks once the window has been oriented, since which
    # side is upstream depends on the anchor's strand.
    flank = max(spec["flank_upstream_bp"], spec["flank_downstream_bp"])
    features = genome_features(cfg, paths, tag)
    anchors = [f for f in features if _group_key(groups, f) == anchor_key]
    if len(anchors) != 1:
        raise ValueError(f"{tag}: expected one {anchor_key}, got {len(anchors)}")
    anchor = anchors[0]
    window = [f for f in features
              if f.end >= anchor.start - flank and f.start <= anchor.end + flank]
    return anchor, window


def ortholog_of(cfg: dict, feature: Feature, key: str = None) -> str:
    """Ortholog group key for a feature, or "" if it is not in the core set."""
    return _group_key(locus(cfg, key)["ortholog_groups"], feature)


def locus_protein(cfg: dict, paths: Paths, tag: str, key: str) -> tuple:
    """(Feature, sequence) for one genome's copy of the σ factor itself."""
    anchor, _ = neighbourhood(cfg, paths, tag, key)
    faa = read_fasta(paths.assembly_dir(tag) / "protein.faa")
    if anchor.protein_id not in faa:
        raise KeyError(f"{tag}: {anchor.protein_id} not in protein.faa")
    return anchor, faa[anchor.protein_id]


# ── Vector sequence and annotation ────────────────────────────

@dataclass(frozen=True)
class VectorFeature:
    label: str
    parts: tuple        # ((start, end), ...) 1-based; more than one when the
                        # feature wraps the origin, as R6K γ ori does
    strand: int

    @property
    def start(self) -> int:
        return self.parts[0][0]

    @property
    def end(self) -> int:
        return self.parts[-1][1]

    @property
    def length(self) -> int:
        return sum(hi - lo + 1 for lo, hi in self.parts)


def vector(cfg: dict, paths: Paths) -> tuple:
    """(sequence, features) for the sequenced pRE118.

    Blocks are the labels named in the config, in config order; the GenBank
    also carries fragment annotations, which are not drawn.
    """
    from Bio import SeqIO

    record = SeqIO.read(paths.data(cfg["vector"]["genbank"]), "genbank")
    wanted = cfg["vector"]["features"]
    by_label = {}
    for feature in record.features:
        for label in feature.qualifiers.get("label", []):
            if label in wanted and label not in by_label:
                by_label[label] = VectorFeature(
                    label=label,
                    parts=tuple((int(part.start) + 1, int(part.end))
                                for part in feature.location.parts),
                    strand=int(feature.location.strand or 1),
                )
    missing = [l for l in wanted if l not in by_label]
    if missing:
        raise ValueError(f"{cfg['vector']['genbank']}: no feature labelled {missing}")
    return str(record.seq).upper(), [by_label[l] for l in wanted]


def construct(cfg: dict, paths: Paths, integration) -> tuple:
    """(size, features) for pRE118 carrying the rpoS arm.

    The construct was never sequenced, but every coordinate in it is known:
    the arm replaced the stretch of pRE118 the integration shows to be absent,
    so the map is the sequenced vector with that stretch swapped out. Drawing
    the empty vector instead would show the arm as the 35 bp site it displaced
    rather than at its own 497 bp.
    """
    sequence, features = vector(cfg, paths)
    size = len(sequence)
    cut_lo, cut_hi = integration.absent_start, integration.absent_end
    arm_bp = integration.arm_bp
    grew = arm_bp - (cut_hi - cut_lo + 1)

    def moved(position):
        return position + grew if position >= cut_hi else position

    placed = [VectorFeature(label=f.label,
                            parts=tuple((moved(lo), moved(hi))
                                        for lo, hi in f.parts),
                            strand=f.strand)
              for f in features]
    placed.append(VectorFeature(label=ARM_LABEL,
                                parts=((cut_lo, cut_lo + arm_bp - 1),),
                                strand=1))
    return size + grew, placed


# The arm's name in both the integration and standalone construct maps.
ARM_LABEL = "rpoS'"


# ── Convenience ───────────────────────────────────────────────

def load(workdir: Path = None) -> tuple:
    """(config, Paths) — the pair every entry point starts from."""
    paths = Paths(workdir) if workdir else Paths()
    return load_config(), paths


@dataclass(frozen=True)
class Gene:
    """One gene arrow in a synteny track, positioned relative to the anchor."""
    start: int          # bp from the anchor start codon, anchor reading left to right
    end: int
    strand: int         # +1 / -1 after the window has been oriented
    ortholog: str       # core-set key, or "" for context
    label: str
    local: str = ""     # what this genome calls a core gene, where that
                        # differs from the group key — algU for PA14's rpoE
    pseudo: bool = False

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2


@dataclass(frozen=True)
class Track:
    tag: str
    label: str
    genes: tuple


def synteny_tracks(cfg: dict, paths: Paths, key: str = None) -> list:
    """One track per comparator, every window oriented so the anchor reads
    rightward.

    Some genomes carry the anchor on the minus strand — rpoS in K-12 and
    V. cholerae, rpoD in HS-3 and PA14 — so those windows are flipped;
    without that the conserved block would read backwards against the others
    and the panel would show a difference that is not there.
    """
    spec = locus(cfg, key)
    groups = spec["ortholog_groups"]
    aliases = spec.get("aliases") or {}
    up, down = spec["flank_upstream_bp"], spec["flank_downstream_bp"]
    tracks = []
    for tag in cfg["comparators"]:
        anchor, window = neighbourhood(cfg, paths, tag, key)
        flip = anchor.strand == "-"
        origin = anchor.end if flip else anchor.start
        genes = []
        for feature in window:
            start, end = feature.start - origin, feature.end - origin
            strand = 1 if feature.strand == "+" else -1
            if flip:
                start, end, strand = -end, -start, -strand
            group = _group_key(groups, feature)
            genes.append(Gene(start, end, strand, group,
                              group or feature.gene or "",
                              aliases.get(tag, {}).get(group, feature.gene)
                              if group else "",
                              feature.pseudo))
        genes = [g for g in genes
                 if g.end > -up and g.start < anchor.length + down]
        tracks.append(Track(tag=tag, label=cfg["assemblies"][tag]["short"],
                            genes=tuple(sorted(genes, key=lambda g: g.start))))
    return tracks
