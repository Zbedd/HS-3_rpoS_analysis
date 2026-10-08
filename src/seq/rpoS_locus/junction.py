"""Where pRE118 sits in the HS-3 rpoS- chromosome, to the base.

A single crossover between the chromosome and a plasmid carrying an internal
fragment of the target gene inserts the whole plasmid and duplicates the
fragment, leaving one copy either side of the vector. This module recovers
that structure by exact matching: anchor a unique k-mer in each flank, walk
outward to the first mismatch, and read the boundaries off. The duplicated
segment is the overlap between what the left flank ends with and what the
right flank starts with.

Exact matching rather than alignment because the recombination is exact —
there is no indel to place — and because the cached BLAST+ under
`data/rpos_sequence/_tools/` is a win64 build that does not run here.

The walk is repeated at every anchor length in `config.junction.anchor_k` and
the calls must agree; a repeat or an assembly artefact that swallows one
anchor length shows up as a disagreement instead of a silent wrong answer.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import io
from .paths import Paths


@dataclass(frozen=True)
class Integration:
    """The insertion, in the coordinates of each sequence it touches."""

    # Chromosome (reference genome, 1-based inclusive)
    arm_start: int          # first base of the duplicated homology arm
    arm_end: int            # last base of it
    # pRE118 (1-based inclusive, circular — start > end means it wraps)
    vector_start: int
    vector_end: int
    vector_bp: int
    absent_start: int       # the stretch of pRE118 not in the chromosome:
    absent_end: int         # the site the homology arm was cloned into
    # Mutant contig (1-based), the two chromosome/vector boundaries
    mutant_left: int
    mutant_right: int
    strand: str             # orientation of the locus on the mutant contig

    @property
    def arm_bp(self) -> int:
        return self.arm_end - self.arm_start + 1

    @property
    def absent_bp(self) -> int:
        return self.absent_end - self.absent_start + 1

    def arm_split(self, gene: io.Feature) -> dict:
        """How the arm sits across `gene`.

        `first_cds_base` is where the arm opens inside the CDS, counting the
        A of the start codon as 1; `covers_start_codon` is what decides
        whether the copy downstream of the vector can initiate at all.
        """
        lo = max(self.arm_start, gene.start)
        hi = min(self.arm_end, gene.end)
        coding = max(0, hi - lo + 1)
        return {
            "upstream": max(0, gene.start - self.arm_start),
            "first_cds_base": lo - gene.start + 1,
            "covers_start_codon": self.arm_start <= gene.start,
            "coding": coding,
            "residues": coding // 3,
            "cds_bp": gene.length,
            "cds_residues": gene.length // 3 - 1,   # less the stop codon
        }


def _unique_index(seq: str, kmer: str) -> int:
    """Index of `kmer` in `seq` when it occurs exactly once, else -1."""
    first = seq.find(kmer)
    if first < 0 or seq.find(kmer, first + 1) >= 0:
        return -1
    return first


def _anchor(window: str, target: str, start: int, k: int, step: int = 50) -> tuple:
    """Walk right from `start` until a window k-mer is unique in `target`.

    Returns (window index, target index).
    """
    for i in range(start, len(window) - k, step):
        at = _unique_index(target, window[i:i + k])
        if at >= 0:
            return i, at
    raise ValueError("no unique anchor found")


def _extend_right(a: str, i: int, b: str, j: int) -> int:
    n = 0
    while i + n < len(a) and j + n < len(b) and a[i + n] == b[j + n]:
        n += 1
    return n


def _extend_left(a: str, i: int, b: str, j: int) -> int:
    n = 0
    while i - n - 1 >= 0 and j - n - 1 >= 0 and a[i - n - 1] == b[j - n - 1]:
        n += 1
    return n


def _resolve_at(chromosome: str, window: str, offset: int, vector: str,
                k: int) -> dict:
    """One pass at anchor length `k`. Coordinates are 0-based indices."""
    doubled = vector + vector
    reverse = io.revcomp(vector)
    doubled_rc = reverse + reverse

    # Left flank: anchor in the chromosome, walk to the first mismatch.
    win_i, chrom_i = _anchor(window, chromosome, 200, k)
    shift = chrom_i - win_i
    j = win_i
    while j + 1 < len(window) and window[j] == chromosome[j + shift]:
        j += 1
    left_break = j                       # first window index that is not chromosome

    # Vector: anchor just inside the break, then extend both ways. The vector
    # is circular, so it is searched doubled; a block that crosses the origin
    # is one run in the doubled sequence.
    probe = window[left_break + 5:left_break + 5 + k]
    for strand, target in (("+", doubled), ("-", doubled_rc)):
        at = target.find(probe)
        if at >= 0:
            break
    else:
        raise ValueError(f"no pRE118 match beyond the junction at k={k}")
    back = _extend_left(window, left_break + 5, target, at)
    forward = _extend_right(window, left_break + 5, target, at)
    vector_from = at - back
    vector_span = back + forward
    right_break = left_break + 5 + forward

    # Right flank: anchor beyond the vector, back-extend to the exact edge.
    win_j, chrom_j = _anchor(window, chromosome, right_break + 200, k)
    shift2 = chrom_j - win_j
    back2 = _extend_left(window, win_j, chromosome, win_j + shift2)
    right_start = win_j - back2

    return {
        "arm_start0": right_start + shift2,
        "arm_end0": left_break - 1 + shift,
        "vector_from0": vector_from % len(vector),
        "vector_span": vector_span,
        "vector_strand": strand,
        "left_break_win": left_break,
        "right_start_win": right_start,
        "offset": offset,
    }


def resolve(cfg: dict, paths: Paths, verbose: bool = False) -> Integration:
    """Locate the integration, agreeing across every configured anchor length."""
    chromosome = io.replicon(cfg, paths, "hs3")
    contig = io.mutant_contig(cfg, paths)
    vector, _ = io.vector(cfg, paths)
    lo, hi = cfg["junction"]["window"]
    ks = list(cfg["junction"]["anchor_k"])

    # Orient the window to the reference. The locus reads on the mutant's
    # minus strand, so the window is reverse-complemented and a window index
    # maps back to `hi - index`.
    k0 = ks[0]
    probe = contig[lo + 200:lo + 200 + k0]
    if chromosome.find(probe) >= 0:
        window, strand = contig[lo:hi], "+"
    elif chromosome.find(io.revcomp(probe)) >= 0:
        window, strand = io.revcomp(contig[lo:hi]), "-"
    else:
        raise ValueError("mutant window does not match the reference chromosome")

    calls = [_resolve_at(chromosome, window, lo, vector, k) for k in ks]
    for key in ("arm_start0", "arm_end0", "vector_from0", "vector_span"):
        values = {c[key] for c in calls}
        if len(values) != 1:
            raise ValueError(
                f"anchor lengths {ks} disagree on {key}: {sorted(values)}")
    call = calls[0]

    span, size = call["vector_span"], len(vector)
    vector_start = call["vector_from0"] + 1
    vector_end = (call["vector_from0"] + span - 1) % size + 1
    absent_start = vector_end % size + 1
    absent_end = (vector_start - 2) % size + 1

    def to_contig(win_index: int) -> int:
        return hi - win_index if strand == "-" else lo + win_index + 1

    result = Integration(
        arm_start=call["arm_start0"] + 1,
        arm_end=call["arm_end0"] + 1,
        vector_start=vector_start,
        vector_end=vector_end,
        vector_bp=span,
        absent_start=absent_start,
        absent_end=absent_end,
        mutant_left=to_contig(call["left_break_win"]),
        mutant_right=to_contig(call["right_start_win"]),
        strand=strand,
    )
    if verbose:
        arm = chromosome[result.arm_start - 1:result.arm_end]
        copies = contig.count(arm) + contig.count(io.revcomp(arm))
        print(f"  homology arm  {result.arm_start:,}..{result.arm_end:,} "
              f"({result.arm_bp:,} bp), {copies} copies in the mutant")
        print(f"  pRE118        {result.vector_start:,}..{result.vector_end:,} "
              f"({result.vector_bp:,} of {size:,} bp), cloning site "
              f"{result.absent_start:,}..{result.absent_end:,} absent")
        print(f"  agreed across anchor lengths {ks}")
    return result


def copies_of_arm(cfg: dict, paths: Paths, integration: Integration) -> int:
    """How many times the duplicated arm occurs in the mutant. Must be 2."""
    chromosome = io.replicon(cfg, paths, "hs3")
    contig = io.mutant_contig(cfg, paths)
    arm = chromosome[integration.arm_start - 1:integration.arm_end]
    return contig.count(arm) + contig.count(io.revcomp(arm))


# ── Layout for panel D ────────────────────────────────────────

@dataclass(frozen=True)
class Segment:
    """One drawn block, in bp relative to the WT rpoS start codon."""
    start: int
    end: int
    strand: int
    label: str
    kind: str           # gene | truncated | arm | vector | backbone

    @property
    def length(self) -> int:
        return self.end - self.start


def locus_layout(cfg: dict, paths: Paths, integration: Integration) -> dict:
    """The three tiers panel D draws, all in one coordinate system.

    Origin is the WT rpoS start codon; the mutant tier continues past the
    insertion on the same scale, so the 6 kb the vector adds is the visible
    difference between the two tiers rather than a number in a caption.
    """
    features = io.genome_features(cfg, paths, "hs3")
    rpoS = io.rpoS_cds(cfg, paths)
    origin = rpoS.start
    context = [f for f in features
               if f.gene != "rpoS"
               and f.end >= rpoS.start - 1500 and f.start <= rpoS.end + 1200]

    def rel(value):
        return value - origin

    arm = (rel(integration.arm_start), rel(integration.arm_end) + 1)
    arm_lo, arm_hi = arm
    # What a single crossover inserts: the vector plus one extra copy of the
    # arm. Everything 3' of the arm's start moves by that much.
    shift = arm_hi + integration.vector_bp - arm_lo

    wt = [Segment(rel(f.start), rel(f.end) + 1,
                  1 if f.strand == "+" else -1,
                  io.ortholog_of(cfg, f) or f.gene, "gene")
          for f in context]
    wt.append(Segment(0, rel(rpoS.end) + 1, 1, "rpoS", "gene"))
    wt.sort(key=lambda s: s.start)

    def moved(seg):
        return Segment(seg.start + shift, seg.end + shift, seg.strand,
                       seg.label, seg.kind)

    mutant = []
    for seg in wt:
        if seg.end <= arm_lo:
            mutant.append(seg)
        elif seg.start >= arm_lo:
            mutant.append(moved(seg))
        else:
            # The arm opens inside this gene, so the crossover leaves two
            # partial copies of it: a 5' one ending where the duplication
            # ends, and a 3' one that begins at the arm's own first base —
            # past the start codon, which is why it cannot initiate.
            mutant.append(Segment(seg.start, arm_hi, seg.strand,
                                  seg.label, "truncated"))
            mutant.append(Segment(arm_lo + shift, seg.end + shift, seg.strand,
                                  seg.label, seg.kind))
    # The vector is not a block on this tier: only the genes it carries are
    # drawn, and everything else it contributes is distance.
    markers = (cfg.get("integration_panel") or {}).get("markers") or []
    for label, start, end in vector_blocks(cfg, paths, integration):
        if label in markers:
            mutant.append(Segment(arm_hi + start, arm_hi + end, 1,
                                  label, "marker"))
    mutant.sort(key=lambda s: s.start)

    mutant, gap = _collapse(cfg, paths, integration, mutant,
                            arm_hi + integration.vector_bp)

    # Design tier: the construct linearised at the cloning site, so its copy
    # of the arm sits directly above the chromosome's.
    construct = [
        Segment(arm_lo, arm_hi, 1, "rpoS", "arm"),
        Segment(arm_hi, arm_hi + integration.vector_bp, 1,
                cfg["vector"]["name"], "backbone"),
    ]
    return {
        "break": gap,
        "construct": construct,
        "wt": wt,
        "mutant": mutant,
        "arm": arm,
        "arm_bp": arm_hi - arm_lo,
        "rpoS": rpoS,
        "shift": shift,
        "junction_left": arm_hi,
        "junction_right": arm_hi + integration.vector_bp,
    }


def _collapse(cfg: dict, paths: Paths, integration: Integration,
              segments: list, vector_to: int) -> tuple:
    """Shorten the vector between the last marker drawn and the 3' copy.

    Returns the segments with everything beyond the cut moved left, and the
    (start, end) the break mark occupies — or the segments untouched and None
    when no marker is configured.
    """
    panel = cfg.get("integration_panel") or {}
    ends = [seg.end for seg in segments if seg.kind == "marker"]
    if not ends:
        return segments, None

    cut_from = max(ends) + panel.get("keep_bp", 260)
    break_bp = panel.get("break_bp", 420)
    removed = vector_to - cut_from - break_bp
    if removed <= 0:
        return segments, None

    moved = []
    for seg in segments:
        moved.append(Segment(seg.start - removed, seg.end - removed,
                             seg.strand, seg.label, seg.kind)
                     if seg.start >= vector_to else seg)
    return moved, (cut_from, cut_from + break_bp)


def vector_blocks(cfg: dict, paths: Paths, integration: Integration) -> list:
    """Where each pRE118 feature lands inside the integrated block.

    Offsets are bp from the first base of the vector as it appears in the
    chromosome, which is pRE118 `vector_start` — not pRE118 position 1.
    """
    sequence, features = io.vector(cfg, paths)
    size = len(sequence)
    blocks = []
    for feature in features:
        for lo, hi in feature.parts:
            start = (lo - integration.vector_start) % size
            blocks.append((feature.label, start, start + (hi - lo + 1)))
    return blocks
