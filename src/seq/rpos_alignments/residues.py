"""Column-by-column RpoS-vs-RpoD residue diagnostic.

Approach (column-based; no pre-published residue list needed):
    1. Anchor on E. coli RpoS (UniProt P13445, 330 aa) — its sequence
       coordinates define the σ region 2 and region 4 ranges (config).
    2. Walk through every alignment column. Skip columns where the
       E. coli RpoS row is a gap. For each remaining column:
         a. Read the residue in each RpoS reference (≥ N) and each RpoD
            reference (≥ N).
         b. If RpoS references reach majority on residue R_S AND RpoD
            references reach majority on residue R_D AND R_S ≠ R_D,
            the column is "diagnostic" for the RpoS-vs-RpoD distinction.
         c. Check what residue the HS-3 candidate carries: R_S, R_D, or
            something else.
    3. Tabulate counts within region 2 (subregion 2.4 sits at the C-end)
       and region 4 (subregion 4.2 sits at the C-end).
"""
from __future__ import annotations

from dataclasses import dataclass

# Position-to-column map sentinel: -1 means "this E. coli RpoS residue
# falls in a gap column" (cannot happen at our anchor since E. coli RpoS
# is itself the reference, but defensive).
GAP = "-"


@dataclass
class DiagColumn:
    column_index: int           # 0-based alignment column
    ecoli_rpos_pos: int         # 1-based residue in E. coli RpoS
    rpos_consensus: str         # majority residue across RpoS panel
    rpod_consensus: str         # majority residue across RpoD panel
    hs3_residue: str            # residue in the HS-3 candidate
    hs3_residue_pos: int        # 1-based residue position in HS-3 candidate
    verdict: str                # "RpoS-like", "RpoD-like", "neither"


@dataclass
class RegionDiagnostic:
    region_name: str
    ecoli_rpos_range: tuple[int, int]
    diag_columns: list[DiagColumn]

    @property
    def n_total(self) -> int:
        return len(self.diag_columns)

    @property
    def n_rpos_like(self) -> int:
        return sum(1 for c in self.diag_columns if c.verdict == "RpoS-like")

    @property
    def n_rpod_like(self) -> int:
        return sum(1 for c in self.diag_columns if c.verdict == "RpoD-like")

    @property
    def n_neither(self) -> int:
        return sum(1 for c in self.diag_columns if c.verdict == "neither")


def _majority(residues: list[str], min_count: int) -> str | None:
    counts: dict[str, int] = {}
    for r in residues:
        if r == GAP:
            continue
        counts[r] = counts.get(r, 0) + 1
    if not counts:
        return None
    res, n = max(counts.items(), key=lambda kv: kv[1])
    return res if n >= min_count else None


def _ungapped_position(aligned_seq: str, column_index: int) -> int:
    """Return the 1-based residue position in the ungapped sequence
    corresponding to ``column_index``. Returns 0 if the column is a gap."""
    if aligned_seq[column_index] == GAP:
        return 0
    return column_index + 1 - aligned_seq[:column_index + 1].count(GAP)


def build_diagnostic(
    aligned: dict[str, str],
    ecoli_rpos_id: str,
    hs3_id: str,
    rpos_ids: list[str],
    rpod_ids: list[str],
    region_name: str,
    ecoli_rpos_range: tuple[int, int],
    consensus_min: int,
) -> RegionDiagnostic:
    """Walk the alignment within the E. coli RpoS region of interest and
    extract diagnostic columns."""
    ecoli_seq = aligned[ecoli_rpos_id]
    hs3_seq = aligned[hs3_id]
    n_cols = len(ecoli_seq)

    diag: list[DiagColumn] = []
    ecoli_resnum = 0
    for col in range(n_cols):
        if ecoli_seq[col] != GAP:
            ecoli_resnum += 1
        if ecoli_seq[col] == GAP:
            continue
        if ecoli_resnum < ecoli_rpos_range[0] or ecoli_resnum > ecoli_rpos_range[1]:
            continue
        rpos_res = [aligned[i][col] for i in rpos_ids]
        rpod_res = [aligned[i][col] for i in rpod_ids]
        rpos_c = _majority(rpos_res, consensus_min)
        rpod_c = _majority(rpod_res, consensus_min)
        if rpos_c is None or rpod_c is None or rpos_c == rpod_c:
            continue
        hs3_r = hs3_seq[col]
        if hs3_r == GAP:
            continue
        if hs3_r == rpos_c:
            verdict = "RpoS-like"
        elif hs3_r == rpod_c:
            verdict = "RpoD-like"
        else:
            verdict = "neither"
        diag.append(DiagColumn(
            column_index=col, ecoli_rpos_pos=ecoli_resnum,
            rpos_consensus=rpos_c, rpod_consensus=rpod_c,
            hs3_residue=hs3_r,
            hs3_residue_pos=_ungapped_position(hs3_seq, col),
            verdict=verdict,
        ))
    return RegionDiagnostic(
        region_name=region_name,
        ecoli_rpos_range=ecoli_rpos_range,
        diag_columns=diag,
    )
