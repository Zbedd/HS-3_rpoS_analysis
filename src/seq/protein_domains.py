"""Protein-domain records and database-version checks shared by sequence analyses."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# Pfam / NCBIfam accessions we care about
PF_R11 = "PF03979"    # Sigma70_r1_1 (Group 1 essential, Group 2 absent)
PF_R12 = "PF00140"    # Sigma70_r1_2
PF_R2  = "PF04542"    # Sigma70_r2
PF_R3  = "PF04539"    # Sigma70_r3
PF_R4  = "PF04545"    # Sigma70_r4
PF_R4_2 = "PF08281"   # Sigma70_r4_2

NF_RPOS_PRK = "PRK05657.1"
NF_RPOD_PRK = "PRK05658.1"
TIGR_RPOS = "TIGR02394"
TIGR_RPOD = "TIGR02393"

SIGMA70_PFAMS_INTEREST = {PF_R11, PF_R12, PF_R2, PF_R3, PF_R4, PF_R4_2}
NCBIFAM_INTEREST = {NF_RPOS_PRK, NF_RPOD_PRK, TIGR_RPOS, TIGR_RPOD}


@dataclass
class DomainArchitecture:
    protein_id: str
    length: int
    # One (start, end) span per InterProScan location: dict values are LISTS,
    # never merged to a single min-to-max span, which would hide the internal
    # gap in a two-piece match.
    pfams: dict[str, list[tuple[int, int]]] = field(default_factory=dict)       # σ⁷⁰-family Pfams only (filtered)
    pfams_full: dict[str, list[tuple[int, int]]] = field(default_factory=dict)  # every Pfam-A hit (unfiltered)
    nfams: dict[str, list[tuple[int, int]]] = field(default_factory=dict)
    raw: list[dict] = field(default_factory=list)                               # full hit list, includes non-Pfam

    @property
    def has_r11(self) -> bool:
        return PF_R11 in self.pfams

    @property
    def carries_full_sigma70_core(self) -> bool:
        return (PF_R12 in self.pfams and PF_R2 in self.pfams
                and PF_R3 in self.pfams
                and (PF_R4 in self.pfams or PF_R4_2 in self.pfams))


    def gevin_call(self) -> str:
        if not self.carries_full_sigma70_core:
            return "not σ⁷⁰ family"
        return "Group 1 / RpoD-class" if self.has_r11 else "Group 2 / RpoS-class"


def verify_versions(data: dict, source: Path, cfg: dict) -> None:
    """Hard-fail if a cached InterProScan JSON's software / library versions
    don't match the pins in ``config.yaml → software_versions``.

    Prevents silent drift when EBI ships a newer InterProScan release or when
    the Pfam / NCBIfam member databases are updated server-side. If a mismatch
    is found the pipeline aborts with a message telling the user which pin
    to update (or which cache to regenerate against the pinned toolchain).
    """
    pinned = (cfg or {}).get("software_versions") or {}
    expected_ipr = pinned.get("interproscan")
    expected_pfam = pinned.get("pfam_a")
    expected_ncbifam = pinned.get("ncbifam")
    got_ipr = data.get("interproscan-version")
    if expected_ipr and got_ipr and got_ipr != expected_ipr:
        raise SystemExit(
            f"[interpro] version drift: cache {source.name} was produced with "
            f"InterProScan {got_ipr!r}, but config.yaml pins {expected_ipr!r}. "
            f"Delete the cache and rerun against the pinned release, or update "
            f"config.yaml → software_versions.interproscan."
        )
    for prot in data.get("results", []):
        for m in prot.get("matches", []):
            lib = m.get("signature", {}).get("signatureLibraryRelease", {})
            name = lib.get("library")
            version = lib.get("version")
            if not (name and version):
                continue
            expected = None
            key = None
            if name == "PFAM":
                expected, key = expected_pfam, "pfam_a"
            elif name == "NCBIFAM":
                expected, key = expected_ncbifam, "ncbifam"
            if expected and version != expected:
                raise SystemExit(
                    f"[interpro] library drift in {source.name}: {name} release "
                    f"{version!r} in cache vs pinned {expected!r} "
                    f"(config.yaml → software_versions.{key})."
                )


def architecture_from_iprscan_record(record: dict) -> DomainArchitecture:
    pid = record.get("xref", [{}])[0].get("id", "?")
    seq = record.get("sequence", "")
    arch = DomainArchitecture(protein_id=pid, length=len(seq))
    for m in record.get("matches", []):
        sig = m.get("signature", {})
        acc = (sig.get("accession") or "").strip()
        if not acc:
            continue
        locs = m.get("locations", [])
        if not locs:
            continue
        # Preserve every location as its own span — no widest-span collapse.
        for loc in locs:
            s, e = int(loc.get("start", 0)), int(loc.get("end", 0))
            if acc.startswith("PF"):
                arch.pfams_full.setdefault(acc, []).append((s, e))
            if acc in SIGMA70_PFAMS_INTEREST:
                arch.pfams.setdefault(acc, []).append((s, e))
            if acc in NCBIFAM_INTEREST or acc.startswith("PRK") or acc.startswith("TIGR"):
                arch.nfams.setdefault(acc, []).append((s, e))
            arch.raw.append({"accession": acc, "name": sig.get("name", ""),
                             "description": sig.get("description", ""),
                             "start": s, "end": e})
    return arch
