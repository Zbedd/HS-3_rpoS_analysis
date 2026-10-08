"""Single source of truth for protein-sequence retrieval.

One fetcher, one cache, one validation path. Every figure in every pipeline
resolves its sequences through :func:`get_protein`, so a given accession maps
to exactly one file on disk and exactly one sequence in memory.

Design rules
------------
1. **Cache-first.** ``get_protein`` returns the cached record if present and
   valid; otherwise it downloads, validates, then caches. Callers never do
   their own existence checks.
2. **Validate on every path.** The same checks run against freshly downloaded
   bytes *and* against bytes loaded from cache, so a file that was poisoned by
   an earlier buggy run cannot survive by sitting in the cache.
3. **Never cache a bad payload.** Validation happens before the write. A
   deleted UniProt entry returns HTTP 200 with an empty body; that used to be
   written to disk as a 0-byte file and silently folded into the BLAST input
   as an empty-sequence record. It now raises.
4. **Keep the description.** The FASTA description carries the protein name
   and organism. It is retained on the record so callers can cross-check an
   accession against the label a config claims for it — see
   :func:`check_label` and ``python -m seq.seqfetch``.
"""
from __future__ import annotations

import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "ProteinRecord",
    "LabelCheck",
    "get_protein",
    "check_label",
    "infer_source",
    "DEFAULT_CACHE",
]

USER_AGENT = "HS3-seqfetch/1.0 (research use)"

# The one canonical cache. Both RQ1 and RQ2 resolve here, so an accession used
# by both pipelines is fetched once and read from the same bytes by both.
DEFAULT_CACHE = Path("data/rpos_sequence/queries")

# RefSeq/GenBank protein accession prefixes -> NCBI efetch. Everything else is
# assumed to be a UniProtKB accession.
_NCBI_PREFIXES = ("WP_", "NP_", "XP_", "YP_", "AP_", "ZP_", "EP_")


@dataclass(frozen=True)
class ProteinRecord:
    """One protein sequence plus the provenance needed to audit it."""

    accession: str      # accession as requested by the caller
    source: str         # "uniprot" | "ncbi"
    header: str         # full FASTA header line, '>' stripped
    description: str    # header with the leading identifier token removed
    sequence: str       # ungapped one-letter sequence
    path: Path          # cache file this record came from

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class LabelCheck:
    """Result of cross-checking a config's human label against the record."""

    accession: str
    declared: str        # the name the config claims
    actual: str          # the description the database returned
    organism_ok: bool | None   # None when the label carries no organism
    gene_ok: bool | None       # None when no gene/protein token was derivable
    verdict: str         # "ok" | "suspect" | "unverifiable"
    reason: str


# --------------------------------------------------------------- HTTP

def _http_get(url: str, retries: int = 4, backoff: float = 2.0) -> bytes:
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError) as e:
            last_err = e
            time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"GET failed after {retries} attempts: {url} :: {last_err}")


def infer_source(accession: str) -> str:
    """Return "ncbi" or "uniprot" based on the accession's shape."""
    return "ncbi" if accession.startswith(_NCBI_PREFIXES) else "uniprot"


def _download(accession: str, source: str) -> bytes:
    if source == "uniprot":
        return _http_get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta")
    if source == "ncbi":
        return _http_get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=protein&id={accession}&rettype=fasta&retmode=text"
        )
    raise ValueError(f"Unknown source {source!r} for {accession}")


# --------------------------------------------------------------- parsing

def _parse_single_fasta(text: str, accession: str, where: str) -> tuple[str, str]:
    """Return (header, sequence). Raise if this is not one usable record."""
    if not text.strip():
        raise ValueError(
            f"{accession}: empty payload from {where}. UniProt returns HTTP 200 "
            f"with an empty body for deleted/demerged entries — check whether "
            f"{accession} still exists."
        )
    headers: list[str] = []
    seq_parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            headers.append(line[1:].strip())
        else:
            seq_parts.append(line.strip())
    if not headers:
        raise ValueError(f"{accession}: no FASTA header in {where}")
    if len(headers) > 1:
        raise ValueError(
            f"{accession}: expected 1 record from {where}, got {len(headers)}. "
            f"An ambiguous accession can expand to several entries; pin a "
            f"specific one."
        )
    sequence = "".join(seq_parts).replace(" ", "")
    if not sequence:
        raise ValueError(f"{accession}: header present but zero-length sequence in {where}")
    bad = set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZUO*")
    if bad:
        raise ValueError(
            f"{accession}: non-amino-acid characters {sorted(bad)!r} in {where}"
        )
    return headers[0], sequence


def _assert_accession_echoed(header: str, accession: str, source: str) -> None:
    """Confirm the record we got back is the accession we asked for.

    UniProt returns ``sp|P13445|RPOS_ECOLI ...``; NCBI returns
    ``WP_003113871.1 ...``. Either way the requested accession must appear in
    the identifier field, otherwise a silent redirect or a server-side
    accession remap has handed us a different protein.
    """
    identifier = header.split()[0]
    core = accession.split(".")[0]
    if core.lower() not in identifier.lower():
        raise ValueError(
            f"{accession}: {source} returned a record whose identifier is "
            f"{identifier!r}, which does not contain the requested accession. "
            f"Refusing to cache a substituted sequence."
        )


def _describe(header: str) -> str:
    """Header minus its leading identifier token."""
    parts = header.split(None, 1)
    return parts[1].strip() if len(parts) > 1 else ""


# --------------------------------------------------------------- public API

def get_protein(
    accession: str,
    source: str | None = None,
    *,
    cache_dir: Path | None = None,
) -> ProteinRecord:
    """Return the protein sequence for ``accession``, fetching only if needed.

    ``source`` is auto-detected from the accession shape when omitted. The
    cached file is validated on read, so a bad file written by an older run is
    rejected rather than trusted.
    """
    source = source or infer_source(accession)
    cache_dir = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE
    target = cache_dir / f"{source}_{accession}.faa"

    if target.exists() and target.stat().st_size > 0:
        text = target.read_text(encoding="utf-8")
        header, sequence = _parse_single_fasta(text, accession, f"cache {target}")
        _assert_accession_echoed(header, accession, source)
        return ProteinRecord(accession, source, header, _describe(header),
                             sequence, target)

    raw = _download(accession, source).decode("utf-8", errors="replace")
    header, sequence = _parse_single_fasta(raw, accession, f"{source} API")
    _assert_accession_echoed(header, accession, source)

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(raw if raw.endswith("\n") else raw + "\n", encoding="utf-8")
    return ProteinRecord(accession, source, header, _describe(header),
                         sequence, target)


# Gene/protein synonyms for reference-label validation.
_SYNONYMS: dict[str, set[str]] = {
    "rpos": {"rpos", "sigma factor rpos", "sigma s", "sigma-38", "katf"},
    "rpod": {"rpod", "sigma factor rpod", "sigma-70", "siga"},
}


def _gene_class(token: str) -> str | None:
    t = token.lower()
    for cls, members in _SYNONYMS.items():
        if t in members:
            return cls
    return None


def check_label(record: ProteinRecord, declared: str) -> LabelCheck:
    """Cross-check a config's human-readable label against the actual record.

    A screen, not a proof: reports ``suspect`` when the organism or the gene
    token in the label is contradicted by the database description. The audit
    CLI prints both strings side by side rather than rewriting the label.
    """
    actual = record.description
    actual_l = actual.lower()

    # Organism: take a capitalised genus from the label, tolerating the
    # "E. coli" / "P. aeruginosa" abbreviations used throughout the configs.
    organism_ok: bool | None = None
    reasons: list[str] = []
    m = re.search(r"\b([A-Z])[a-z]*\.?\s+([a-z]{3,})\b", declared)
    if m:
        initial, species = m.group(1).lower(), m.group(2).lower()
        sp_m = re.search(r"OS=([A-Z][a-z]+)\s+([a-z]+)", actual)
        if sp_m:
            a_genus, a_species = sp_m.group(1).lower(), sp_m.group(2).lower()
            organism_ok = a_genus.startswith(initial) and a_species == species
            if not organism_ok:
                reasons.append(
                    f"organism: label says '{m.group(0)}', record says "
                    f"'{sp_m.group(1)} {sp_m.group(2)}'"
                )
        elif "[" in actual:  # NCBI style: "... [Pseudomonas]"
            bracket = actual[actual.rfind("[") + 1:actual.rfind("]")].lower()
            organism_ok = bracket.startswith(initial)
            if not organism_ok:
                reasons.append(f"organism: label says '{m.group(0)}', record says '{bracket}'")

    # Gene / protein identity.
    gene_ok: bool | None = None
    declared_classes = {
        c for c in (_gene_class(tok) for tok in re.findall(r"[A-Za-z]+", declared)) if c
    }
    if declared_classes:
        actual_tokens = set(re.findall(r"[A-Za-z0-9]+", actual_l))
        gn = re.search(r"GN=(\w+)", actual)
        if gn:
            actual_tokens.add(gn.group(1).lower())
        actual_classes = {c for c in (_gene_class(t) for t in actual_tokens) if c}
        # "RNA polymerase sigma factor RpoS" is caught by token match above;
        # Gene aliases are recognized by the synonym set.
        gene_ok = bool(declared_classes & actual_classes)
        if not gene_ok:
            reasons.append(
                f"gene: label implies {sorted(declared_classes)}, record reads {actual!r}"
            )

    if organism_ok is None and gene_ok is None:
        return LabelCheck(record.accession, declared, actual, None, None,
                          "unverifiable", "no organism or gene token derivable from label")
    verdict = "ok" if (organism_ok is not False and gene_ok is not False) else "suspect"
    return LabelCheck(record.accession, declared, actual, organism_ok, gene_ok,
                      verdict, "; ".join(reasons) or "label consistent with record")
