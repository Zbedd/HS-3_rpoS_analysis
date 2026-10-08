"""Build publication-style organism / gene labels from cached FASTA headers.

Sources (all cached locally; no network):
  * UniProt SwissProt/TrEMBL FASTAs — headers carry ``OS=`` (organism), ``GN=``
    (gene), and the entry's swissprot description.
  * NCBI RefSeq FASTAs (``WP_*``, ``NP_*``) — headers carry the organism in
    square brackets and the product description before the bracket.

A ``Label`` is a small data structure:
  - ``short``    : ``"E. coli rpoS"``      (used as MSA row id)
  - ``display``  : ``"$\\it{E.\\ coli}$ rpoS"`` (matplotlib-ready, italicised
                   binomial nomenclature)
  - ``species``  : ``"Escherichia coli"``  (full binomial for tooltips/manifests)
  - ``gene``     : ``"rpoS"``
  - ``acc``      : the original accession (kept for cross-reference)
  - ``is_hs3``   : ``True`` for HS-3 candidate proteins
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Label:
    short: str        # plain text with spaces (used in matplotlib axis labels)
    msa_id: str       # FASTA-safe id, no whitespace (pyMSAviz truncates at spaces)
    display: str      # matplotlib-ready with italic species
    species: str      # full binomial
    gene: str         # gene symbol (lowercase, e.g. "rpoS")
    acc: str          # source accession
    is_hs3: bool


# Hardcoded fallback mapping for NCBI MULTISPECIES "[Pseudomonas]" entries that
# we know come from PA14 because we pinned the accession in config.yaml.
NCBI_PA14_HINT = {
    "WP_003113871.1": ("Pseudomonas aeruginosa", "rpoS"),
    "WP_003085035.1": ("Pseudomonas aeruginosa", "rpoD"),
}

# Genus → genus-letter mapping for unusual cases (e.g. ambiguous abbreviations).
SPECIES_OVERRIDE: dict[str, str] = {
    # Salmonella typhimurium is conventionally written "S. Typhimurium"
    # (capitalised serovar) — but for our axis labels we use lowercase
    # to match journal convention for the genus species pair.
}


def _abbrev_genus(species: str) -> str:
    """`Escherichia coli` → `E. coli`. Drops strain info in parentheses.

    Salmonella entries collapse to ``S. enterica`` for stylistic uniformity
    (the species name is the taxonomic standard; "Typhimurium" is the
    serovar).
    """
    cleaned = re.sub(r"\s*\(.*?\)\s*$", "", species).strip()
    parts = cleaned.split()
    if not parts:
        return species
    if len(parts) == 1:
        return cleaned
    genus_letter = parts[0][0].upper()
    epithet = parts[1].lower()
    return f"{genus_letter}. {epithet}"


def _format_display(short_species: str, gene: str) -> str:
    """Matplotlib-friendly string with italicised binomial.

    ``E. coli rpoS`` → ``"$\\mathit{E.\\ coli}$ rpoS"``.
    """
    italic = short_species.replace(" ", r"\ ")
    return f"$\\mathit{{{italic}}}$ {gene}"


# ----------------------------------------------------------------- parsers

UNIPROT_OS_RE = re.compile(r"OS=([^=]+?)\s+OX=")
UNIPROT_GN_RE = re.compile(r"GN=(\S+)")
NCBI_BRACKETS_RE = re.compile(r"\[([^\]]+)\]\s*$")


def parse_fasta_header(header: str) -> tuple[str, str] | None:
    """Return (species, gene) from a single FASTA header line.

    Returns ``None`` if the header has neither form we recognise.
    """
    header = header.lstrip(">").strip()
    # UniProt: "sp|P13445|RPOS_ECOLI ... OS=Escherichia coli (...) ... GN=rpoS ..."
    m_os = UNIPROT_OS_RE.search(header)
    m_gn = UNIPROT_GN_RE.search(header)
    if m_os and m_gn:
        return m_os.group(1).strip(), m_gn.group(1).strip()
    # NCBI RefSeq: "WP_... product description [Organism name]"
    m_br = NCBI_BRACKETS_RE.search(header)
    if m_br:
        organism = m_br.group(1).strip()
        # Gene name not in header — caller must supply via NCBI_PA14_HINT or
        # a side channel. Return organism with empty gene to signal that.
        return organism, ""
    return None


def label_from_cached_fasta(
    cached_fasta: Path,
    role_hint: str,
    is_hs3: bool = False,
) -> Label:
    """Read the cached FASTA's first record and produce a Label.

    ``role_hint`` is used as the gene-name fallback when the FASTA header has none (NCBI
    MULTISPECIES entries) and to disambiguate gene synonyms
    when needed.
    """
    head = cached_fasta.read_text(encoding="utf-8").splitlines()[0]
    acc = head.lstrip(">").split()[0].split("|")[1] if head.startswith(">sp|") or head.startswith(">tr|") else head.lstrip(">").split()[0]
    parsed = parse_fasta_header(head)
    species, gene = parsed if parsed else ("?", "")

    # NCBI fallback: look up known PA14 accession → (species, gene)
    if (not gene) and acc in NCBI_PA14_HINT:
        species, gene = NCBI_PA14_HINT[acc]
    if not gene:
        gene = role_hint

    # HS-3 special-case: collapse "Jeongeupia sp. HS-3" + "Jeongeupia sacculi"
    # to "Jeongeupia sacculi" (the published-genome species name)
    if is_hs3:
        species = "Jeongeupia sacculi"

    short_species = _abbrev_genus(species)
    if is_hs3:
        short_text = f"HS-3 {short_species} {gene}"
        msa_id = f"HS3.{short_species.replace(' ', '').replace('.', '')}.{gene}"
        display = f"$\\mathbf{{HS\\text{{-}}3}}$ $\\mathit{{{short_species.replace(' ', chr(92)+' ')}}}$ {gene}"
    else:
        short_text = f"{short_species} {gene}"
        # MSA id: "E.coli.rpoS" — pyMSAviz tokenises at spaces so we use
        # dots / underscores only. Strip the trailing period from the genus
        # abbreviation to keep things compact.
        msa_id = f"{short_species.replace(' ', '').replace('.', '')}.{gene}"
        display = _format_display(short_species, gene)

    return Label(
        short=short_text, msa_id=msa_id, display=display, species=species,
        gene=gene, acc=acc, is_hs3=is_hs3,
    )
