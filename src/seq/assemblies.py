"""NCBI assembly and annotation downloads shared by the sequence workflows."""
from __future__ import annotations

import io
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


USER_AGENT = "HS3-sequence/1.0 (research use; biopython client)"


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


# ---------------------------------------------------------------- assemblies

def _assembly_zip(accession: str, include_gff: bool = False) -> bytes:
    parts = ["GENOME_FASTA", "PROT_FASTA"]
    if include_gff:
        parts.append("GENOME_GFF")    # NCBI Datasets v2 enum name
    qs = "&".join(f"include_annotation_type={p}" for p in parts)
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{accession}/download?{qs}"
    )
    return _http_get(url)


def ensure_assembly(paths, tag: str, accession: str) -> tuple[Path, Path]:
    """Fetch genome FASTA (.fna) and proteome FASTA (.faa) for an assembly.

    Returns (genome_fna_path, proteome_faa_path).
    """
    out_dir = paths.assembly_dir(tag)
    genome = out_dir / "genomic.fna"
    proteome = out_dir / "protein.faa"
    if genome.exists() and proteome.exists():
        return genome, proteome

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [fetch] NCBI Datasets: {accession} ({tag})")
    zip_bytes = _assembly_zip(accession)
    found_genome = found_prot = False
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        for member in zf.namelist():
            lower = member.lower()
            # NCBI Datasets layout: ncbi_dataset/data/<acc>/<acc>_<asm>_genomic.fna
            # and a protein.faa inside the same folder.
            if lower.endswith("_genomic.fna") and not found_genome:
                genome.write_bytes(zf.read(member))
                found_genome = True
            elif lower.endswith("protein.faa") and not found_prot:
                proteome.write_bytes(zf.read(member))
                found_prot = True
    if not (found_genome and found_prot):
        raise RuntimeError(
            f"Could not extract genome+proteome from NCBI Datasets zip "
            f"for {accession}; got genome={found_genome}, proteome={found_prot}"
        )
    return genome, proteome


def ensure_gff(paths, tag: str, accession: str) -> Path:
    """Fetch GFF3 annotation for an assembly. Cache-or-fetch."""
    out_dir = paths.assembly_dir(tag)
    gff = out_dir / "genomic.gff"
    if gff.exists() and gff.stat().st_size > 0:
        return gff

    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [fetch] NCBI Datasets GFF3: {accession} ({tag})")
    zip_bytes = _assembly_zip(accession, include_gff=True)
    found = False
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        for member in zf.namelist():
            if member.lower().endswith(".gff") and "genomic" in member.lower():
                gff.write_bytes(zf.read(member))
                found = True
                break
    if not found:
        raise RuntimeError(f"Could not extract GFF3 from NCBI Datasets zip for {accession}")
    return gff
