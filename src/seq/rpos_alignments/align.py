"""Multiple-sequence alignment via EBI Clustal Omega REST.

submit/poll/fetch flow identical to RQ1's InterProScan client. Result is
cached as Stockholm + aligned FASTA so reruns are offline-reproducible.
"""
from __future__ import annotations

import time
import urllib.error
import os
import urllib.parse
import urllib.request
from pathlib import Path


CLUSTALO_BASE = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
POLL_INTERVAL_SEC = 5
POLL_TIMEOUT_SEC = 1200
USER_AGENT = "HS3-rq2/1.0 (research; biopython)"
# EBI's REST submit endpoint requires a non-empty email field. Default to a
# generic placeholder; users who want failure notifications can set the
# EBI_NOTIFY_EMAIL environment variable.
CONTACT_EMAIL = os.environ.get("EBI_NOTIFY_EMAIL", "noreply@example.com")


def align(
    fasta_text: str, cache_dir: Path, tag: str,
) -> str:
    """Return the Clustal-formatted alignment text (cached)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"{tag}.clustal"
    if cache.exists() and cache.stat().st_size > 0:
        return cache.read_text(encoding="utf-8")

    print(f"  [align] submitting {tag} to EBI Clustal Omega")
    job_id = _submit(fasta_text)
    print(f"  [align]   job {job_id} — polling")
    _wait(job_id)
    aln = _fetch_result(job_id, "aln-clustal_num")
    cache.write_text(aln, encoding="utf-8")
    # also cache the fasta-aln view for downstream parsing
    aln_fa = _fetch_result(job_id, "fa")
    (cache_dir / f"{tag}.aln.fasta").write_text(aln_fa, encoding="utf-8")
    return aln


def aligned_fasta(cache_dir: Path, tag: str) -> Path:
    return cache_dir / f"{tag}.aln.fasta"


def _submit(fasta_text: str) -> str:
    payload = urllib.parse.urlencode({
        "stype": "protein",
        "sequence": fasta_text,
        "email": CONTACT_EMAIL,
        "outfmt": "clustal_num",
    }).encode("utf-8")
    req = urllib.request.Request(
        f"{CLUSTALO_BASE}/run", data=payload,
        headers={"User-Agent": USER_AGENT,
                 "Content-Type": "application/x-www-form-urlencoded",
                 "Accept": "text/plain"},
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.read().decode("utf-8").strip()
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(
            f"Clustal Omega submit failed: HTTP {e.code} {e.reason}\n  body: {body[:500]}"
        )


def _wait(job_id: str) -> None:
    waited = 0
    while waited < POLL_TIMEOUT_SEC:
        req = urllib.request.Request(
            f"{CLUSTALO_BASE}/status/{job_id}",
            headers={"User-Agent": USER_AGENT, "Accept": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                status = r.read().decode("utf-8").strip()
        except urllib.error.URLError:
            time.sleep(POLL_INTERVAL_SEC); waited += POLL_INTERVAL_SEC
            continue
        if status == "FINISHED":
            return
        if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(f"Clustal Omega job {job_id} ended in status {status}")
        time.sleep(POLL_INTERVAL_SEC); waited += POLL_INTERVAL_SEC
    raise RuntimeError(f"Clustal Omega job {job_id} did not finish within {POLL_TIMEOUT_SEC}s")


def _fetch_result(job_id: str, result_type: str) -> str:
    req = urllib.request.Request(
        f"{CLUSTALO_BASE}/result/{job_id}/{result_type}",
        headers={"User-Agent": USER_AGENT, "Accept": "text/plain"},
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8")


def parse_aligned_fasta(path: Path) -> dict[str, str]:
    """Return {id: gapped_sequence}. id is the first token of the header."""
    records: dict[str, str] = {}
    cur: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur is not None:
                records[cur] = "".join(chunks)
            cur = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if cur is not None:
        records[cur] = "".join(chunks)
    return records


def build_fasta_input(entries: list[tuple[str, str]]) -> str:
    """Build a FASTA string from a list of (id, seq) tuples."""
    return "".join(f">{i}\n{s}\n" for i, s in entries)
