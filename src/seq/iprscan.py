"""Cached EBI InterProScan submissions for protein-domain searches."""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
import os
from dataclasses import dataclass
from pathlib import Path


IPRSCAN_BASE = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
POLL_INTERVAL_SEC = 8
POLL_TIMEOUT_SEC = 1200
USER_AGENT = "HS3-sequence/1.0 (research; biopython)"
# EBI's REST submit endpoints require a non-empty email field for courtesy
# notifications if a job fails. We default to a generic placeholder so no
# personal address is committed; users who want failure notifications can
# set the EBI_NOTIFY_EMAIL environment variable.
CONTACT_EMAIL = os.environ.get("EBI_NOTIFY_EMAIL", "noreply@example.com")


@dataclass
class HmmerHit:
    candidate_id: str   # HS-3 protein accession we asked about
    target: str         # HMM that hit (e.g., PF00669)
    description: str
    evalue: float
    bitscore: float


def scan_candidates(
    candidate_proteins: dict[str, str],
    cache_dir: Path,
    tag: str,
    target_profile: str,
) -> list[HmmerHit]:
    """Submit candidate proteins to InterProScan; return matches to target_profile.

    candidate_proteins: dict of {accession: protein_sequence}
    cache_dir: where to cache the raw InterProScan JSON
    tag: label used in cache filenames
    target_profile: the HMM accession we filter for (e.g., PF00669 for
        a configured domain profile)
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    if not candidate_proteins:
        return []

    cache = cache_dir / f"{tag}__iprscan.json"
    if cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
    else:
        fasta = _format_fasta(candidate_proteins)
        print(f"  [hmmer] submitting {len(candidate_proteins)} candidates from {tag} to InterProScan")
        job_id = _submit(fasta)
        print(f"  [hmmer]   job {job_id} — polling")
        _wait_until_done(job_id)
        data = _fetch_result(job_id)
        cache.write_text(json.dumps(data), encoding="utf-8")

    return _extract_target_hits(data, target_profile)


def _format_fasta(proteins: dict[str, str]) -> str:
    chunks = []
    for acc, seq in proteins.items():
        chunks.append(f">{acc}\n{seq}\n")
    return "".join(chunks)


def _submit(fasta: str) -> str:
    # InterProScan v5 REST: only stype/appl/sequence/goterms/pathways accepted;
    # Query the configured annotation libraries.
    payload = urllib.parse.urlencode({
        "stype": "p",
        "sequence": fasta,
        "appl": "PfamA,NCBIfam",
        "goterms": "false",
        "pathways": "false",
        "email": CONTACT_EMAIL,
    }).encode("utf-8")
    req = urllib.request.Request(
        f"{IPRSCAN_BASE}/run", data=payload,
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
            f"InterProScan submit failed: HTTP {e.code} {e.reason}\n  body: {body[:500]}"
        )


def _wait_until_done(job_id: str) -> None:
    waited = 0
    while waited < POLL_TIMEOUT_SEC:
        req = urllib.request.Request(
            f"{IPRSCAN_BASE}/status/{job_id}",
            headers={"User-Agent": USER_AGENT, "Accept": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                status = r.read().decode("utf-8").strip()
        except urllib.error.URLError as e:
            time.sleep(POLL_INTERVAL_SEC)
            waited += POLL_INTERVAL_SEC
            continue
        if status == "FINISHED":
            return
        if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
            raise RuntimeError(f"InterProScan job {job_id} ended in status {status}")
        time.sleep(POLL_INTERVAL_SEC)
        waited += POLL_INTERVAL_SEC
    raise RuntimeError(f"InterProScan job {job_id} did not finish within {POLL_TIMEOUT_SEC}s")


def _fetch_result(job_id: str) -> dict:
    req = urllib.request.Request(
        f"{IPRSCAN_BASE}/result/{job_id}/json",
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode("utf-8"))


def _extract_target_hits(data: dict, target_profile: str) -> list[HmmerHit]:
    """Find any match in any submitted protein whose signature accession or
    name contains target_profile (case-insensitive)."""
    target = target_profile.upper()
    hits: list[HmmerHit] = []
    for prot in data.get("results", []):
        candidate_id = prot.get("xref", [{}])[0].get("id") or prot.get("xref", [{}])[0].get("name", "?")
        for match in prot.get("matches", []):
            sig = match.get("signature", {})
            acc = (sig.get("accession") or "").upper()
            name = (sig.get("name") or "").upper()
            if target in acc or target in name:
                # InterProScan reports per-location scores; take the best.
                best_e = min(
                    (loc.get("evalue", 1e300) for loc in match.get("locations", [])),
                    default=1e300,
                )
                best_score = max(
                    (loc.get("score", 0.0) for loc in match.get("locations", [])),
                    default=0.0,
                )
                desc = sig.get("description") or sig.get("name") or ""
                hits.append(HmmerHit(
                    candidate_id=candidate_id,
                    target=sig.get("accession", target_profile),
                    description=desc,
                    evalue=float(best_e),
                    bitscore=float(best_score),
                ))
    return hits
