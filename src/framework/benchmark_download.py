"""Download and verify original articles used for benchmark expansion."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import requests


PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EUROPE_PMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
PMC_FULLTEXT_HTML = "https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/?report=xml"


def sha256(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest of file ``path``."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_xml(session: requests.Session, url: str, path: Path, expected_tag: bytes, force: bool) -> dict:
    """Download XML from ``url`` and return an auditable file-status record.

    Parameters
    ----------
    session:
        Configured HTTP session.
    url:
        Authoritative PubMed or Europe PMC endpoint.
    path:
        Destination path for the original XML artifact.
    expected_tag:
        XML byte sequence required in the response payload.
    force:
        Replace an existing valid file when ``True``.

    Returns
    -------
    dict
        Download status, URL, local path, byte count, and SHA-256 digest.
    """
    if path.exists() and not force:
        payload = path.read_bytes()
        if expected_tag in payload:
            return {
                "status": "exists",
                "url": url,
                "local_path": str(path),
                "bytes": len(payload),
                "sha256": sha256(path),
            }

    response = session.get(url, timeout=120)
    response.raise_for_status()
    payload = response.content
    if expected_tag not in payload:
        raise ValueError(f"Downloaded payload from {url} lacks expected XML tag {expected_tag!r}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {
        "status": "downloaded",
        "url": url,
        "local_path": str(path),
        "bytes": len(payload),
        "sha256": sha256(path),
    }


def _download_pmc_html(session: requests.Session, url: str, path: Path, force: bool) -> dict:
    """Download a PMC author-manuscript HTML page and return its audit record."""
    if path.exists() and not force:
        payload = path.read_bytes()
        lower = payload.lower()
        if b"<html" in lower and b"<article" in lower and b"citation_title" in lower:
            return {
                "status": "exists",
                "url": url,
                "local_path": str(path),
                "bytes": len(payload),
                "sha256": sha256(path),
            }
    response = session.get(url, timeout=120)
    response.raise_for_status()
    payload = response.content
    lower = payload.lower()
    # Challenge/error pages may also be HTML; require article markup and the
    # bibliographic title metadata emitted by a genuine PMC full-text page.
    if b"<html" not in lower or b"<article" not in lower or b"citation_title" not in lower:
        raise ValueError(f"Downloaded payload from {url} is not a PMC full-text article page")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {
        "status": "downloaded",
        "url": url,
        "local_path": str(path),
        "bytes": len(payload),
        "sha256": sha256(path),
    }


def download_expanded_benchmark(module_root: Path, force: bool = False) -> Path:
    """Download all candidate original records listed by one cancer module.

    Parameters
    ----------
    module_root:
        Path to the `PDAC` or `LUDA` package directory.
    force:
        Redownload valid existing files when ``True``.

    Returns
    -------
    pathlib.Path
        Path to the generated article-level download manifest.
    """
    benchmark_root = module_root / "datasources" / "benchmark"
    candidate_path = benchmark_root / "expanded_benchmark_articles.tsv"
    candidates = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    out_dir = benchmark_root / "papers" / "expanded"
    out_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": "evidence-traceable-benchmark/1.0 (academic reproducibility)"})

    records: list[dict[str, object]] = []
    for row in candidates.to_dict("records"):
        key = row["article_key"]
        if row.get("pmid"):
            url = f"{PUBMED_EFETCH}?db=pubmed&id={row['pmid']}&retmode=xml"
            try:
                status = _download_xml(session, url, out_dir / f"{key}_pubmed.xml", b"<PubmedArticle", force)
                records.append({**row, "artifact_type": "pubmed_xml", **status})
            except Exception as exc:
                records.append({**row, "artifact_type": "pubmed_xml", "status": "failed", "url": url, "error": str(exc)})
        if row.get("pmcid"):
            url = EUROPE_PMC_FULLTEXT.format(pmcid=row["pmcid"])
            try:
                status = _download_xml(session, url, out_dir / f"{key}_fulltext.xml", b"<article", force)
                records.append({**row, "artifact_type": "pmc_fulltext_xml", **status})
            except Exception as exc:
                records.append({**row, "artifact_type": "pmc_fulltext_xml", "status": "failed", "url": url, "error": str(exc)})
                fallback_url = PMC_FULLTEXT_HTML.format(pmcid=row["pmcid"])
                try:
                    status = _download_pmc_html(
                        session,
                        fallback_url,
                        out_dir / f"{key}_fulltext.html",
                        force,
                    )
                    records.append({**row, "artifact_type": "pmc_fulltext_html", **status})
                except Exception as fallback_exc:
                    records.append({
                        **row,
                        "artifact_type": "pmc_fulltext_html",
                        "status": "failed",
                        "url": fallback_url,
                        "error": str(fallback_exc),
                    })

    manifest_path = benchmark_root / "papers" / "expanded_article_download_manifest.tsv"
    pd.DataFrame(records).to_csv(manifest_path, sep="\t", index=False)
    return manifest_path

