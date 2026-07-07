"""Download supplementary files and create templates for literature benchmarks."""

from __future__ import annotations

import csv
import re
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from PDAC.datasources.download_lib.common import download_file, request_session
from PDAC.datasources.pdac_io import append_missing_data, write_json


def _europe_pmc_bin_url(pmcid: str, filename: str) -> str:
    """Build an NCBI PMC supplementary-file URL."""
    instance_id = pmcid.replace("PMC", "")
    return f"https://pmc.ncbi.nlm.nih.gov/articles/instance/{instance_id}/bin/{filename}"


def _write_manual_template(path: Path, columns: list[str]) -> None:
    """Create an empty manual target-list template if it does not exist."""
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(columns)


def _has_expected_magic(path: Path) -> bool:
    """Reject downloaded challenge/error HTML that was saved as a supplement."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    prefix = path.read_bytes()[:16].lstrip()
    suffix = path.suffix.lower()
    if suffix == ".pdf":
        return prefix.startswith(b"%PDF")
    if suffix == ".xlsx":
        return prefix.startswith(b"PK")
    if suffix in {".csv", ".tsv", ".txt"}:
        return not prefix.lower().startswith((b"<html", b"<!doctype"))
    return not prefix.lower().startswith((b"<html", b"<!doctype"))


def _find_valid_file(source_dir: Path, filename: str) -> Path | None:
    """Find a previously downloaded or extracted valid supplement by filename."""
    for path in source_dir.rglob(filename):
        if _has_expected_magic(path):
            return path
    return None


def _download_pmc_oa_package(project_root: Path, pmcid: str, source_dir: Path, force: bool = False) -> dict[str, Any]:
    """Download and extract the PMC Open Access package when available."""
    api_url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}"
    session = request_session()
    response = session.get(api_url, timeout=60)
    response.raise_for_status()
    root = ET.fromstring(response.text)
    link = root.find(".//link[@format='tgz']")
    if link is None or not link.attrib.get("href"):
        append_missing_data(project_root, f"Literature benchmark {pmcid}: no PMC OA tgz package found.")
        return {"status": "not_available", "api_url": api_url}

    href = link.attrib["href"].replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")
    archive = source_dir / f"{pmcid}_oa_package.tar.gz"
    result = download_file(href, archive, force=force, timeout=300)
    result.update({"url": href, "api_url": api_url})

    extract_dir = source_dir / "oa_package"
    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(extract_dir)
    result["extracted_to"] = str(extract_dir)
    result["status"] = "downloaded" if result["status"] in {"downloaded", "exists"} else result["status"]
    return result


def _safe_name(value: str) -> str:
    """Create a readable filesystem name for article titles."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")[:120]


def _download_figshare_project(project_id: int, source_dir: Path, force: bool = False) -> dict[str, Any]:
    """Download all public files from a figshare project."""
    session = request_session()
    articles_url = f"https://api.figshare.com/v2/projects/{project_id}/articles"
    articles = session.get(articles_url, timeout=60).json()
    figshare_dir = source_dir / "figshare_project"
    figshare_dir.mkdir(parents=True, exist_ok=True)
    status: dict[str, Any] = {"project_id": project_id, "articles_url": articles_url, "articles": []}
    for article in articles:
        detail_url = f"https://api.figshare.com/v2/articles/{article['id']}"
        detail = session.get(detail_url, timeout=60).json()
        article_dir = figshare_dir / f"{article['id']}_{_safe_name(article.get('title', 'article'))}"
        article_dir.mkdir(parents=True, exist_ok=True)
        article_status = {
            "id": article["id"],
            "title": article.get("title", ""),
            "doi": article.get("doi", ""),
            "url_public_html": article.get("url_public_html", ""),
            "files": [],
        }
        for file_info in detail.get("files", []):
            result = download_file(file_info["download_url"], article_dir / file_info["name"], force=force, timeout=600)
            result.update({
                "filename": file_info["name"],
                "download_url": file_info["download_url"],
                "expected_size": file_info.get("size"),
            })
            article_status["files"].append(result)
        status["articles"].append(article_status)
    return status


def download_literature_benchmark(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download article supplementary material and create benchmark target-list templates."""
    out_dir = project_root / "data" / "raw" / "literature_benchmark"
    out_dir.mkdir(parents=True, exist_ok=True)
    status: dict[str, Any] = {"sources": []}

    template_columns = {
        "guo2025": ["gene_symbol", "target_category", "reported_evidence_type", "reported_rank_or_score", "notes"],
        "nwosu2025": ["gene_symbol", "target_category", "reported_evidence_type", "reported_rank_or_score", "notes"],
        "tcga_depmap": ["gene_symbol", "cancer_type_or_context", "reported_score", "notes"],
        "savage2024": ["gene_symbol", "target_category", "reported_score", "notes"],
    }
    manual_files = {
        "guo2025": out_dir / "guo2025_targets.tsv",
        "nwosu2025": out_dir / "nwosu2025_targets.tsv",
        "tcga_depmap": out_dir / "tcga_depmap_targets.tsv",
        "savage2024": out_dir / "savage2024_targets.tsv",
    }
    for key, path in manual_files.items():
        _write_manual_template(path, template_columns[key])

    for source in config["downloads"].get("literature_benchmark", {}).get("sources", []):
        key = source["key"]
        source_dir = out_dir / key
        source_dir.mkdir(parents=True, exist_ok=True)
        source_status: dict[str, Any] = {"key": key, "pmcid": source.get("pmcid"), "files": []}
        if source.get("pmcid"):
            try:
                source_status["oa_package"] = _download_pmc_oa_package(project_root, source["pmcid"], source_dir, force=force)
            except Exception as exc:
                append_missing_data(project_root, f"Literature benchmark {key}: failed to download PMC OA package: {exc}")
                source_status["oa_package"] = {"status": "failed", "error": str(exc)}
        if source.get("direct_files"):
            source_status["direct_files"] = []
            for file_info in source["direct_files"]:
                try:
                    result = download_file(file_info["url"], source_dir / file_info["filename"], force=force, timeout=600)
                    result.update({"filename": file_info["filename"], "url": file_info["url"]})
                    if not _has_expected_magic(Path(result["path"])):
                        raise ValueError("downloaded file is not a valid supplement payload")
                except Exception as exc:
                    append_missing_data(project_root, f"Literature benchmark {key}: failed direct download {file_info['filename']}: {exc}")
                    result = {"filename": file_info["filename"], "url": file_info["url"], "status": "failed", "error": str(exc)}
                source_status["direct_files"].append(result)
        if source.get("figshare_project_id"):
            try:
                source_status["figshare_project"] = _download_figshare_project(source["figshare_project_id"], source_dir, force=force)
            except Exception as exc:
                append_missing_data(project_root, f"Literature benchmark {key}: failed figshare project download: {exc}")
                source_status["figshare_project"] = {"status": "failed", "error": str(exc)}
        for filename in source.get("files", []):
            valid_existing = _find_valid_file(source_dir, filename)
            if valid_existing and not force:
                source_status["files"].append({
                    "filename": filename,
                    "status": "exists",
                    "path": str(valid_existing),
                    "bytes": valid_existing.stat().st_size,
                    "source": "pmc_oa_package",
                })
                continue
            url = _europe_pmc_bin_url(source["pmcid"], filename)
            try:
                result = download_file(url, source_dir / filename, force=force, timeout=180)
                result.update({"filename": filename, "url": url})
                path = Path(result["path"])
                if not _has_expected_magic(path):
                    path.unlink(missing_ok=True)
                    raise ValueError("downloaded file is not a valid supplement payload")
            except Exception as exc:
                append_missing_data(project_root, f"Literature benchmark {key}: failed to download {filename}: {exc}")
                result = {"filename": filename, "url": url, "status": "failed", "error": str(exc)}
            source_status["files"].append(result)
        if source.get("note"):
            source_status["note"] = source["note"]
            append_missing_data(project_root, f"Literature benchmark {key}: {source['note']}")
        status["sources"].append(source_status)

    readme = out_dir / "README.md"
    readme.write_text(
        "This directory stores downloaded article supplementary files and manual target-list templates. "
        "The *_targets.tsv files must be filled or programmatically extracted from the supplementary files before benchmark analysis.\n",
        encoding="utf-8",
    )
    write_json(out_dir / "download_status.json", status)
    return status


