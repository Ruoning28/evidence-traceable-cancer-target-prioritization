"""Download public safety evidence sources."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PDAC.datasources.download_lib.common import download_file
from PDAC.datasources.pdac_io import write_json


def download_safety(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download GTEx, Human Protein Atlas, and gnomAD constraint files."""
    cfg = config["downloads"]["safety"]
    status: dict[str, Any] = {"gtex": [], "hpa": [], "gnomad": []}

    gtex_dir = project_root / "data" / "raw" / "gtex"
    hpa_dir = project_root / "data" / "raw" / "hpa"
    gnomad_dir = project_root / "data" / "raw" / "gnomad"
    for directory in [gtex_dir, hpa_dir, gnomad_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    for item in cfg["gtex_files"]:
        result = download_file(item["url"], gtex_dir / item["filename"], force=force, timeout=300)
        result["filename"] = item["filename"]
        status["gtex"].append(result)

    for item in cfg["hpa_files"]:
        result = download_file(item["url"], hpa_dir / item["filename"], force=force, timeout=300)
        result["filename"] = item["filename"]
        status["hpa"].append(result)

    item = cfg["gnomad_constraint"]
    result = download_file(item["url"], gnomad_dir / item["filename"], force=force, timeout=300)
    result["filename"] = item["filename"]
    status["gnomad"].append(result)

    write_json(gtex_dir / "download_status.json", status)
    write_json(hpa_dir / "download_status.json", status)
    write_json(gnomad_dir / "download_status.json", status)
    return status


