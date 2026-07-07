"""Download pathway and network evidence sources."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PDAC.datasources.download_lib.common import download_file
from PDAC.datasources.pdac_io import write_json


def download_network(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download Reactome and STRING data for pathway/network interpretation."""
    cfg = config["downloads"]["network"]
    reactome_dir = project_root / "data" / "raw" / "reactome"
    string_dir = project_root / "data" / "raw" / "string"
    biogrid_dir = project_root / "data" / "raw" / "biogrid"
    reactome_dir.mkdir(parents=True, exist_ok=True)
    string_dir.mkdir(parents=True, exist_ok=True)
    biogrid_dir.mkdir(parents=True, exist_ok=True)

    status: dict[str, Any] = {"reactome": [], "string": [], "biogrid": []}
    for item in cfg["reactome_files"]:
        result = download_file(item["url"], reactome_dir / item["filename"], force=force, timeout=180)
        result["filename"] = item["filename"]
        status["reactome"].append(result)

    for item in cfg["string_files"]:
        result = download_file(item["url"], string_dir / item["filename"], force=force, timeout=300)
        result["filename"] = item["filename"]
        status["string"].append(result)

    for item in cfg.get("biogrid_files", []):
        result = download_file(item["url"], biogrid_dir / item["filename"], force=force, timeout=300)
        result["filename"] = item["filename"]
        status["biogrid"].append(result)

    write_json(project_root / "data" / "raw" / "reactome" / "download_status.json", status)
    write_json(project_root / "data" / "raw" / "string" / "download_status.json", status)
    write_json(project_root / "data" / "raw" / "biogrid" / "download_status.json", status)
    return status


