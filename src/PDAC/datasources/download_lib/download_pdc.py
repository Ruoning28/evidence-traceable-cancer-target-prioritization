"""Download accessible CPTAC-PDAC files from LinkedOmics."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any
from urllib.parse import quote

from PDAC.datasources.download_lib.common import download_file, request_session
from PDAC.datasources.pdac_io import append_missing_data, write_json


BASE_URL = "https://linkedomics.org/data_download/CPTAC-PDAC/"


def _quote_url(url: str) -> str:
    return quote(url, safe=":/?&=%")


def list_linkedomics_files() -> list[str]:
    """Return data file URLs linked from the CPTAC-PDAC download page."""
    session = request_session()
    response = session.get(BASE_URL, timeout=60)
    response.raise_for_status()
    urls = re.findall(r'href="([^"]+)"', response.text)
    return [url for url in urls if url.startswith(BASE_URL)]


def download_cptac_pdac(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download configured LinkedOmics CPTAC-PDAC files."""
    out_dir = project_root / "data" / "raw" / "pdc_cptac_pdac"
    out_dir.mkdir(parents=True, exist_ok=True)
    available = list_linkedomics_files()
    write_json(out_dir / "linkedomics_file_listing.json", available)

    available_by_name = {Path(url).name: url for url in available}
    status: dict[str, Any] = {"source": "LinkedOmics CPTAC-PDAC", "files": []}
    for filename in config["downloads"]["cptac_pdac"]["files"]:
        url = available_by_name.get(filename, BASE_URL + filename)
        if filename not in available_by_name:
            append_missing_data(project_root, f"LinkedOmics CPTAC-PDAC: configured file not found in listing: {filename}")
        try:
            result = download_file(_quote_url(url), out_dir / filename, force=force, timeout=180)
            result["filename"] = filename
        except Exception as exc:
            append_missing_data(project_root, f"LinkedOmics CPTAC-PDAC: failed to download {filename}: {exc}")
            result = {"filename": filename, "status": "failed", "error": str(exc)}
        status["files"].append(result)

    write_json(out_dir / "download_status.json", status)
    return status


