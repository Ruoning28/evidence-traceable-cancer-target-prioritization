"""Download external GEO validation cohorts."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PDAC.datasources.download_lib.common import download_file
from PDAC.datasources.pdac_io import write_json


def _geo_series_matrix_url(gse: str) -> str:
    """Build the standard NCBI GEO FTP series-matrix URL for a GSE accession."""
    prefix = f"{gse[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz"


def download_external_geo(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download configured GEO series matrix files."""
    out_dir = project_root / "data" / "raw" / "external_geo"
    out_dir.mkdir(parents=True, exist_ok=True)
    status: dict[str, Any] = {"source": "NCBI GEO", "series": []}
    for gse in config["downloads"]["external_geo"]["series"]:
        url = _geo_series_matrix_url(gse)
        try:
            result = download_file(url, out_dir / f"{gse}_series_matrix.txt.gz", force=force, timeout=180)
            result.update({"gse": gse, "url": url})
        except Exception as exc:
            result = {"gse": gse, "url": url, "status": "failed", "error": str(exc)}
        status["series"].append(result)
    write_json(project_root / "data" / "raw" / "external_geo" / "download_status.json", status)
    return status


