"""Download DepMap public release files."""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

from PDAC.datasources.download_lib.common import download_file, md5_file, request_session
from PDAC.datasources.pdac_io import append_missing_data, write_json


FILE_LISTING_URL = "https://depmap.org/portal/api/download/files"


def _load_file_listing() -> list[dict[str, str]]:
    session = request_session()
    response = session.get(FILE_LISTING_URL, timeout=120)
    response.raise_for_status()
    return list(csv.DictReader(io.StringIO(response.text)))


def _latest_public_release(rows: list[dict[str, str]]) -> str:
    for row in rows:
        if row.get("release", "").lower().startswith("depmap public"):
            return row["release"]
    raise ValueError("Could not identify a DepMap Public release from /download/files")


def download_depmap(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download configured DepMap files from the latest or configured release."""
    out_dir = project_root / "data" / "raw" / "depmap"
    supplemental_dir = out_dir / "supplemental"
    out_dir.mkdir(parents=True, exist_ok=True)
    supplemental_dir.mkdir(parents=True, exist_ok=True)
    rows = _load_file_listing()

    listing_path = out_dir / "depmap_download_files_listing.csv"
    with listing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    release = config["downloads"]["depmap"].get("release", "latest")
    if release == "latest":
        release = _latest_public_release(rows)

    selected_files = config["downloads"]["depmap"]["files"]
    status: dict[str, Any] = {"source": "DepMap", "release": release, "files": [], "supplemental_files": []}
    by_filename = {(row["release"], row["filename"]): row for row in rows}

    for filename in selected_files:
        row = by_filename.get((release, filename))
        if row is None:
            message = f"DepMap {release}: missing configured file {filename}"
            append_missing_data(project_root, message)
            status["files"].append({"filename": filename, "status": "missing"})
            continue
        result = download_file(row["url"], out_dir / filename, force=force, timeout=180)
        result.update({
            "filename": filename,
            "release": release,
            "release_date": row.get("release_date"),
            "expected_md5": row.get("md5_hash"),
        })
        if row.get("md5_hash") and Path(result["path"]).exists():
            result["observed_md5"] = md5_file(Path(result["path"]))
            result["md5_ok"] = result["observed_md5"] == row.get("md5_hash")
        status["files"].append(result)

    for item in config["downloads"]["depmap"].get("supplemental_files", []):
        supplemental_release = item["release"]
        filename = item["filename"]
        row = by_filename.get((supplemental_release, filename))
        if row is None:
            message = f"DepMap {supplemental_release}: missing supplemental file {filename}"
            append_missing_data(project_root, message)
            status["supplemental_files"].append({"filename": filename, "release": supplemental_release, "status": "missing"})
            continue
        out_name = f"{supplemental_release.replace(' ', '_').replace('/', '-')}_{filename}"
        result = download_file(row["url"], supplemental_dir / out_name, force=force, timeout=300)
        result.update({
            "filename": filename,
            "release": supplemental_release,
            "release_date": row.get("release_date"),
            "expected_md5": row.get("md5_hash"),
        })
        if row.get("md5_hash") and Path(result["path"]).exists():
            result["observed_md5"] = md5_file(Path(result["path"]))
            result["md5_ok"] = result["observed_md5"] == row.get("md5_hash")
        status["supplemental_files"].append(result)

    write_json(out_dir / "download_status.json", status)
    return status


