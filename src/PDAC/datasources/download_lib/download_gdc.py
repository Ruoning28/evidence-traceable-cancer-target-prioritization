"""Download TCGA-PAAD metadata and open data from the GDC API."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from PDAC.datasources.download_lib.common import request_session
from PDAC.datasources.pdac_io import append_missing_data, write_json


GDC_API = "https://api.gdc.cancer.gov"


def _post_json(endpoint: str, payload: dict[str, Any]) -> dict[str, Any]:
    session = request_session()
    response = session.post(f"{GDC_API}{endpoint}", json=payload, timeout=120)
    response.raise_for_status()
    return response.json()


def _filters(project_id: str, extra: list[dict[str, Any]]) -> dict[str, Any]:
    content = [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
        *extra,
    ]
    return {"op": "and", "content": content}


def query_files(project_id: str, extra_filters: list[dict[str, Any]], size: int = 5000) -> list[dict[str, Any]]:
    """Query GDC file records for a project and filter set."""
    payload = {
        "filters": _filters(project_id, extra_filters),
        "fields": "file_id,file_name,data_type,data_format,experimental_strategy,analysis.workflow_type,cases.submitter_id",
        "format": "JSON",
        "size": str(size),
    }
    return _post_json("/files", payload).get("data", {}).get("hits", [])


def write_manifest(path: Path, records: list[dict[str, Any]]) -> None:
    """Write a compact GDC manifest TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["id", "filename", "data_type", "data_format", "workflow_type", "case_submitter_id"])
        for record in records:
            cases = record.get("cases") or [{}]
            workflow = (record.get("analysis") or {}).get("workflow_type", "")
            writer.writerow([
                record.get("file_id", ""),
                record.get("file_name", ""),
                record.get("data_type", ""),
                record.get("data_format", ""),
                workflow,
                cases[0].get("submitter_id", ""),
            ])


def download_gdc_data(file_ids: list[str], out_path: Path, force: bool = False) -> dict[str, Any]:
    """Download multiple GDC file IDs through the GDC data endpoint."""
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        return {"status": "exists", "path": str(out_path), "file_count": len(file_ids)}
    session = request_session()
    response = session.post(f"{GDC_API}/data", json={"ids": file_ids}, stream=True, timeout=120)
    response.raise_for_status()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    with tmp_path.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    tmp_path.replace(out_path)
    return {"status": "downloaded", "path": str(out_path), "file_count": len(file_ids), "bytes": out_path.stat().st_size}


def download_gdc_data_chunked(file_ids: list[str], out_dir: Path, prefix: str, chunk_size: int = 100, force: bool = False) -> dict[str, Any]:
    """Download GDC file IDs in smaller archives to avoid oversized API requests."""
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks = []
    failures = []
    for index in range(0, len(file_ids), chunk_size):
        chunk_ids = file_ids[index:index + chunk_size]
        chunk_number = index // chunk_size + 1
        archive = out_dir / f"{prefix}_chunk_{chunk_number:03d}.tar.gz"
        try:
            result = download_gdc_data(chunk_ids, archive, force=force)
            result["chunk_number"] = chunk_number
            chunks.append(result)
        except Exception as exc:
            failure = {"chunk_number": chunk_number, "file_count": len(chunk_ids), "status": "failed", "error": str(exc)}
            failures.append(failure)
    status = "downloaded" if not failures else "partial"
    if len(failures) == len(range(0, len(file_ids), chunk_size)):
        status = "failed"
    return {"status": status, "chunks": chunks, "failures": failures, "file_count": len(file_ids)}


def download_tcga_paad(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download TCGA-PAAD project summary, clinical records, manifests, and RNA counts."""
    project_id = config["project"]["tcga_project"]
    out_dir = project_root / "data" / "raw" / "gdc_tcga_paad"
    out_dir.mkdir(parents=True, exist_ok=True)
    status: dict[str, Any] = {"source": "GDC", "project_id": project_id}

    session = request_session()
    summary = session.get(f"{GDC_API}/projects/{project_id}", timeout=60)
    summary.raise_for_status()
    write_json(out_dir / "project_summary.json", summary.json())

    cases_payload = {
        "filters": {"op": "in", "content": {"field": "project.project_id", "value": [project_id]}},
        "expand": "diagnoses,demographic,exposures",
        "format": "JSON",
        "size": "2000",
    }
    clinical = _post_json("/cases", cases_payload)
    write_json(out_dir / "clinical_cases.json", clinical)
    status["clinical_cases"] = len(clinical.get("data", {}).get("hits", []))

    expression_records = query_files(
        project_id,
        [
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": [config["downloads"]["gdc"]["expression_workflow"]]}},
            {"op": "in", "content": {"field": "data_format", "value": ["TSV"]}},
        ],
    )
    write_json(out_dir / "gdc_star_counts_files.json", expression_records)
    write_manifest(out_dir / "gdc_star_counts_manifest.tsv", expression_records)
    status["expression_file_count"] = len(expression_records)

    mutation_records = query_files(
        project_id,
        [
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "data_format", "value": ["MAF"]}},
        ],
    )
    write_manifest(out_dir / "gdc_masked_somatic_mutation_manifest.tsv", mutation_records)
    status["mutation_file_count"] = len(mutation_records)

    cnv_records = query_files(
        project_id,
        [
            {"op": "in", "content": {"field": "data_type", "value": ["Copy Number Segment"]}},
        ],
    )
    write_manifest(out_dir / "gdc_copy_number_segment_manifest.tsv", cnv_records)
    status["copy_number_file_count"] = len(cnv_records)

    if config["downloads"]["gdc"].get("download_expression_counts", True) and expression_records:
        max_files = int(config["downloads"]["gdc"].get("expression_max_files", 0) or 0)
        selected = expression_records[:max_files] if max_files > 0 else expression_records
        file_ids = [record["file_id"] for record in selected]
        archive = out_dir / "tcga_paad_star_counts_gdc_download.tar.gz"
        status["expression_download"] = download_gdc_data(file_ids, archive, force=force)
    elif not expression_records:
        append_missing_data(project_root, "GDC TCGA-PAAD STAR Counts: no matching expression files found.")

    if config["downloads"]["gdc"].get("download_mutation_maf", True) and mutation_records:
        archive = out_dir / "tcga_paad_masked_somatic_mutation_gdc_download.tar.gz"
        try:
            status["mutation_download"] = download_gdc_data([record["file_id"] for record in mutation_records], archive, force=force)
        except Exception as exc:
            append_missing_data(project_root, f"GDC TCGA-PAAD mutation download failed: {exc}")
            status["mutation_download"] = {"status": "failed", "error": str(exc)}
    elif not mutation_records:
        append_missing_data(project_root, "GDC TCGA-PAAD mutation: no matching MAF files found.")

    if config["downloads"]["gdc"].get("download_copy_number_segments", True) and cnv_records:
        status["copy_number_download"] = download_gdc_data_chunked(
            [record["file_id"] for record in cnv_records],
            out_dir / "copy_number_segment_chunks",
            "tcga_paad_copy_number_segments",
            chunk_size=75,
            force=force,
        )
        if status["copy_number_download"]["failures"]:
            append_missing_data(project_root, "GDC TCGA-PAAD copy number: one or more chunk downloads failed.")
    elif not cnv_records:
        append_missing_data(project_root, "GDC TCGA-PAAD copy number: no matching segment files found.")

    write_json(out_dir / "download_status.json", status)
    return status


