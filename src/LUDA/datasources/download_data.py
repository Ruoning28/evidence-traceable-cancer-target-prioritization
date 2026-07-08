from __future__ import annotations

import csv
import time
from pathlib import Path
from urllib.parse import quote

import requests

from LUDA.datasources import download_opentargets


PROJECT_ROOT = Path(__file__).resolve().parents[2]
LUDA_ROOT = PROJECT_ROOT / "LUDA"
RAW = LUDA_ROOT / "data" / "raw"

LINKEDOMICS_URLS = [
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_acetylproteome_ratio_norm_NArm_NORMAL.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_acetylproteome_ratio_norm_NArm_TUMOR.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_phosphoproteome_ratio_norm_NArm_NORMAL.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_phosphoproteome_ratio_norm_NArm_TUMOR.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_proteome_ratio_NArm_NORMAL.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_NORMAL.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_TUMOR.cct",
    "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_somatic_mutation_gene.cbt",
]

CBIO_STUDY = "luad_tcga_pan_can_atlas_2018"
CBIO_FILES = [
    "case_lists/cases_all.txt",
    "case_lists/cases_cnaseq.txt",
    "case_lists/cases_rppa.txt",
    "case_lists/cases_sequenced.txt",
    "data_clinical_patient.txt",
    "data_clinical_sample.txt",
    "data_cna.txt",
    "data_log2_cna.txt",
    "data_mrna_seq_v2_rsem.txt",
    "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
    "data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",
    "data_mutations.txt",
    "data_rppa.txt",
    "data_rppa_zscores.txt",
    "meta_study.txt",
]


def download_file(url: str, out_path: Path, force: bool = False, retries: int = 3) -> dict[str, object]:
    """Download one public file and return a small status record.

    Parameters
    ----------
    url:
        HTTP URL to download.
    out_path:
        Destination path below ``LUDA/data/raw``.
    force:
        Redownload the file even when it already exists.
    retries:
        Number of attempts for transient network failures.

    Returns
    -------
    dict[str, object]
        File name, URL, status and byte count when available.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and not force:
        return {"file": out_path.name, "url": url, "status": "exists", "bytes": out_path.stat().st_size}

    last_error = ""
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=180)
            response.raise_for_status()
            out_path.write_bytes(response.content)
            return {"file": out_path.name, "url": url, "status": "downloaded", "bytes": out_path.stat().st_size}
        except requests.RequestException as exc:
            last_error = str(exc)
            time.sleep(3 * (attempt + 1))
    return {"file": out_path.name, "url": url, "status": "failed", "error": last_error}


def download_linkedomics(force: bool = False) -> list[dict[str, object]]:
    """Download CPTAC-LUAD omics tables from LinkedOmics.

    Parameters
    ----------
    force:
        If ``True``, existing files are downloaded again.

    Returns
    -------
    list[dict[str, object]]
        One status dictionary per requested LinkedOmics file.
    """
    out_dir = RAW / "linkedomics_cptac_luad"
    rows = [download_file(url, out_dir / Path(url).name, force=force) for url in LINKEDOMICS_URLS]
    write_manifest(out_dir / "download_manifest.tsv", rows)
    return rows


def download_cbioportal(force: bool = False) -> list[dict[str, object]]:
    """Download TCGA-LUAD PanCancer Atlas files from cBioPortal DataHub.

    Parameters
    ----------
    force:
        If ``True``, existing files are downloaded again.

    Returns
    -------
    list[dict[str, object]]
        One status dictionary per requested cBioPortal file.
    """
    out_dir = RAW / "cbioportal_luad"
    rows: list[dict[str, object]] = []
    for item in CBIO_FILES:
        encoded = "/".join(quote(part) for part in item.split("/"))
        url = f"https://raw.githubusercontent.com/cBioPortal/datahub/master/public/{CBIO_STUDY}/{encoded}"
        out_name = item.replace("/", "__")
        rows.append(download_file(url, out_dir / out_name, force=force))
    write_manifest(out_dir / "download_manifest.tsv", rows)
    return rows


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a tabular download manifest.

    Parameters
    ----------
    path:
        Destination TSV path.
    rows:
        Status dictionaries produced by one of the download helpers.

    Returns
    -------
    None
        The manifest is written to disk.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def download_all(force: bool = False) -> dict[str, object]:
    """Download LUAD-specific raw sources.

    Parameters
    ----------
    force:
        If ``True``, existing LUAD raw files are downloaded again.

    Returns
    -------
    dict[str, object]
        Status records for cBioPortal, LinkedOmics and Open Targets.
    """
    status = {
        "cbioportal_luad": download_cbioportal(force=force),
        "linkedomics_cptac_luad": download_linkedomics(force=force),
    }
    download_opentargets.main()
    status["opentargets"] = "see LUDA/data/raw/opentargets/download_manifest.tsv"
    return status


def main(force: bool = False) -> None:
    """CLI-friendly wrapper for LUAD raw-data downloads.

    Parameters
    ----------
    force:
        Whether to overwrite existing LUAD raw files.

    Returns
    -------
    None
        Downloaded files and manifests are written below ``LUDA/data/raw``.
    """
    status = download_all(force=force)
    print(list(status))


if __name__ == "__main__":
    main()

