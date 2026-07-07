"""Download public druggability evidence sources."""

from __future__ import annotations

import csv
import gzip
import time
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlencode, urlparse

from PDAC.datasources.download_lib.common import download_file, request_session
from PDAC.datasources.pdac_io import append_missing_data, write_json


OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
DGIDB_BASE = "https://dgidb.org/data/latest/"
CHEMBL_BASE = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
UNIPROT_SEARCH_API = "https://rest.uniprot.org/uniprotkb/search"


def _download_opentargets_disease(disease_id: str, label: str, out_dir: Path, page_size: int = 1000) -> dict[str, Any]:
    """Download all Open Targets associated target rows for a disease."""
    query = """
    query DiseaseAssociatedTargets($diseaseId: String!, $index: Int!, $size: Int!) {
      disease(efoId: $diseaseId) {
        id
        name
        associatedTargets(page: {index: $index, size: $size}) {
          count
          rows {
            score
            target {
              id
              approvedSymbol
              approvedName
            }
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }
    """
    session = request_session()
    rows: list[dict[str, Any]] = []
    disease_name = ""
    total_count = None
    page = 0
    while total_count is None or len(rows) < total_count:
        response = session.post(
            OPENTARGETS_API,
            json={"query": query, "variables": {"diseaseId": disease_id, "index": page, "size": page_size}},
            timeout=120,
        )
        response.raise_for_status()
        payload = response.json()
        disease = payload.get("data", {}).get("disease")
        if not disease:
            raise ValueError(f"Open Targets did not return disease {disease_id}: {payload}")
        disease_name = disease.get("name", "")
        associated = disease["associatedTargets"]
        total_count = associated["count"]
        rows.extend(associated["rows"])
        page += 1
        time.sleep(0.2)
        if not associated["rows"]:
            break

    all_datatypes = sorted({score["id"] for row in rows for score in row.get("datatypeScores", [])})
    out_path = out_dir / f"opentargets_{label}_associated_targets.tsv"
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["disease_id", "disease_name", "target_id", "gene_symbol", "target_name", "overall_score", *all_datatypes]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            datatype_scores = {score["id"]: score["score"] for score in row.get("datatypeScores", [])}
            target = row["target"]
            writer.writerow({
                "disease_id": disease_id,
                "disease_name": disease_name,
                "target_id": target.get("id", ""),
                "gene_symbol": target.get("approvedSymbol", ""),
                "target_name": target.get("approvedName", ""),
                "overall_score": row.get("score", ""),
                **{datatype: datatype_scores.get(datatype, 0.0) for datatype in all_datatypes},
            })
    return {"disease_id": disease_id, "label": label, "rows": len(rows), "path": str(out_path)}


def _download_uniprot_tsv_gz(url: str, out_path: Path, force: bool = False, page_size: int = 500) -> dict[str, Any]:
    """Download UniProt TSV data page by page and gzip it locally."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        return {"status": "exists", "path": str(out_path), "bytes": out_path.stat().st_size}

    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    if tmp_path.exists():
        tmp_path.unlink()

    parsed = urlparse(url)
    params = {key: values[-1] for key, values in parse_qs(parsed.query).items()}
    params["format"] = "tsv"
    params["compressed"] = "false"
    params["size"] = str(page_size)
    params.pop("includeIsoform", None)

    session = request_session()
    next_url = UNIPROT_SEARCH_API + "?" + urlencode(params)
    rows = 0
    wrote_header = False

    with gzip.open(tmp_path, "wt", encoding="utf-8", newline="") as handle:
        while next_url:
            response = session.get(next_url, timeout=180)
            response.raise_for_status()
            lines = response.text.splitlines()
            if lines:
                if not wrote_header:
                    handle.write(lines[0] + "\n")
                    wrote_header = True
                for line in lines[1:]:
                    handle.write(line + "\n")
                    rows += 1
            next_url = None
            for link in response.links.values():
                if link.get("rel") == "next":
                    next_url = link.get("url")
                    break
            time.sleep(0.1)

    tmp_path.replace(out_path)
    return {"status": "downloaded", "path": str(out_path), "bytes": out_path.stat().st_size, "rows": rows}


def download_druggability(project_root: Path, config: dict[str, Any], force: bool = False) -> dict[str, Any]:
    """Download public druggability sources and flag licensed manual sources."""
    cfg = config["downloads"]["druggability"]
    opentargets_dir = project_root / "data" / "raw" / "opentargets"
    chembl_dir = project_root / "data" / "raw" / "chembl"
    dgidb_dir = project_root / "data" / "raw" / "dgidb"
    ttd_dir = project_root / "data" / "raw" / "ttd"
    uniprot_dir = project_root / "data" / "raw" / "uniprot"
    drugcentral_dir = project_root / "data" / "raw" / "drugcentral"
    iuphar_dir = project_root / "data" / "raw" / "iuphar_gtopdb"
    drugbank_dir = project_root / "data" / "raw" / "drugbank_manual"
    for directory in [opentargets_dir, chembl_dir, dgidb_dir, ttd_dir, uniprot_dir, drugcentral_dir, iuphar_dir, drugbank_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    status: dict[str, Any] = {
        "opentargets": [],
        "dgidb": [],
        "chembl": [],
        "ttd": [],
        "uniprot": [],
        "drugcentral": [],
        "iuphar_gtopdb": [],
        "drugbank": "access_denied_or_manual_required",
    }
    for disease in cfg["opentargets_diseases"]:
        out_file = opentargets_dir / f"opentargets_{disease['label']}_associated_targets.tsv"
        if out_file.exists() and out_file.stat().st_size > 0 and not force:
            status["opentargets"].append({"disease_id": disease["id"], "label": disease["label"], "status": "exists", "path": str(out_file)})
        else:
            status["opentargets"].append(_download_opentargets_disease(disease["id"], disease["label"], opentargets_dir))

    for filename in cfg["dgidb_files"]:
        result = download_file(DGIDB_BASE + filename, dgidb_dir / filename, force=force, timeout=180)
        result["filename"] = filename
        status["dgidb"].append(result)

    if cfg.get("download_chembl_sqlite", True):
        for filename in ["chembl_37_release_notes.txt", "chembl_37_sqlite.tar.gz"]:
            result = download_file(CHEMBL_BASE + filename, chembl_dir / filename, force=force, timeout=300)
            result["filename"] = filename
            status["chembl"].append(result)

    for filename in cfg.get("ttd_files", []):
        result = download_file(cfg["ttd_base_url"] + filename, ttd_dir / filename, force=force, timeout=300)
        result["filename"] = filename
        status["ttd"].append(result)

    for item in cfg.get("drugcentral_files", []):
        result = download_file(item["url"], drugcentral_dir / item["filename"], force=force, timeout=600)
        result["filename"] = item["filename"]
        status["drugcentral"].append(result)

    for item in cfg.get("iuphar_gtopdb_files", []):
        result = download_file(item["url"], iuphar_dir / item["filename"], force=force, timeout=600)
        result["filename"] = item["filename"]
        status["iuphar_gtopdb"].append(result)

    uniprot = cfg.get("uniprot_human_proteome")
    if uniprot:
        result = _download_uniprot_tsv_gz(uniprot["url"], uniprot_dir / uniprot["filename"], force=force)
        result["filename"] = uniprot["filename"]
        status["uniprot"].append(result)

    readme = drugbank_dir / "README_manual_required.txt"
    readme.write_text(
        "DrugBank is licensed. Official HTTP Basic Auth release downloads were attempted, "
        "but the current account received HTTP 403 for the full database and downloadable release packages. "
        "If available, place approved DrugBank target/drug export files in this directory.\n",
        encoding="utf-8",
    )
    append_missing_data(project_root, "DrugBank: unavailable because Academic Downloads are closed or this account lacks release download permission.")
    write_json(opentargets_dir / "download_status.json", status)
    write_json(chembl_dir / "download_status.json", status)
    write_json(dgidb_dir / "download_status.json", status)
    write_json(ttd_dir / "download_status.json", status)
    write_json(uniprot_dir / "download_status.json", status)
    write_json(drugcentral_dir / "download_status.json", status)
    write_json(iuphar_dir / "download_status.json", status)
    return status


