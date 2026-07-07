from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PROJECT_ROOT / "LUDA" / "data" / "raw" / "opentargets"
API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
DISEASE_ID = "EFO_0000571"

QUERY = """
query AssociatedTargets($efoId: String!, $index: Int!, $size: Int!) {
  disease(efoId: $efoId) {
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
        datasourceScores {
          id
          score
        }
      }
    }
  }
}
"""


def post_graphql(variables: dict, retries: int = 5) -> dict:
    """Send one Open Targets GraphQL request with conservative retries.

    Parameters
    ----------
    variables:
        GraphQL variables containing the disease EFO identifier and page
        controls.
    retries:
        Maximum number of attempts for transient server or rate-limit errors.

    Returns
    -------
    dict
        Parsed JSON payload returned by the Open Targets API.
    """
    for attempt in range(retries):
        response = requests.post(
            API_URL,
            json={"query": QUERY, "variables": variables},
            timeout=90,
            headers={"Content-Type": "application/json", "User-Agent": "Codex local reproducibility audit"},
        )
        if response.status_code == 200:
            payload = response.json()
            if "errors" in payload:
                raise RuntimeError(payload["errors"])
            return payload
        if response.status_code in {429, 500, 502, 503, 504}:
            time.sleep(5 * (attempt + 1))
            continue
        response.raise_for_status()
    response.raise_for_status()
    raise RuntimeError("Open Targets query failed after retries")


def main() -> None:
    """Download LUAD-associated target scores into ``LUDA/data/raw/opentargets``."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    page_size = 500
    first = post_graphql({"efoId": DISEASE_ID, "index": 0, "size": page_size})
    disease = first["data"]["disease"]
    assoc = disease["associatedTargets"]
    count = int(assoc["count"])
    rows = []

    def add_rows(payload: dict) -> None:
        for row in payload["data"]["disease"]["associatedTargets"]["rows"]:
            target = row["target"]
            rows.append({
                "disease_id": disease["id"],
                "disease_name": disease["name"],
                "target_id": target["id"],
                "gene_symbol": target["approvedSymbol"],
                "approved_name": target["approvedName"],
                "opentargets_overall_score": row["score"],
                "datasource_scores_json": json.dumps(row.get("datasourceScores") or [], sort_keys=True),
            })

    add_rows(first)
    for page_index in range(1, (count + page_size - 1) // page_size):
        payload = post_graphql({"efoId": DISEASE_ID, "index": page_index, "size": page_size})
        add_rows(payload)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "opentargets_luad_associated_targets.tsv", sep="\t", index=False)
    manifest = pd.DataFrame([{
        "source": "Open Targets Platform GraphQL API",
        "api_url": API_URL,
        "disease_id": DISEASE_ID,
        "disease_name": disease["name"],
        "target_count_reported": count,
        "target_count_downloaded": len(df),
    }])
    manifest.to_csv(OUT_DIR / "download_manifest.tsv", sep="\t", index=False)
    print(f"Downloaded {len(df)} Open Targets associations for {disease['name']} ({DISEASE_ID}).")


if __name__ == "__main__":
    main()

