from __future__ import annotations

import json
from pathlib import Path

import yaml

from PDAC.datasources.download_lib.download_depmap import download_depmap
from PDAC.datasources.download_lib.download_druggability import download_druggability
from PDAC.datasources.download_lib.download_external_geo import download_external_geo
from PDAC.datasources.download_lib.download_gdc import download_tcga_paad
from PDAC.datasources.download_lib.download_literature_benchmark import download_literature_benchmark
from PDAC.datasources.download_lib.download_network import download_network
from PDAC.datasources.download_lib.download_pdc import download_cptac_pdac
from PDAC.datasources.download_lib.download_safety import download_safety
from PDAC.datasources.pdac_io import write_json


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PDAC_ROOT = PROJECT_ROOT / "PDAC"
CONFIG_PATH = PDAC_ROOT / "datasources" / "data_sources.yaml"


def download_all(force: bool = False) -> dict[str, object]:
    """Download all public source datasets used to rebuild the PDAC inputs.

    Parameters
    ----------
    force:
        If ``True``, existing files in ``PDAC/data/raw`` are downloaded again.

    Returns
    -------
    dict[str, object]
        Nested download status records keyed by source family.
    """
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    status: dict[str, object] = {}

    # Each source writes its own status file under PDAC/data/raw/<source>.
    status["gdc"] = download_tcga_paad(PDAC_ROOT, config, force=force)
    status["depmap"] = download_depmap(PDAC_ROOT, config, force=force)
    status["cptac"] = download_cptac_pdac(PDAC_ROOT, config, force=force)
    status["geo"] = download_external_geo(PDAC_ROOT, config, force=force)
    status["druggability"] = download_druggability(PDAC_ROOT, config, force=force)
    status["safety"] = download_safety(PDAC_ROOT, config, force=force)
    status["network"] = download_network(PDAC_ROOT, config, force=force)
    status["literature"] = download_literature_benchmark(PDAC_ROOT, config, force=force)

    out_path = PDAC_ROOT / "result" / "logs" / "pdac_download_status.json"
    write_json(out_path, status)
    return status


def main(force: bool = False) -> None:
    """CLI-friendly wrapper for ``download_all``.

    Parameters
    ----------
    force:
        Whether to redownload existing public data files.

    Returns
    -------
    None
        The function writes raw files and a JSON status manifest to disk.
    """
    status = download_all(force=force)
    print(json.dumps({"downloaded_source_groups": list(status)}, indent=2))


if __name__ == "__main__":
    main()

