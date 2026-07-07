"""Input/output helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

from .utils import find_project_root


def load_config(project_root: Path | None = None) -> dict[str, Any]:
    """Load the main YAML configuration."""
    root = project_root or find_project_root()
    with (root / "config" / "config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_json(path: Path, payload: Any) -> None:
    """Write JSON with stable formatting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def read_json(path: Path) -> Any:
    """Read JSON from disk."""
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def append_missing_data(project_root: Path, message: str) -> None:
    """Append a line to the missing data report."""
    log_path = project_root / "result" / "logs" / "missing_data_report.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")

