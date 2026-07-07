"""Logging helpers for command-line scripts."""

from __future__ import annotations

import logging
from pathlib import Path


def configure_logging(project_root: Path, name: str) -> logging.Logger:
    """Configure a logger that writes both to console and results/logs."""
    log_dir = project_root / "results" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    file_handler = logging.FileHandler(log_dir / f"{name}.log", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger

