"""Shared project utilities."""

from __future__ import annotations

from pathlib import Path


def find_project_root(start: Path | None = None) -> Path:
    """Return the project root by walking upward to the config directory."""
    current = Path(start or Path.cwd()).resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "config" / "config.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate project root containing config/config.yaml")


def ensure_directories(paths: list[Path]) -> None:
    """Create each directory in a list if it does not already exist."""
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def human_bytes(num_bytes: int | float | None) -> str:
    """Convert a byte count to a compact human-readable string."""
    if num_bytes is None:
        return "unknown"
    value = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024.0:
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{value:.1f} PB"

