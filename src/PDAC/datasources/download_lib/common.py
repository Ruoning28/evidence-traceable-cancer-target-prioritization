"""Common download helpers."""

from __future__ import annotations

import hashlib
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry


USER_AGENT = "pdac-target-prioritization/0.1"


def request_session() -> requests.Session:
    """Create a requests session with a project user agent."""
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1.0,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=None,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def download_file(url: str, out_path: Path, force: bool = False, timeout: int = 60) -> dict:
    """Download a URL to disk with a progress bar and skip existing files."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        return {"status": "exists", "path": str(out_path), "bytes": out_path.stat().st_size}

    session = request_session()
    with session.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length") or 0)
        tmp_path = out_path.with_suffix(out_path.suffix + ".part")
        with tmp_path.open("wb") as handle, tqdm(
            total=total,
            unit="B",
            unit_scale=True,
            desc=out_path.name,
            leave=False,
        ) as bar:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
                    bar.update(len(chunk))
        tmp_path.replace(out_path)
    return {"status": "downloaded", "path": str(out_path), "bytes": out_path.stat().st_size}


def md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute an MD5 checksum for a file."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


