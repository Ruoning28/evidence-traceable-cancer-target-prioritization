"""Download LUAD benchmark-expansion candidate articles."""

from pathlib import Path

from framework.benchmark_download import download_expanded_benchmark


def main(force: bool = False) -> Path:
    """Download LUAD candidate originals and return the manifest path."""
    module_root = Path(__file__).resolve().parents[1]
    return download_expanded_benchmark(module_root, force=force)


if __name__ == "__main__":
    print(main())

