from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
import argparse

from framework.experiment import run_experiment
from framework.verify import verify_results, write_verification_report
from PDAC.datasources.download_data import main as download_pdac_data
from LUDA.datasources.download_data import main as download_luda_data
from LUDA.datasources.build_inputs import main as build_luda_inputs


def parse_args() -> argparse.Namespace:
    """Parse command-line options for the reproducibility workflow.

    Parameters
    ----------
    None

    Returns
    -------
    argparse.Namespace
        Parsed module, download, analysis and verification switches.
    """
    parser = argparse.ArgumentParser(description="Run the evidence-traceable PDAC/LUDA reproducibility bundle.")
    parser.add_argument("--module", default="all", choices=["all", "PDAC", "LUDA", "LUAD"], help="Which module to run.")
    parser.add_argument("--download-data", action="store_true", help="Download public raw data before analysis.")
    parser.add_argument("--force-download", action="store_true", help="Redownload raw files that already exist.")
    parser.add_argument("--build-luda-inputs", action="store_true", help="Rebuild LUDA processed inputs from raw data.")
    parser.add_argument("--skip-analysis", action="store_true", help="Reuse existing result tables.")
    parser.add_argument("--skip-figures", action="store_true", help="Skip module and paper figure regeneration.")
    parser.add_argument("--verify-only", action="store_true", help="Only run the manuscript-number verification checks.")
    return parser.parse_args()


def main() -> int:
    """Run data download, analysis, visualization and verification.

    Parameters
    ----------
    None

    Returns
    -------
    int
        Process exit code. ``0`` means all requested checks passed.
    """
    args = parse_args()
    if args.download_data:
        if args.module in {"all", "PDAC"}:
            download_pdac_data(force=args.force_download)
        if args.module in {"all", "LUDA", "LUAD"}:
            download_luda_data(force=args.force_download)
    if args.build_luda_inputs:
        # LUDA input rebuilding needs LUDA raw files plus PDAC/DepMap raw data.
        build_luda_inputs()
    if not args.verify_only:
        run_experiment(args.module, skip_analysis=args.skip_analysis, skip_figures=args.skip_figures)

    report_path = write_verification_report()
    report = verify_results()
    failures = report[report["status"].ne("PASS")]
    print(f"Wrote verification report: {report_path}")
    if not failures.empty:
        print(failures.to_string(index=False))
        return 1
    print("Verification status: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


