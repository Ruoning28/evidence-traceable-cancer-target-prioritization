from __future__ import annotations

from framework.experiment import run_analysis, run_module_figures


def main() -> None:
    """Run the PDAC analysis and module-level figures.

    Parameters
    ----------
    None

    Returns
    -------
    None
        PDAC result tables and figures are written below ``PDAC/result``.
    """
    run_analysis("PDAC")
    run_module_figures("PDAC")


if __name__ == "__main__":
    main()

