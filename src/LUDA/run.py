from __future__ import annotations

from framework.experiment import run_analysis, run_module_figures


def main() -> None:
    """Run the LUDA package analysis for the LUAD cancer module.

    Parameters
    ----------
    None

    Returns
    -------
    None
        LUAD result tables and figures are written below ``LUDA/result``.
    """
    run_analysis("LUDA")
    run_module_figures("LUDA")


if __name__ == "__main__":
    main()

