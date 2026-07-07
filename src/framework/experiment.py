from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from framework.pipeline import run_framework
from framework import module_figures, paper_figures, paper_tables


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ModuleSpec:
    """Filesystem specification for one cancer experiment module."""

    package: str
    cancer_code: str
    config_path: Path
    output_root: Path


MODULES = {
    "PDAC": ModuleSpec(
        package="PDAC",
        cancer_code="PDAC",
        config_path=PROJECT_ROOT / "PDAC" / "config" / "framework_config.json",
        output_root=PROJECT_ROOT / "PDAC" / "result",
    ),
    "LUDA": ModuleSpec(
        package="LUDA",
        cancer_code="LUAD",
        config_path=PROJECT_ROOT / "LUDA" / "config" / "framework_config.json",
        output_root=PROJECT_ROOT / "LUDA" / "result",
    ),
}


def normalize_modules(selection: str) -> list[str]:
    """Normalize a CLI module selector into package keys.

    Parameters
    ----------
    selection:
        ``all``, ``PDAC``, ``LUDA`` or ``LUAD``.

    Returns
    -------
    list[str]
        Ordered package keys to run.
    """
    value = selection.upper()
    if value == "ALL":
        return ["PDAC", "LUDA"]
    if value == "LUAD":
        return ["LUDA"]
    if value in MODULES:
        return [value]
    raise ValueError(f"Unknown module selector: {selection}")


def run_analysis(module_key: str) -> None:
    """Run the framework analysis for one module.

    Parameters
    ----------
    module_key:
        Normalized module key, either ``PDAC`` or ``LUDA``.

    Returns
    -------
    None
        Tables are written below ``<module>/result/tables``.
    """
    spec = MODULES[module_key]
    run_framework(PROJECT_ROOT, config_path=spec.config_path, output_root=spec.output_root)


def run_module_figures(module_key: str) -> None:
    """Regenerate module-level PNG/PDF figures for one module.

    Parameters
    ----------
    module_key:
        ``PDAC`` or ``LUDA``.

    Returns
    -------
    None
        Figures are written below ``<module>/result/figures``.
    """
    spec = MODULES[module_key]
    module_figures.setup_style()
    _, scores, review, features, comparator, perm = module_figures.read_tables(spec.package)
    out_dir = spec.output_root / "figures"
    # These three figures document score space, representative evidence and
    # validation behaviour for the individual module output.
    module_figures.figure_triage(spec.cancer_code, out_dir, scores, review)
    module_figures.figure_evidence(spec.cancer_code, out_dir, features, review)
    module_figures.figure_validation(spec.cancer_code, out_dir, comparator, perm)


def run_paper_outputs() -> None:
    """Regenerate manuscript-level figures, table fragments and consistency log.

    Parameters
    ----------
    None

    Returns
    -------
    None
        Files are written below ``result/paper_figures`` and
        ``result/paper_tables``.
    """
    paper_figures.main()
    paper_tables.main()


def run_experiment(selection: str = "all", skip_analysis: bool = False, skip_figures: bool = False) -> list[str]:
    """Run the complete reproducibility workflow.

    Parameters
    ----------
    selection:
        Module selector passed from the command line.
    skip_analysis:
        If ``True``, existing result tables are reused.
    skip_figures:
        If ``True``, figure and manuscript table regeneration is skipped.

    Returns
    -------
    list[str]
        Normalized module keys that were included in this run.
    """
    modules = normalize_modules(selection)
    for module_key in modules:
        if not skip_analysis:
            run_analysis(module_key)
        if not skip_figures:
            run_module_figures(module_key)

    # Paper-level figures combine both cancers, so regenerate them only when
    # both module result directories are available.
    if not skip_figures and all((MODULES[key].output_root / "tables").exists() for key in MODULES):
        run_paper_outputs()
    return modules

