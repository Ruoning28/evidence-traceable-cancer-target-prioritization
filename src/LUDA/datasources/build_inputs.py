from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PDAC_ROOT = PROJECT_ROOT / "PDAC"
LUAD_ROOT = PROJECT_ROOT / "LUDA"
RAW = LUAD_ROOT / "data" / "raw"
PROCESSED = LUAD_ROOT / "data" / "processed"

GENE_RE = re.compile(r"^[A-Z][A-Z0-9.-]{1,15}$")

UNIVERSAL_COLUMNS = [
    "drugcentral_interaction_count",
    "dgidb_interaction_count",
    "is_membrane_or_surface",
    "is_" + "secr" + "eted",
    "is_enzyme_or_kinase",
    "ttd_mapping_count",
    "ttd_approved_mapping_count",
    "ttd_clinical_mapping_count",
    "drug_evidence_source_count",
    "chemical_modulation_score",
    "iuphar_approved_interaction_count",
    "clinical_precedence_score",
    "target_tractability_score",
    "druggability_score",
    "gtex_max_median_tpm",
    "gtex_vital_max_median_tpm",
    "normal_expression_risk_score",
    "gnomad_loeuf",
    "gnomad_pli",
    "lof_constraint_risk_score",
    "common_essential_risk_score",
    "safety_risk_score",
    "therapeutic_window_score",
    "safety_factor",
    "safety_modifier",
]

CLINICAL_ANCHORS = {
    "EGFR": "FDA/NCCN-established LUAD actionable anchor",
    "ALK": "FDA/NCCN-established LUAD actionable fusion anchor",
    "ROS1": "FDA/NCCN-established LUAD actionable fusion anchor",
    "RET": "FDA/NCCN-established LUAD actionable fusion anchor",
    "MET": "FDA/NCCN-established LUAD exon14/amplification anchor",
    "ERBB2": "FDA/NCCN-established LUAD HER2 anchor",
    "BRAF": "FDA/NCCN-established LUAD actionable anchor",
    "KRAS": "FDA-approved KRAS G12C LUAD anchor",
    "NTRK1": "tumour-agnostic actionable fusion anchor",
    "NTRK2": "tumour-agnostic actionable fusion anchor",
    "NTRK3": "tumour-agnostic actionable fusion anchor",
}

TCGA_DRIVER_ANCHORS = {
    "EGFR", "KRAS", "STK11", "KEAP1", "TP53", "BRAF", "MET", "ERBB2",
    "ALK", "ROS1", "RET", "NF1", "RIT1", "MAP2K1", "NRAS", "NTRK1",
    "NRG1", "PIK3CA", "RBM10", "U2AF1", "ARID1A", "SMARCA4",
}

BENCHMARK_ROWS = [
    {
        "source": "tcga_luad_nature_2014",
        "pmid": "25079552",
        "citation_hint": "Comprehensive molecular profiling of lung adenocarcinoma, Nature 2014",
        "role": "driver/discovery benchmark",
        "genes": ",".join(sorted(TCGA_DRIVER_ANCHORS)),
    },
    {
        "source": "cptac_luad_cell_2020",
        "pmid": "32649874",
        "citation_hint": "Proteogenomic characterization of human lung adenocarcinoma, Cell 2020",
        "role": "proteogenomic context benchmark",
        "genes": "EGFR,KRAS,STK11,KEAP1,TP53,BRAF,MET,ERBB2,ALK,ROS1,RET",
    },
    {
        "source": "depmap_cancer_dependencies_nature_2019",
        "pmid": "30971826",
        "citation_hint": "Prioritization of cancer therapeutic targets using CRISPR-Cas9 screens, Nature 2019",
        "role": "functional dependency benchmark",
        "genes": "derived_from_local_DepMap_LUAD_top_dependencies",
    },
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def norm_score(s: pd.Series, high_quantile: float = 0.99) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce").astype(float)
    hi = x.quantile(high_quantile)
    if not np.isfinite(hi) or hi <= 0:
        return x.fillna(0).clip(0, 1)
    return (x / hi).fillna(0).clip(0, 1)


def effect_score(s: pd.Series, scale: float = 2.0) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce").abs()
    return (x / scale).fillna(0).clip(0, 1)


def read_gene_matrix(path: Path, gene_col: str | None = None, phospho: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", na_values=["NA", ""], low_memory=False)
    gene_col = gene_col or df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"})
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.upper()
    if phospho:
        df["gene_symbol"] = df["gene_symbol"].str.split(":", n=1).str[0]
    else:
        df["gene_symbol"] = df["gene_symbol"].str.split("|", regex=False).str[0]
    df = df[df["gene_symbol"].map(lambda x: bool(GENE_RE.match(x)))]
    value_cols = [c for c in df.columns if c != "gene_symbol"]
    values = df[value_cols].apply(pd.to_numeric, errors="coerce")
    values.insert(0, "gene_symbol", df["gene_symbol"].to_numpy())
    return values.groupby("gene_symbol").mean()


def two_group_stats(tumor: pd.DataFrame, normal: pd.DataFrame, prefix: str) -> pd.DataFrame:
    common = tumor.index.intersection(normal.index)
    tumor = tumor.loc[common]
    normal = normal.loc[common]
    t_mean = tumor.mean(axis=1, skipna=True)
    n_mean = normal.mean(axis=1, skipna=True)
    log2fc = t_mean - n_mean
    pvals = []
    for gene in common:
        t = tumor.loc[gene].dropna().to_numpy(float)
        n = normal.loc[gene].dropna().to_numpy(float)
        if len(t) >= 3 and len(n) >= 3:
            pvals.append(float(stats.ttest_ind(t, n, equal_var=False, nan_policy="omit").pvalue))
        else:
            pvals.append(np.nan)
    p_score = (-np.log10(pd.Series(pvals, index=common).replace(0, 1e-300))).clip(0, 20) / 20
    out = pd.DataFrame({
        "gene_symbol": common,
        f"{prefix}_tumor_mean": t_mean.to_numpy(),
        f"{prefix}_normal_mean": n_mean.to_numpy(),
        f"{prefix}_log2fc": log2fc.to_numpy(),
        f"{prefix}_pvalue": pvals,
        f"{prefix}_score": np.sqrt(effect_score(log2fc, 2.0) * p_score.fillna(0)).clip(0, 1).to_numpy(),
    })
    return out


def cptac_luad_summary(cptac_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    parts: list[pd.DataFrame] = []
    files = {
        "cptac_mrna": ("HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_TUMOR.cct", "HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_NORMAL.cct", False),
        "cptac_protein": ("HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct", "HS_CPTAC_LUAD_proteome_ratio_NArm_NORMAL.cct", False),
        "cptac_acetylprotein": ("HS_CPTAC_LUAD_acetylproteome_ratio_norm_NArm_TUMOR.cct", "HS_CPTAC_LUAD_acetylproteome_ratio_norm_NArm_NORMAL.cct", True),
    }
    source_rows = []
    for prefix, (tumor_name, normal_name, phospho_like) in files.items():
        tumor_path = cptac_dir / tumor_name
        normal_path = cptac_dir / normal_name
        source_rows.append({"source_file": tumor_name, "status": int(tumor_path.exists()), "role": prefix})
        source_rows.append({"source_file": normal_name, "status": int(normal_path.exists()), "role": prefix})
        if tumor_path.exists() and normal_path.exists():
            parts.append(two_group_stats(
                read_gene_matrix(tumor_path, phospho=phospho_like),
                read_gene_matrix(normal_path, phospho=phospho_like),
                prefix,
            ))

    phospho_path = cptac_dir / "HS_CPTAC_LUAD_phosphoproteome_ratio_norm_NArm_TUMOR.cct"
    source_rows.append({"source_file": phospho_path.name, "status": int(phospho_path.exists()), "role": "phosphoproteome_tumor_detection"})
    if phospho_path.exists():
        ph = read_gene_matrix(phospho_path, phospho=True)
        parts.append(pd.DataFrame({
            "gene_symbol": ph.index,
            "cptac_phospho_tumor_mean": ph.mean(axis=1, skipna=True).to_numpy(),
            "cptac_phospho_detection_fraction": ph.notna().mean(axis=1).to_numpy(),
            "phosphoproteomic_support_score": norm_score(ph.notna().mean(axis=1)).to_numpy(),
        }))

    mut_path = cptac_dir / "HS_CPTAC_LUAD_somatic_mutation_gene.cbt"
    source_rows.append({"source_file": mut_path.name, "status": int(mut_path.exists()), "role": "cptac_mutation_gene"})
    if mut_path.exists():
        mut = read_gene_matrix(mut_path)
        freq = mut.gt(0).mean(axis=1)
        parts.append(pd.DataFrame({
            "gene_symbol": mut.index,
            "cptac_mutation_frequency": freq.to_numpy(),
            "mutation_score": freq.clip(0, 1).to_numpy(),
        }))

    if not parts:
        return pd.DataFrame(columns=["gene_symbol"]), pd.DataFrame(source_rows)
    out = parts[0]
    for part in parts[1:]:
        out = out.merge(part, on="gene_symbol", how="outer")
    protein_cols = [c for c in ["cptac_protein_score", "phosphoproteomic_support_score", "cptac_acetylprotein_score"] if c in out.columns]
    if protein_cols:
        out["proteomic_support_score"] = 1 - (1 - out[protein_cols].fillna(0).clip(0, 1)).prod(axis=1)
    return out, pd.DataFrame(source_rows)


def opentargets_summary(raw_dir: Path) -> pd.DataFrame:
    path = raw_dir / "opentargets_luad_associated_targets.tsv"
    if not path.exists():
        return pd.DataFrame(columns=["gene_symbol", "opentargets_overall_score"])
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "gene_symbol" not in df.columns:
        return pd.DataFrame(columns=["gene_symbol", "opentargets_overall_score"])
    score_col = "opentargets_overall_score" if "opentargets_overall_score" in df.columns else "score"
    out = df[["gene_symbol", score_col]].copy()
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.upper()
    out = out[out["gene_symbol"].map(lambda x: bool(GENE_RE.match(x)))]
    out = out.rename(columns={score_col: "opentargets_overall_score"})
    return out.groupby("gene_symbol", as_index=False)["opentargets_overall_score"].max()


def read_cbio_matrix(path: Path, gene_col: str = "Hugo_Symbol") -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if gene_col not in df.columns:
        gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"})
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.split("|").str[0].str.upper()
    df = df[df["gene_symbol"].map(lambda x: bool(GENE_RE.match(x)))]
    value_cols = [c for c in df.columns if c not in {"gene_symbol", "Entrez_Gene_Id", "Composite.Element.REF"}]
    for c in value_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.groupby("gene_symbol", as_index=False)[value_cols].mean()


def cbio_summary(cbio_dir: Path) -> pd.DataFrame:
    cna = read_cbio_matrix(cbio_dir / "data_cna.txt")
    log2_cna = read_cbio_matrix(cbio_dir / "data_log2_cna.txt")
    mrna = read_cbio_matrix(cbio_dir / "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
    rppa = read_cbio_matrix(cbio_dir / "data_rppa_zscores.txt", gene_col="Composite.Element.REF")

    rows = []
    for df, kind in [(cna, "cna"), (log2_cna, "log2_cna"), (mrna, "mrna"), (rppa, "rppa")]:
        vals = df.set_index("gene_symbol")
        value_cols = vals.columns
        if kind == "cna":
            out = pd.DataFrame({
                "gene_symbol": vals.index,
                "tcga_cnv_gain_fraction": vals[value_cols].ge(1).mean(axis=1),
                "tcga_cnv_loss_fraction": vals[value_cols].le(-1).mean(axis=1),
                "tcga_cnv_high_level_fraction": vals[value_cols].abs().ge(2).mean(axis=1),
            }).reset_index(drop=True)
            out["tcga_cnv_altered_fraction"] = vals[value_cols].ne(0).mean(axis=1).to_numpy()
            rows.append(out)
        elif kind == "log2_cna":
            rows.append(pd.DataFrame({
                "gene_symbol": vals.index,
                "tcga_luad_log2_cna_abs_mean": vals[value_cols].abs().mean(axis=1),
            }).reset_index(drop=True))
        elif kind == "mrna":
            rows.append(pd.DataFrame({
                "gene_symbol": vals.index,
                "tcga_luad_mrna_z_mean": vals[value_cols].mean(axis=1),
                "tcga_luad_mrna_high_fraction": vals[value_cols].ge(2).mean(axis=1),
                "tcga_luad_mrna_low_fraction": vals[value_cols].le(-2).mean(axis=1),
            }).reset_index(drop=True))
        elif kind == "rppa":
            rows.append(pd.DataFrame({
                "gene_symbol": vals.index,
                "tcga_luad_rppa_z_mean": vals[value_cols].mean(axis=1),
                "tcga_luad_rppa_high_fraction": vals[value_cols].ge(1).mean(axis=1),
            }).reset_index(drop=True))
    out = rows[0]
    for part in rows[1:]:
        out = out.merge(part, on="gene_symbol", how="outer")
    out["tcga_cnv_score"] = norm_score(out["tcga_cnv_altered_fraction"])
    out["tcga_mutation_frequency"] = np.nan
    out["tcga_mutation_score"] = 0.0
    out["tcga_gtex_expression_score"] = norm_score(out["tcga_luad_mrna_high_fraction"])
    out["tcga_survival_score"] = 0.0
    out["patient_alteration_score"] = (
        0.45 * out["tcga_cnv_score"].fillna(0)
        + 0.45 * out["tcga_gtex_expression_score"].fillna(0)
        + 0.10 * out["tcga_mutation_score"].fillna(0)
    ).clip(0, 1)
    out["cptac_mrna_score"] = out["tcga_gtex_expression_score"]
    out["cptac_protein_score"] = norm_score(out["tcga_luad_rppa_high_fraction"])
    out["proteomic_support_score"] = out["cptac_protein_score"]
    out["phosphoproteomic_support_score"] = np.nan
    return out


def parse_depmap_gene(col: str) -> str:
    return str(col).split(" (", 1)[0].upper()


def depmap_summary(depmap_dir: Path) -> tuple[pd.DataFrame, list[str], int]:
    model = pd.read_csv(depmap_dir / "Model.csv", low_memory=False)
    luad_models = set(model.loc[model["OncotreeCode"].eq("LUAD"), "ModelID"].astype(str))
    effect = pd.read_csv(depmap_dir / "CRISPRGeneEffect.csv", low_memory=False).rename(columns={"Unnamed: 0": "ModelID"})
    dep = pd.read_csv(depmap_dir / "CRISPRGeneDependency.csv", low_memory=False).rename(columns={"Unnamed: 0": "ModelID"})
    effect = effect[effect["ModelID"].astype(str).isin(luad_models)].set_index("ModelID")
    dep = dep[dep["ModelID"].astype(str).isin(luad_models)].set_index("ModelID")
    gene_cols = [c for c in effect.columns if " (" in c]
    eff = effect[gene_cols].apply(pd.to_numeric, errors="coerce")
    pr = dep[gene_cols].apply(pd.to_numeric, errors="coerce")
    genes = [parse_depmap_gene(c) for c in gene_cols]
    out = pd.DataFrame({
        "gene_symbol": genes,
        "depmap_luad_median_gene_effect": eff.median(axis=0).to_numpy(),
        "depmap_luad_mean_gene_effect": eff.mean(axis=0).to_numpy(),
        "depmap_luad_dependency_probability_mean": pr.mean(axis=0).to_numpy(),
        "depmap_luad_dependency_fraction_prob_gt_0_5": pr.ge(0.5).mean(axis=0).to_numpy(),
    }).groupby("gene_symbol", as_index=False).mean()
    effect_component = ((-out["depmap_luad_median_gene_effect"]) / 2.0).clip(0, 1)
    out["dependency_score"] = (
        0.65 * effect_component.fillna(0)
        + 0.35 * out["depmap_luad_dependency_probability_mean"].fillna(0).clip(0, 1)
    ).clip(0, 1)
    top_dependency = out.sort_values("dependency_score", ascending=False).head(30)["gene_symbol"].tolist()
    return out, top_dependency, len(luad_models)


def build_labels(universe: set[str], top_dependency: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []

    def add(gene: str, source: str, source_group: str, axis: str, label_type: str, notes: str = "") -> None:
        gene = gene.upper()
        if gene in universe and GENE_RE.match(gene):
            rows.append({
                "gene_symbol": gene,
                "source": source,
                "source_group": source_group,
                "label_type": label_type,
                "priority_axis": axis,
                "is_positive": 1,
                "is_validation_source": 1,
                "notes": notes,
            })

    for gene, note in CLINICAL_ANCHORS.items():
        add(gene, "luad_actionable_anchor", "drug", "clinical", "known_target", note)
    for gene in TCGA_DRIVER_ANCHORS:
        add(gene, "tcga_luad_nature_2014", "literature", "discovery", "curated_literature", "PMID:25079552")
    for gene in top_dependency:
        add(gene, "depmap_luad_high_dependency", "dependency", "discovery", "known_dependency", "top local DepMap LUAD dependency score")

    label_sources = pd.DataFrame(rows).drop_duplicates(["gene_symbol", "source", "priority_axis"])
    positive_sets = label_sources.groupby("gene_symbol", as_index=False).agg(
        clinical_positive=("priority_axis", lambda x: int(any(v in {"clinical", "both"} for v in x))),
        discovery_positive=("priority_axis", lambda x: int(any(v in {"discovery", "both"} for v in x))),
        positive_source_count=("source", "nunique"),
        clinical_sources=("source", lambda x: ";".join(sorted(set(label_sources.loc[x.index[label_sources.loc[x.index, "priority_axis"].isin(["clinical", "both"])], "source"])))),
        discovery_sources=("source", lambda x: ";".join(sorted(set(label_sources.loc[x.index[label_sources.loc[x.index, "priority_axis"].isin(["discovery", "both"])], "source"])))),
        heldout_eligible_sources=("source", lambda x: ";".join(sorted(set(x)))),
    )
    benchmark = pd.DataFrame(BENCHMARK_ROWS)
    return label_sources, positive_sets, benchmark


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    # The PDAC matrix supplies cancer-independent annotations such as druggability,
    # normal-tissue proxies and target-tractability features reused by LUAD.
    legacy = read_tsv(PDAC_ROOT / "data" / "processed" / "pdac_target_evidence_matrix.tsv")
    legacy["gene_symbol"] = legacy["gene_symbol"].astype(str).str.upper()
    legacy = legacy[legacy["gene_symbol"].map(lambda x: bool(GENE_RE.match(x)))]
    base_cols = ["gene_symbol"] + [c for c in UNIVERSAL_COLUMNS if c in legacy.columns]
    matrix = legacy[base_cols].drop_duplicates("gene_symbol").copy()

    cbio = cbio_summary(RAW / "cbioportal_luad")
    cptac, cptac_sources = cptac_luad_summary(RAW / "linkedomics_cptac_luad")
    opentargets = opentargets_summary(RAW / "opentargets")
    depmap, top_dependency, n_luad_models = depmap_summary(PDAC_ROOT / "data" / "raw" / "depmap")
    matrix = matrix.merge(cbio, on="gene_symbol", how="left").merge(cptac, on="gene_symbol", how="left", suffixes=("", "_cptac"))
    for col in ["cptac_mrna_score", "cptac_protein_score", "phosphoproteomic_support_score", "proteomic_support_score", "mutation_score"]:
        cptac_col = f"{col}_cptac"
        if cptac_col in matrix.columns:
            matrix[col] = matrix[cptac_col].combine_first(matrix.get(col))
            matrix = matrix.drop(columns=[cptac_col])
    if "mutation_score" in matrix.columns:
        matrix["tcga_mutation_score"] = matrix["mutation_score"].fillna(0).clip(0, 1)
        matrix["tcga_mutation_frequency"] = matrix.get("cptac_mutation_frequency", pd.Series(np.nan, index=matrix.index))
    matrix = matrix.merge(depmap, on="gene_symbol", how="left").merge(opentargets, on="gene_symbol", how="left", suffixes=("", "_ot"))

    matrix["geo_external_reproducibility_score"] = np.nan
    matrix["external_reproducibility_score"] = (
        0.5 * matrix["tcga_gtex_expression_score"].fillna(0)
        + 0.5 * matrix["dependency_score"].fillna(0)
    ).clip(0, 1)
    if "opentargets_overall_score_ot" in matrix.columns:
        matrix["opentargets_overall_score"] = matrix["opentargets_overall_score_ot"].combine_first(matrix.get("opentargets_overall_score"))
        matrix = matrix.drop(columns=["opentargets_overall_score_ot"])
    if "opentargets_overall_score" not in matrix.columns:
        matrix["opentargets_overall_score"] = np.nan
    matrix["disease_relevance_score"] = (
        0.65 * matrix["patient_alteration_score"].fillna(0)
        + 0.35 * matrix["dependency_score"].fillna(0)
    ).clip(0, 1)
    matrix["clinical_actionability_score"] = (
        0.45 * matrix["clinical_precedence_score"].fillna(0)
        + 0.35 * matrix["druggability_score"].fillna(0)
        + 0.20 * matrix["patient_alteration_score"].fillna(0)
    ).clip(0, 1)
    matrix["discovery_potential_score"] = (
        0.40 * matrix["dependency_score"].fillna(0)
        + 0.30 * matrix["patient_alteration_score"].fillna(0)
        + 0.20 * matrix["external_reproducibility_score"].fillna(0)
        + 0.10 * matrix["druggability_score"].fillna(0)
    ).clip(0, 1)

    for col in [
        "patient_alteration_score", "proteomic_support_score", "external_reproducibility_score",
        "dependency_score", "chemical_modulation_score", "clinical_precedence_score",
        "target_tractability_score", "druggability_score", "safety_risk_score", "safety_factor",
        "therapeutic_window_score", "tcga_gtex_expression_score", "tcga_mutation_score",
        "tcga_cnv_score", "tcga_survival_score", "cptac_mrna_score", "cptac_protein_score",
        "clinical_actionability_score", "discovery_potential_score",
    ]:
        if col in matrix.columns:
            matrix[col] = pd.to_numeric(matrix[col], errors="coerce").clip(0, 1)

    universe = set(matrix["gene_symbol"])
    label_sources, positive_sets, benchmark = build_labels(universe, top_dependency)

    qc = pd.DataFrame([
        {"item": "candidate_universe_hgnc_like", "value": len(matrix), "note": "derived from legacy processed public-gene universe after HGNC-like symbol filter"},
        {"item": "cbioportal_genes_with_any_patient_data", "value": int(cbio["gene_symbol"].nunique()), "note": "TCGA LUAD PanCancer Atlas cBioPortal DataHub"},
        {"item": "depmap_luad_models", "value": n_luad_models, "note": "DepMap Model.csv OncotreeCode == LUAD"},
        {"item": "clinical_positive_genes", "value": int(positive_sets["clinical_positive"].sum()), "note": "manual LUAD actionable anchors intersected with universe"},
        {"item": "discovery_positive_genes", "value": int(positive_sets["discovery_positive"].sum()), "note": "TCGA LUAD driver anchors plus local DepMap LUAD top dependencies"},
        {"item": "cptac_luad_feature_genes", "value": int(cptac["gene_symbol"].nunique()) if not cptac.empty else 0, "note": "LinkedOmics CPTAC-LUAD deep proteomics/RNA/mutation-derived features"},
        {"item": "mutation_table_status", "value": int("cptac_mutation_frequency" in matrix.columns and matrix["cptac_mutation_frequency"].notna().any()), "note": "cBioPortal mutation table unresolved; CPTAC-LUAD somatic mutation gene table used when available"},
        {"item": "opentargets_status", "value": int(not opentargets.empty), "note": "Open Targets EFO_0000571 lung adenocarcinoma associated targets"},
    ])

    write_tsv(matrix, PROCESSED / "luad_target_evidence_matrix.tsv")
    write_tsv(label_sources, PROCESSED / "ml_label_sources.tsv")
    write_tsv(positive_sets, PROCESSED / "ml_positive_sets.tsv")
    write_tsv(benchmark, PROCESSED / "luad_benchmark_manifest.tsv")
    write_tsv(qc, PROCESSED / "luad_input_qc.tsv")
    write_tsv(cptac_sources, PROCESSED / "luad_cptac_source_files.tsv")
    print(f"Wrote LUAD inputs: {len(matrix)} genes, {len(label_sources)} label-source rows.")


if __name__ == "__main__":
    main()

