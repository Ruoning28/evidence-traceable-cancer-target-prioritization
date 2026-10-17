"""End-to-end V4 PDAC target triage analysis.

This module implements a practical first complete pass over the downloaded
public datasets. It intentionally keeps raw feature columns beside normalized
scores so every ranking can be audited.
"""

from __future__ import annotations

import math
import re
import gzip
import io
import tarfile
import zipfile
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats

from pdac_tp.io import load_config, write_json


GENE_RE = re.compile(r"^[A-Z0-9][A-Z0-9.-]{1,14}$")
VITAL_GTEX_TERMS = ("Heart", "Brain", "Liver", "Kidney", "Lung", "Pancreas", "Bone_Marrow")
KNOWN_TARGET_DEFAULTS = {
    "clinically_actionable_biomarkers": ["BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "NTRK1", "NTRK2", "NTRK3"],
    "direct_drug_targets": ["EGFR", "ERBB2", "PARP1", "PARP2", "BRAF", "NTRK1", "NTRK2", "NTRK3"],
    "pdac_drivers_hard_to_drug": ["KRAS", "TP53", "SMAD4", "CDKN2A"],
    "synthetic_lethal_or_dependency_targets": ["CHEK1", "CHEK2", "ATR", "WEE1", "PRMT5", "AURKA", "AURKB"],
}


def _clean_gene(value: Any) -> str:
    """Return an HGNC-like gene symbol from common matrix headers."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "n/a", "none", "null", "---", "--"}:
        return ""
    text = re.sub(r"\s+\(\d+\)$", "", text)
    text = text.split("|")[0].strip()
    return text.upper()


def _score_from_effect(effect: pd.Series, cap: float) -> pd.Series:
    """Convert an absolute effect size into a bounded 0-1 score."""
    return effect.abs().fillna(0).clip(0, cap) / cap


def _safe_minmax(series: pd.Series, invert: bool = False) -> pd.Series:
    """Scale a series to 0-1, returning zeros for constant or empty input."""
    values = pd.to_numeric(series, errors="coerce")
    if values.notna().sum() == 0:
        return pd.Series(0.0, index=series.index)
    lo = values.min()
    hi = values.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        scaled = pd.Series(0.0, index=series.index)
    else:
        scaled = (values - lo) / (hi - lo)
    scaled = scaled.fillna(0).clip(0, 1)
    return 1.0 - scaled if invert else scaled


def _geomean_frame(df: pd.DataFrame, columns: list[str], epsilon: float = 1e-6) -> pd.Series:
    """Compute row-wise geometric mean over score columns."""
    arr = df[columns].fillna(0).clip(epsilon, 1.0).to_numpy(float)
    return pd.Series(np.exp(np.mean(np.log(arr), axis=1)), index=df.index)


def _mean_frame(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    """Compute row-wise mean over available bounded score columns."""
    if not columns:
        return pd.Series(0.0, index=df.index)
    return df[columns].fillna(0).clip(0.0, 1.0).mean(axis=1).clip(0.0, 1.0)


def _noisy_or_frame(df: pd.DataFrame, columns: list[str], weights: list[float] | None = None) -> pd.Series:
    """Aggregate substitutable evidence as bounded cumulative support.

    Each component can independently support the same evidence axis. The
    noisy-OR form preserves a strong single signal while allowing multiple
    moderate signals to accumulate without exceeding one.
    """
    if not columns:
        return pd.Series(0.0, index=df.index)
    arr = df[columns].fillna(0).clip(0.0, 1.0).to_numpy(float)
    if weights is None:
        w = np.ones(arr.shape[1], dtype=float)
    else:
        w = np.asarray(weights, dtype=float)
        if len(w) != arr.shape[1]:
            raise ValueError("weights length must match columns length")
    support = 1.0 - np.prod(1.0 - np.clip(arr * w, 0.0, 1.0), axis=1)
    return pd.Series(np.clip(support, 0.0, 1.0), index=df.index)


def _read_gene_matrix(path: Path) -> pd.DataFrame:
    """Read a gene-by-sample table whose first column contains gene symbols."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    gene_col = df.columns[0]
    df[gene_col] = df[gene_col].map(_clean_gene)
    df = df[df[gene_col].astype(bool)].drop_duplicates(gene_col).set_index(gene_col)
    return df.apply(pd.to_numeric, errors="coerce")


def _read_gtex_pancreas(root: Path) -> pd.Series:
    """Read GTEx pancreas median TPM as log2(TPM+1) by gene symbol."""
    gtex = root / "data" / "raw" / "gtex" / "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
    if not gtex.exists():
        return pd.Series(dtype=float)
    gt = pd.read_csv(gtex, sep="\t", skiprows=2, compression="gzip", low_memory=False)
    if "Pancreas" not in gt.columns:
        return pd.Series(dtype=float)
    gt["gene_symbol"] = gt["Description"].map(_clean_gene)
    out = pd.Series(np.log2(pd.to_numeric(gt["Pancreas"], errors="coerce").fillna(0).to_numpy() + 1), index=gt["gene_symbol"].to_numpy())
    return out.groupby(out.index).max()


def _parse_gtf_attrs(attr_text: str) -> dict[str, str]:
    """Parse the semicolon-delimited GTF attributes field."""
    attrs: dict[str, str] = {}
    for part in attr_text.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def _read_gencode_gene_coordinates(root: Path) -> pd.DataFrame:
    """Read protein-coding gene coordinates from a fixed GRCh38 GENCODE GTF."""
    gtf = root / "data" / "raw" / "reference" / "gencode.v36.annotation.gtf.gz"
    if not gtf.exists():
        return pd.DataFrame(columns=["chromosome", "start", "end", "midpoint", "gene_symbol"])
    rows: list[dict[str, Any]] = []
    with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_gtf_attrs(parts[8])
            if attrs.get("gene_type") != "protein_coding":
                continue
            gene = _clean_gene(attrs.get("gene_name"))
            if not gene:
                continue
            chrom = parts[0].replace("chr", "")
            if chrom not in {str(i) for i in range(1, 23)} | {"X", "Y"}:
                continue
            start = int(parts[3])
            end = int(parts[4])
            rows.append({"chromosome": chrom, "start": start, "end": end, "midpoint": (start + end) // 2, "gene_symbol": gene})
    coords = pd.DataFrame(rows).drop_duplicates("gene_symbol")
    coords.to_csv(root / "data" / "processed" / "gdc_tcga_paad" / "gencode_v36_gene_coordinates.tsv", sep="\t", index=False)
    return coords


def _extract_tcga_survival(root: Path) -> pd.DataFrame:
    """Extract overall-survival time and event from GDC clinical cases JSON."""
    path = root / "data" / "raw" / "gdc_tcga_paad" / "clinical_cases.json"
    if not path.exists():
        return pd.DataFrame(columns=["sample_id", "os_days", "os_event"])
    import json

    payload = json.loads(path.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for hit in payload.get("data", {}).get("hits", []):
        case = str(hit.get("submitter_id", "")).strip()
        demo = hit.get("demographic") or {}
        vital = str(demo.get("vital_status", "")).lower()
        event = 1 if vital == "dead" else 0
        days_to_death = pd.to_numeric(demo.get("days_to_death"), errors="coerce")
        followups = []
        for diag in hit.get("diagnoses", []) or []:
            val = pd.to_numeric(diag.get("days_to_last_follow_up"), errors="coerce")
            if pd.notna(val):
                followups.append(float(val))
        os_days = float(days_to_death) if event and pd.notna(days_to_death) else (max(followups) if followups else np.nan)
        if case and np.isfinite(os_days) and os_days > 0:
            rows.append({"sample_id": case, "os_days": os_days, "os_event": event})
    out = pd.DataFrame(rows).drop_duplicates("sample_id")
    out.to_csv(root / "data" / "processed" / "gdc_tcga_paad" / "clinical_survival.tsv", sep="\t", index=False)
    return out


def _logrank_pvalue(time: np.ndarray, event: np.ndarray, group: np.ndarray) -> float:
    """Two-group log-rank p value for median-split survival screening."""
    valid = np.isfinite(time) & np.isfinite(event) & np.isfinite(group)
    time = time[valid].astype(float)
    event = event[valid].astype(int)
    group = group[valid].astype(int)
    if len(time) < 20 or group.sum() < 10 or (len(group) - group.sum()) < 10 or event.sum() < 5:
        return np.nan
    observed = expected = variance = 0.0
    for t in np.sort(np.unique(time[event == 1])):
        at_risk = time >= t
        n = int(at_risk.sum())
        n1 = int((at_risk & (group == 1)).sum())
        d = int(((time == t) & (event == 1)).sum())
        d1 = int(((time == t) & (event == 1) & (group == 1)).sum())
        if n <= 1:
            continue
        expected += d * n1 / n
        variance += n1 * (n - n1) * d * (n - d) / ((n ** 2) * (n - 1))
        observed += d1
    if variance <= 0:
        return np.nan
    chi2 = (observed - expected) ** 2 / variance
    return float(stats.chi2.sf(chi2, 1))


def _build_tcga_survival_features(root: Path, expr: pd.DataFrame) -> pd.DataFrame:
    """Screen expression-survival associations using median-split log-rank tests."""
    clinical = _extract_tcga_survival(root)
    if clinical.empty or expr.empty:
        return pd.DataFrame()
    expr_case = expr.T.groupby(level=0).mean().T
    common_samples = [c for c in expr_case.columns if c in set(clinical["sample_id"])]
    if len(common_samples) < 50:
        return pd.DataFrame()
    clin = clinical.set_index("sample_id").loc[common_samples]
    times = clin["os_days"].to_numpy(float)
    events = clin["os_event"].to_numpy(int)
    rows: list[dict[str, Any]] = []
    for gene, values in expr_case[common_samples].iterrows():
        x = values.to_numpy(float)
        if np.isfinite(x).sum() < 50 or np.nanstd(x) == 0:
            continue
        median = np.nanmedian(x)
        group = (x >= median).astype(int)
        pvalue = _logrank_pvalue(times, events, group)
        if not np.isfinite(pvalue):
            continue
        high = group == 1
        low = group == 0
        high_rate = events[high].sum() / max(times[high].sum(), 1.0)
        low_rate = events[low].sum() / max(times[low].sum(), 1.0)
        hr_proxy = (high_rate + 1e-9) / (low_rate + 1e-9)
        log2_hr = math.log2(hr_proxy)
        p_score = min(-math.log10(max(pvalue, 1e-300)) / 10.0, 1.0)
        direction_score = max(0.0, min(log2_hr / 2.0, 1.0))
        rows.append({
            "gene_symbol": gene,
            "tcga_survival_logrank_pvalue": pvalue,
            "tcga_survival_high_expression_hr_proxy": hr_proxy,
            "tcga_survival_log2_hr_proxy": log2_hr,
            "tcga_survival_score": math.sqrt(p_score * direction_score) if direction_score > 0 else 0.0,
        })
    out = pd.DataFrame(rows).set_index("gene_symbol")
    out.to_csv(root / "data" / "processed" / "gdc_tcga_paad" / "survival_expression_gene_level.tsv", sep="\t")
    return out


def _build_tcga_cnv_features(root: Path) -> pd.DataFrame:
    """Convert TCGA copy-number segments to gene-level alteration features."""
    raw = root / "data" / "raw" / "gdc_tcga_paad"
    processed = root / "data" / "processed" / "gdc_tcga_paad"
    coords = _read_gencode_gene_coordinates(root)
    manifest_path = raw / "gdc_copy_number_segment_manifest.tsv"
    chunk_dir = raw / "copy_number_segment_chunks"
    if coords.empty or not manifest_path.exists() or not chunk_dir.exists():
        return pd.DataFrame()
    manifest = pd.read_csv(manifest_path, sep="\t")
    filename_to_case = {str(row["filename"]).split("/")[-1]: str(row["case_submitter_id"]) for _, row in manifest.iterrows()}
    coords = coords.sort_values(["chromosome", "midpoint"]).reset_index(drop=True)
    genes = coords["gene_symbol"].tolist()
    chrom_maps = {
        chrom: (sub["midpoint"].to_numpy(int), sub.index.to_numpy(int))
        for chrom, sub in coords.groupby("chromosome", sort=False)
    }
    sample_vectors: list[pd.Series] = []
    for archive in sorted(chunk_dir.glob("*.tar.gz")):
        with tarfile.open(archive, "r:gz") as tf:
            for member in tf.getmembers():
                basename = Path(member.name).name
                if not (member.isfile() and basename.endswith(".txt") and basename != "MANIFEST.txt"):
                    continue
                case = filename_to_case.get(basename)
                if not case:
                    continue
                handle = tf.extractfile(member)
                if handle is None:
                    continue
                try:
                    seg = pd.read_csv(handle, sep="\t", low_memory=False)
                except Exception:
                    continue
                needed = {"Chromosome", "Start", "End", "Segment_Mean"}
                if not needed.issubset(seg.columns):
                    continue
                vals = np.full(len(genes), np.nan)
                for _, row in seg.iterrows():
                    chrom = str(row["Chromosome"]).replace("chr", "")
                    if chrom not in chrom_maps:
                        continue
                    start = pd.to_numeric(row["Start"], errors="coerce")
                    end = pd.to_numeric(row["End"], errors="coerce")
                    value = pd.to_numeric(row["Segment_Mean"], errors="coerce")
                    if not (np.isfinite(start) and np.isfinite(end) and np.isfinite(value)):
                        continue
                    mids, idx = chrom_maps[chrom]
                    lo = np.searchsorted(mids, int(start), side="left")
                    hi = np.searchsorted(mids, int(end), side="right")
                    if hi > lo:
                        vals[idx[lo:hi]] = float(value)
                sample_vectors.append(pd.Series(vals, index=genes, name=case))
    if not sample_vectors:
        return pd.DataFrame()
    matrix = pd.concat(sample_vectors, axis=1)
    matrix = matrix.T.groupby(level=0).mean().T
    matrix.to_csv(processed / "copy_number_gene_matrix.tsv", sep="\t")
    gain = (matrix >= 0.2).mean(axis=1)
    loss = (matrix <= -0.2).mean(axis=1)
    high_level = ((matrix >= 0.6) | (matrix <= -0.6)).mean(axis=1)
    out = pd.DataFrame(index=matrix.index)
    out["tcga_cnv_gain_fraction"] = gain
    out["tcga_cnv_loss_fraction"] = loss
    out["tcga_cnv_altered_fraction"] = (gain + loss).clip(0, 1)
    out["tcga_cnv_high_level_fraction"] = high_level
    out["tcga_cnv_score"] = ((out["tcga_cnv_altered_fraction"].clip(0, 0.30) / 0.30) * 0.75 + (out["tcga_cnv_high_level_fraction"].clip(0, 0.10) / 0.10) * 0.25).clip(0, 1)
    out.to_csv(processed / "copy_number_gene_level.tsv", sep="\t")
    return out


def build_tcga_features(root: Path) -> pd.DataFrame:
    """Build TCGA RNA-vs-GTEx, mutation, CNV and survival evidence features."""
    raw = root / "data" / "raw" / "gdc_tcga_paad"
    processed = root / "data" / "processed" / "gdc_tcga_paad"
    processed.mkdir(parents=True, exist_ok=True)
    features: list[pd.DataFrame] = []

    expr_tar = raw / "tcga_paad_star_counts_gdc_download.tar.gz"
    manifest = raw / "gdc_star_counts_manifest.tsv"
    if expr_tar.exists() and manifest.exists():
        mf = pd.read_csv(manifest, sep="\t")
        basename_to_case = {
            str(row["filename"]).split("/")[-1]: str(row.get("case_submitter_id", row.get("cases.0.submitter_id", "")))
            for _, row in mf.iterrows()
        }
        series: list[pd.Series] = []
        sample_rows: list[dict[str, str]] = []
        with tarfile.open(expr_tar, "r:gz") as tf:
            for member in tf.getmembers():
                if not (member.isfile() and member.name.endswith(".tsv")):
                    continue
                basename = Path(member.name).name
                sample_id = basename_to_case.get(basename, basename.replace(".rna_seq.augmented_star_gene_counts.tsv", ""))
                handle = tf.extractfile(member)
                if handle is None:
                    continue
                df = pd.read_csv(handle, sep="\t", comment="#", low_memory=False)
                df = df[df["gene_type"].eq("protein_coding")]
                vals = pd.to_numeric(df["tpm_unstranded"], errors="coerce").fillna(0)
                s = pd.Series(np.log2(vals + 1).to_numpy(), index=df["gene_name"].map(_clean_gene), name=sample_id)
                s = s[s.index.astype(bool)]
                series.append(s.groupby(s.index).max())
                sample_rows.append({"sample_id": sample_id, "source_file": basename})
        if series:
            expr = pd.concat(series, axis=1)
            expr.to_csv(processed / "rna_expression_log2tpm.tsv", sep="\t")
            pd.DataFrame(sample_rows).to_csv(processed / "sample_annotation.tsv", sep="\t", index=False)
            median = expr.median(axis=1, skipna=True)
            gtex_pancreas = _read_gtex_pancreas(root)
            common = median.index.intersection(gtex_pancreas.index)
            tmp = pd.DataFrame(index=common)
            tmp["tcga_paad_median_log2tpm"] = median.loc[common]
            tmp["gtex_pancreas_log2tpm"] = gtex_pancreas.loc[common]
            tmp["tcga_gtex_pancreas_log2fc"] = tmp["tcga_paad_median_log2tpm"] - tmp["gtex_pancreas_log2tpm"]
            tmp["tcga_gtex_expression_score"] = _score_from_effect(tmp["tcga_gtex_pancreas_log2fc"], 3.0)
            tmp.to_csv(processed / "tcga_vs_gtex_pancreas_expression.tsv", sep="\t")
            features.append(tmp)
            survival = _build_tcga_survival_features(root, expr)
            if not survival.empty:
                features.append(survival)

    maf_tar = raw / "tcga_paad_masked_somatic_mutation_gdc_download.tar.gz"
    if maf_tar.exists():
        records: list[pd.DataFrame] = []
        with tarfile.open(maf_tar, "r:gz") as tf:
            for member in tf.getmembers():
                if not (member.isfile() and member.name.endswith(".maf.gz")):
                    continue
                raw_handle = tf.extractfile(member)
                if raw_handle is None:
                    continue
                try:
                    with gzip.open(raw_handle, "rt", encoding="utf-8", errors="replace") as handle:
                        m = pd.read_csv(handle, sep="\t", comment="#", low_memory=False)
                    keep = [c for c in ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"] if c in m.columns]
                    if keep:
                        records.append(m[keep])
                except Exception:
                    continue
        if records:
            mut = pd.concat(records, ignore_index=True)
            mut["gene_symbol"] = mut["Hugo_Symbol"].map(_clean_gene)
            nonsilent_exclude = {"Silent", "Intron", "IGR", "3'UTR", "5'UTR", "RNA"}
            if "Variant_Classification" in mut.columns:
                mut = mut[~mut["Variant_Classification"].isin(nonsilent_exclude)]
            denom = max(1, mut["Tumor_Sample_Barcode"].astype(str).str[:12].nunique() if "Tumor_Sample_Barcode" in mut.columns else 185)
            freq = mut[mut["gene_symbol"].astype(bool)].groupby("gene_symbol").size() / denom
            tmp = pd.DataFrame(index=freq.index)
            tmp["tcga_mutation_frequency"] = freq.clip(0, 1)
            tmp["tcga_mutation_score"] = tmp["tcga_mutation_frequency"].clip(0, 1)
            tmp.to_csv(processed / "mutation_gene_level.tsv", sep="\t")
            features.append(tmp)

    cnv = _build_tcga_cnv_features(root)
    if not cnv.empty:
        features.append(cnv)

    if not features:
        return pd.DataFrame()
    out = pd.concat(features, axis=1)
    tcga_patient_cols = [c for c in ["tcga_gtex_expression_score", "tcga_mutation_score", "tcga_cnv_score", "tcga_survival_score"] if c in out.columns]
    out["tcga_patient_score"] = _noisy_or_frame(out, tcga_patient_cols)
    return out


def _two_group_stats(tumor: pd.DataFrame, normal: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Calculate mean difference, t-test p value and a bounded evidence score."""
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
    out = pd.DataFrame(index=common)
    out[f"{prefix}_tumor_mean"] = t_mean
    out[f"{prefix}_normal_mean"] = n_mean
    out[f"{prefix}_log2fc"] = log2fc
    out[f"{prefix}_pvalue"] = pvals
    p_score = (-np.log10(pd.Series(pvals, index=common).replace(0, 1e-300))).clip(0, 20) / 20
    out[f"{prefix}_score"] = np.sqrt(_score_from_effect(log2fc, 2.0) * p_score.fillna(0)).clip(0, 1)
    return out


def build_patient_and_protein_features(root: Path) -> pd.DataFrame:
    """Build patient expression, proteomic, mutation and CNV features."""
    raw = root / "data" / "raw" / "pdc_cptac_pdac"
    features: list[pd.DataFrame] = []
    mrna_t = raw / "mRNA_RSEM_UQ_log2_Tumor.cct"
    mrna_n = raw / "mRNA_RSEM_UQ_log2_Normal.cct"
    if mrna_t.exists() and mrna_n.exists():
        features.append(_two_group_stats(_read_gene_matrix(mrna_t), _read_gene_matrix(mrna_n), "cptac_mrna"))
    prot_t = raw / "proteomics_gene_level_MD_abundance_tumor.cct"
    prot_n = raw / "proteomics_gene_level_MD_abundance_normal.cct"
    if prot_t.exists() and prot_n.exists():
        features.append(_two_group_stats(_read_gene_matrix(prot_t), _read_gene_matrix(prot_n), "cptac_protein"))
    phospho_t = raw / "phosphoproteomics_gene_level_MD_abundance_tumor.cct"
    if phospho_t.exists():
        ph = _read_gene_matrix(phospho_t)
        tmp = pd.DataFrame(index=ph.index)
        tmp["cptac_phospho_tumor_mean"] = ph.mean(axis=1, skipna=True)
        tmp["cptac_phospho_detection_fraction"] = ph.notna().mean(axis=1)
        tmp["phosphoproteomic_support_score"] = _safe_minmax(tmp["cptac_phospho_detection_fraction"])
        features.append(tmp)
    scna = raw / "SCNA_log2_gene_level.cct"
    if scna.exists():
        cn = _read_gene_matrix(scna)
        tmp = pd.DataFrame(index=cn.index)
        tmp["cptac_cnv_abs_mean"] = cn.abs().mean(axis=1, skipna=True)
        tmp["cptac_cnv_altered_fraction"] = (cn.abs() > 0.3).mean(axis=1)
        tmp["cnv_score"] = tmp["cptac_cnv_altered_fraction"].clip(0, 1)
        features.append(tmp)
    maf = raw / "Mutation_results.maf"
    if maf.exists():
        m = pd.read_csv(maf, sep="\t", low_memory=False)
        gene_col = "Hugo_Symbol" if "Hugo_Symbol" in m.columns else m.columns[0]
        case_col = next((c for c in ["Tumor_Sample_Barcode", "case_id", "Sample"] if c in m.columns), None)
        m[gene_col] = m[gene_col].map(_clean_gene)
        denom = max(1, m[case_col].nunique() if case_col else 140)
        freq = m[m[gene_col].astype(bool)].groupby(gene_col).size() / denom
        tmp = pd.DataFrame(index=freq.index)
        tmp["cptac_mutation_frequency"] = freq.clip(0, 1)
        tmp["mutation_score"] = tmp["cptac_mutation_frequency"].clip(0, 1)
        features.append(tmp)
    if not features:
        return pd.DataFrame()
    out = pd.concat(features, axis=1)
    patient_cols = [c for c in ["cptac_mrna_score", "cnv_score", "mutation_score"] if c in out.columns]
    out["patient_alteration_score"] = _noisy_or_frame(out, patient_cols)
    protein_cols = [c for c in ["cptac_protein_score", "phosphoproteomic_support_score"] if c in out.columns]
    out["proteomic_support_score"] = _noisy_or_frame(out, protein_cols) if protein_cols else 0.0
    return out


def build_dependency_features(root: Path) -> pd.DataFrame:
    """Build DepMap PDAC dependency and specificity features."""
    dep = root / "data" / "raw" / "depmap"
    model_path = dep / "Model.csv"
    if not model_path.exists():
        return pd.DataFrame()
    models = pd.read_csv(model_path, low_memory=False)
    text_cols = [c for c in models.columns if any(t in c.lower() for t in ["lineage", "disease", "subtype", "oncotree"])]
    search = models[text_cols].astype(str).agg(" ".join, axis=1).str.lower()
    pdac_ids = set(models.loc[search.str.contains("pancre", na=False), "ModelID"])
    if not pdac_ids:
        return pd.DataFrame()
    effect = pd.read_csv(dep / "CRISPRGeneEffect.csv", low_memory=False).rename(columns={"Unnamed: 0": "ModelID"})
    prob = pd.read_csv(dep / "CRISPRGeneDependency.csv", low_memory=False).rename(columns={"Unnamed: 0": "ModelID"})
    pdac_effect = effect[effect["ModelID"].isin(pdac_ids)].set_index("ModelID")
    other_effect = effect[~effect["ModelID"].isin(pdac_ids)].set_index("ModelID")
    pdac_prob = prob[prob["ModelID"].isin(pdac_ids)].set_index("ModelID")
    genes = [_clean_gene(c) for c in pdac_effect.columns]
    pdac_effect.columns = genes
    other_effect.columns = genes
    pdac_prob.columns = [_clean_gene(c) for c in pdac_prob.columns]
    pdac_effect = pdac_effect.loc[:, pdac_effect.columns.astype(bool)]
    other_effect = other_effect.loc[:, other_effect.columns.astype(bool)]
    pdac_prob = pdac_prob.loc[:, pdac_prob.columns.astype(bool)]
    out = pd.DataFrame(index=pdac_effect.columns)
    out["depmap_pdac_median_gene_effect"] = pdac_effect.median(axis=0, skipna=True)
    out["depmap_pdac_mean_gene_effect"] = pdac_effect.mean(axis=0, skipna=True)
    out["depmap_other_median_gene_effect"] = other_effect.median(axis=0, skipna=True)
    out["depmap_pdac_specificity_delta"] = out["depmap_other_median_gene_effect"] - out["depmap_pdac_median_gene_effect"]
    out["depmap_pdac_dependency_probability_mean"] = pdac_prob.mean(axis=0, skipna=True)
    out["depmap_pdac_dependency_fraction_prob_gt_0_5"] = (pdac_prob > 0.5).mean(axis=0)
    common_path = dep / "AchillesCommonEssentialControls.csv"
    common = set()
    if common_path.exists():
        common = set(pd.read_csv(common_path)["Gene"].map(_clean_gene))
    out["is_common_essential"] = out.index.to_series().isin(common).astype(int)
    effect_score = (-out["depmap_pdac_median_gene_effect"]).clip(0, 2) / 2
    prob_score = out["depmap_pdac_dependency_probability_mean"].fillna(0).clip(0, 1)
    frac_score = out["depmap_pdac_dependency_fraction_prob_gt_0_5"].fillna(0).clip(0, 1)
    specificity_score = out["depmap_pdac_specificity_delta"].clip(0, 1.5) / 1.5
    dep_components = pd.DataFrame({
        "effect_score": effect_score,
        "probability_score": prob_score,
        "fraction_score": frac_score,
        "specificity_score": specificity_score,
    }, index=out.index)
    out["dependency_score"] = _geomean_frame(dep_components, list(dep_components.columns)).clip(0, 1)
    return out.groupby(out.index).max()


def build_druggability_features(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build public druggability features and target-drug edges."""
    raw = root / "data" / "raw"
    frames: list[pd.DataFrame] = []
    edges: list[pd.DataFrame] = []
    ot_files = list((raw / "opentargets").glob("*.tsv"))
    if ot_files:
        ot = pd.concat([pd.read_csv(p, sep="\t", low_memory=False) for p in ot_files], ignore_index=True)
        ot["gene_symbol"] = ot["gene_symbol"].map(_clean_gene)
        agg = ot.groupby("gene_symbol").agg(
            opentargets_overall_score=("overall_score", "max"),
            opentargets_clinical_score=("clinical", "max"),
            opentargets_literature_score=("literature", "max"),
            opentargets_somatic_mutation_score=("somatic_mutation", "max"),
        )
        frames.append(agg)
    dc_path = raw / "drugcentral" / "drug.target.interaction.tsv.gz"
    if dc_path.exists():
        dc = pd.read_csv(dc_path, sep="\t", compression="gzip", low_memory=False)
        dc["gene_symbol"] = dc["GENE"].map(_clean_gene)
        dc = dc[dc["gene_symbol"].astype(bool)]
        frames.append(pd.DataFrame({"drugcentral_interaction_count": dc.groupby("gene_symbol").size()}))
        edges.append(dc[["gene_symbol", "DRUG_NAME"]].dropna().drop_duplicates().rename(columns={"DRUG_NAME": "drug_name"}).assign(source="DrugCentral"))
    dgidb_path = raw / "dgidb" / "interactions.tsv"
    if dgidb_path.exists():
        dg = pd.read_csv(dgidb_path, sep="\t", low_memory=False)
        gene_col = "gene_claim_name" if "gene_claim_name" in dg.columns else dg.columns[0]
        drug_col = "drug_claim_name" if "drug_claim_name" in dg.columns else dg.columns[-1]
        dg["gene_symbol"] = dg[gene_col].map(_clean_gene)
        dg = dg[dg["gene_symbol"].astype(bool)]
        frames.append(pd.DataFrame({"dgidb_interaction_count": dg.groupby("gene_symbol").size()}))
        edges.append(dg[["gene_symbol", drug_col]].dropna().drop_duplicates().rename(columns={drug_col: "drug_name"}).assign(source="DGIdb"))
    iuphar_dir = raw / "iuphar_gtopdb"
    map_path = iuphar_dir / "GtP_to_HGNC_mapping.tsv"
    approved_path = iuphar_dir / "approved_drug_primary_target_interactions.tsv"
    if map_path.exists():
        mp = pd.read_csv(map_path, sep="\t", comment="#", low_memory=False)
        id_col = next((c for c in mp.columns if c.lower().replace(" ", "") in {"iupharid", "targetid"}), None)
        gene_col = next((c for c in mp.columns if "hgnc symbol" in c.lower()), None)
        if id_col and gene_col:
            mp["gene_symbol"] = mp[gene_col].map(_clean_gene)
            if approved_path.exists():
                ap = pd.read_csv(approved_path, sep="\t", comment="#", low_memory=False)
                ap = ap.rename(columns={next((c for c in ap.columns if c.lower() == "target id"), "Target ID"): "Target ID"})
                merged = ap.merge(mp[[id_col, "gene_symbol"]], left_on="Target ID", right_on=id_col, how="left")
                merged = merged[merged["gene_symbol"].astype(str).astype(bool)]
                frames.append(pd.DataFrame({"iuphar_approved_interaction_count": merged.groupby("gene_symbol").size()}))
                drug_col = "Ligand" if "Ligand" in merged.columns else merged.columns[0]
                edges.append(merged[["gene_symbol", drug_col]].dropna().drop_duplicates().rename(columns={drug_col: "drug_name"}).assign(source="IUPHAR_GtoPdb"))
            frames.append(pd.DataFrame({"iuphar_target_family_member": mp.groupby("gene_symbol").size()}))
    hpa_zip = raw / "hpa" / "proteinatlas.tsv.zip"
    if hpa_zip.exists():
        with zipfile.ZipFile(hpa_zip) as zf:
            with zf.open(zf.namelist()[0]) as handle:
                hpa = pd.read_csv(handle, sep="\t", low_memory=False, usecols=lambda c: c in ["Gene", "Protein class", "Subcellular location", "Secretome location"])
        hpa["gene_symbol"] = hpa["Gene"].map(_clean_gene)
        text = hpa.fillna("").astype(str).agg(" ".join, axis=1).str.lower()
        hpa_feat = pd.DataFrame(index=hpa["gene_symbol"])
        hpa_feat["is_membrane_or_surface"] = text.str.contains("membrane|cell surface|plasma membrane").astype(int).to_numpy()
        hpa_feat["is_secreted"] = text.str.contains("secreted|extracellular|secretome").astype(int).to_numpy()
        hpa_feat["is_enzyme_or_kinase"] = text.str.contains("enzyme|kinase|protease|transferase|phosphatase").astype(int).to_numpy()
        frames.append(hpa_feat.groupby(hpa_feat.index).max())
    ttd_dir = raw / "ttd"
    target_path = ttd_dir / "P1-01-TTD_target_download.txt"
    map_path = ttd_dir / "P1-07-Drug-TargetMapping.xlsx"
    drug_path = ttd_dir / "P1-02-TTD_drug_download.txt"
    if target_path.exists() and map_path.exists():
        target_gene: dict[str, str] = {}
        current = ""
        with target_path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                if parts[1] == "TARGETID":
                    current = parts[2]
                elif current and parts[1] == "GENENAME":
                    gene = _clean_gene(parts[2].split(";")[0].split(",")[0])
                    if gene:
                        target_gene[current] = gene
        drug_names: dict[str, str] = {}
        if drug_path.exists():
            current_drug = ""
            with drug_path.open("r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    if parts[1] == "DRUG__ID":
                        current_drug = parts[2]
                    elif current_drug and parts[1] == "TRADNAME":
                        drug_names[current_drug] = parts[2]
        try:
            tm = pd.read_excel(map_path, engine="openpyxl")
            tm["gene_symbol"] = tm["TargetID"].map(target_gene).map(_clean_gene)
            tm = tm[tm["gene_symbol"].astype(bool)]
            tm["drug_name"] = tm["DrugID"].map(drug_names).fillna(tm["DrugID"])
            status = tm["Highest_status"].astype(str).str.lower()
            approved = status.str.contains("approved")
            clinical = status.str.contains("approved|phase|clinical")
            frames.append(pd.DataFrame({
                "ttd_mapping_count": tm.groupby("gene_symbol").size(),
                "ttd_approved_mapping_count": tm[approved].groupby("gene_symbol").size(),
                "ttd_clinical_mapping_count": tm[clinical].groupby("gene_symbol").size(),
            }).fillna(0))
            edges.append(tm[["gene_symbol", "drug_name"]].drop_duplicates().assign(source="TTD"))
        except Exception as exc:
            raise RuntimeError(
                f"Could not read TTD mapping workbook {map_path}. "
                "Install the dependencies in requirements.txt and verify the workbook format and columns."
            ) from exc
    if not frames:
        return pd.DataFrame(), pd.DataFrame()
    out = pd.concat(frames, axis=1).fillna(0)
    source_cols = [c for c in out.columns if c.endswith("_count") or c.endswith("_score") or c.endswith("_member")]
    out["drug_evidence_source_count"] = sum((out[c] > 0).astype(int) for c in source_cols)
    out["chemical_modulation_score"] = (out["drug_evidence_source_count"] / 5).clip(0, 1)
    for col in ["opentargets_clinical_score", "drugcentral_interaction_count", "iuphar_approved_interaction_count", "dgidb_interaction_count", "ttd_approved_mapping_count", "ttd_clinical_mapping_count"]:
        if col not in out.columns:
            out[col] = 0
    clinical_precedence_cols = {
        "opentargets_clinical_indicator": (out["opentargets_clinical_score"] > 0).astype(float),
        "drugcentral_indicator": (out["drugcentral_interaction_count"] > 0).astype(float),
        "iuphar_approved_indicator": (out["iuphar_approved_interaction_count"] > 0).astype(float),
        "ttd_approved_indicator": (out["ttd_approved_mapping_count"] > 0).astype(float),
        "dgidb_indicator": (out["dgidb_interaction_count"] > 0).astype(float),
    }
    out["clinical_precedence_score"] = pd.DataFrame(clinical_precedence_cols, index=out.index).mean(axis=1).clip(0, 1)
    tract_cols = [c for c in ["is_membrane_or_surface", "is_secreted", "is_enzyme_or_kinase", "iuphar_target_family_member"] if c in out.columns]
    out["target_tractability_score"] = out[tract_cols].gt(0).mean(axis=1) if tract_cols else 0.0
    out["druggability_score"] = _mean_frame(out.assign(
        chemical_modulation_score=out["chemical_modulation_score"],
        clinical_precedence_score=out["clinical_precedence_score"],
        target_tractability_score=out["target_tractability_score"],
    ), ["chemical_modulation_score", "clinical_precedence_score", "target_tractability_score"])
    edge_df = pd.concat(edges, ignore_index=True) if edges else pd.DataFrame(columns=["gene_symbol", "drug_name", "source"])
    return out.groupby(out.index).max(), edge_df


def build_safety_features(root: Path, common_essential: pd.Series | None = None) -> pd.DataFrame:
    """Build normal tissue, gnomAD and common essentiality safety features."""
    raw = root / "data" / "raw"
    frames: list[pd.DataFrame] = []
    gtex = raw / "gtex" / "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
    if gtex.exists():
        gt = pd.read_csv(gtex, sep="\t", skiprows=2, compression="gzip", low_memory=False)
        gt["gene_symbol"] = gt["Description"].map(_clean_gene)
        tissue_cols = [c for c in gt.columns if c not in {"Name", "Description", "gene_symbol"}]
        vital_cols = [c for c in tissue_cols if any(term.lower() in c.lower() for term in VITAL_GTEX_TERMS)]
        tmp = pd.DataFrame(index=gt["gene_symbol"])
        tmp["gtex_max_median_tpm"] = gt[tissue_cols].max(axis=1).to_numpy()
        tmp["gtex_vital_max_median_tpm"] = (gt[vital_cols].max(axis=1) if vital_cols else gt[tissue_cols].max(axis=1)).to_numpy()
        tmp["normal_expression_risk_score"] = np.log2(tmp["gtex_vital_max_median_tpm"].fillna(0) + 1).clip(0, 8) / 8
        frames.append(tmp.groupby(tmp.index).max())
    gn = raw / "gnomad" / "gnomad.v4.1.1.constraint_metrics.tsv.bgz"
    if gn.exists():
        cols = ["gene", "canonical", "lof.oe_ci.upper", "lof.pLI"]
        gd = pd.read_csv(gn, sep="\t", compression="gzip", usecols=lambda c: c in cols, low_memory=False)
        if "canonical" in gd.columns:
            gd = gd[gd["canonical"].astype(str).str.lower().isin(["true", "1"])]
        gd["gene_symbol"] = gd["gene"].map(_clean_gene)
        tmp = gd.drop_duplicates("gene_symbol").set_index("gene_symbol")
        loeuf = pd.to_numeric(tmp.get("lof.oe_ci.upper"), errors="coerce")
        pli = pd.to_numeric(tmp.get("lof.pLI"), errors="coerce")
        frames.append(pd.DataFrame({
            "gnomad_loeuf": loeuf,
            "gnomad_pli": pli,
            "lof_constraint_risk_score": pd.concat([_safe_minmax(loeuf, invert=True), pli.fillna(0).clip(0, 1)], axis=1).max(axis=1),
        }))
    if common_essential is not None:
        ce = common_essential.astype(float).rename("common_essential_risk_score").to_frame()
        frames.append(ce)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, axis=1).fillna(0)
    risk_cols = [c for c in ["normal_expression_risk_score", "lof_constraint_risk_score", "common_essential_risk_score"] if c in out.columns]
    out["safety_risk_score"] = out[risk_cols].mean(axis=1).clip(0, 1) if risk_cols else 0.0
    out["therapeutic_window_score"] = (1.0 - out["safety_risk_score"]).clip(0, 1)
    return out.groupby(out.index).max()


def extract_literature_benchmark(root: Path, gene_universe: set[str]) -> pd.DataFrame:
    """Extract usable benchmark target rows from templates and supplements."""
    lit = root / "data" / "raw" / "literature_benchmark"
    rows: list[dict[str, Any]] = []
    for tsv in lit.glob("*_targets.tsv"):
        source = tsv.name.replace("_targets.tsv", "")
        try:
            df = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        if "gene_symbol" in df.columns and len(df):
            for _, row in df.iterrows():
                gene = _clean_gene(row.get("gene_symbol"))
                if gene:
                    rows.append({"gene_symbol": gene, "source": source, "target_category": row.get("target_category", ""), "reported_score": row.get("reported_rank_or_score", row.get("reported_score", "")), "extraction_method": "manual_tsv"})
    for source_dir in ["nwosu2025", "tcga_depmap_shi2024"]:
        base = lit / source_dir
        if not base.exists():
            continue
        found: dict[str, int] = {}
        files = list(base.rglob("*.xlsx")) + list(base.rglob("*.csv"))
        for path in files:
            try:
                if path.suffix.lower() == ".xlsx":
                    with pd.ExcelFile(path, engine="openpyxl") as xls:
                        sheets = xls.sheet_names[:20]
                        for sheet in sheets:
                            df = pd.read_excel(xls, sheet_name=sheet, header=None, nrows=5000)
                            for val in pd.unique(df.astype(str).values.ravel()):
                                gene = _clean_gene(val)
                                if gene in gene_universe and GENE_RE.match(gene):
                                    found[gene] = found.get(gene, 0) + 1
                else:
                    df = pd.read_csv(path, low_memory=False, nrows=50000)
                    for val in pd.unique(df.astype(str).values.ravel()):
                        gene = _clean_gene(val)
                        if gene in gene_universe and GENE_RE.match(gene):
                            found[gene] = found.get(gene, 0) + 1
            except Exception as exc:
                if path.suffix.lower() == ".xlsx":
                    raise RuntimeError(
                        f"Could not read literature workbook {path}. "
                        "Install the dependencies in requirements.txt and verify the workbook."
                    ) from exc
                continue
        for gene, count in found.items():
            rows.append({"gene_symbol": gene, "source": source_dir, "target_category": "extracted_supplement_gene", "reported_score": count, "extraction_method": "supplement_symbol_scan"})
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["gene_symbol", "source", "target_category", "reported_score", "extraction_method"])
    return out.drop_duplicates(["gene_symbol", "source", "target_category"])


def write_geo_qc(root: Path) -> pd.DataFrame:
    """Write external GEO cohort QC without forcing probe IDs into gene scores."""
    rows: list[dict[str, Any]] = []
    raw = root / "data" / "raw" / "external_geo"
    out_dir = root / "data" / "processed" / "external_geo"
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in sorted(raw.glob("GSE*_series_matrix.txt.gz")):
        metadata: dict[str, list[str]] = {}
        n_rows = 0
        n_samples = 0
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            in_table = False
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith(("!Sample_title", "!Sample_source_name_ch1", "!Sample_characteristics_ch1", "!Sample_geo_accession", "!Series_platform_id")):
                    parts = [p.strip().strip('"') for p in line.split("\t")]
                    metadata.setdefault(parts[0], []).append("|".join(parts[1:]))
                elif line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    header = next(handle).rstrip("\n").split("\t")
                    n_samples = max(0, len(header) - 1)
                elif in_table:
                    if line.startswith("!series_matrix_table_end"):
                        break
                    n_rows += 1
        text = " ".join(" ".join(v) for v in metadata.values()).lower()
        tumor_mentions = len(re.findall(r"\btumou?r\b|pancreatic tumor|cancer", text))
        normal_mentions = len(re.findall(r"\bnormal\b|non-tumou?r|adjacent", text))
        platform = metadata.get("!Series_platform_id", [""])[0].replace('"', "")
        gene_level_likely = "GPL20769" in platform
        rows.append({
            "series": path.name.split("_")[0],
            "platform": platform,
            "n_samples": n_samples,
            "n_matrix_rows": n_rows,
            "tumor_metadata_mentions": tumor_mentions,
            "normal_metadata_mentions": normal_mentions,
            "gene_level_likely": gene_level_likely,
            "usable_for_gene_level_scoring": bool(gene_level_likely and normal_mentions > 0 and tumor_mentions > 0),
            "note": "Probe annotation required before gene-level scoring" if not gene_level_likely else "Gene-symbol-like rows; sample design still requires manual review",
        })
    qc = pd.DataFrame(rows)
    qc.to_csv(out_dir / "geo_cohort_qc.tsv", sep="\t", index=False)
    return qc


def _read_geo_platform_mapping(path: Path) -> pd.DataFrame:
    """Read GEO platform annotation into probe-to-gene mapping."""
    if not path.exists():
        return pd.DataFrame(columns=["probe_id", "gene_symbol"])
    rows: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                rows.append(line)
    if not rows:
        return pd.DataFrame(columns=["probe_id", "gene_symbol"])
    annot = pd.read_csv(io.StringIO("".join(rows)), sep="\t", low_memory=False)
    id_col = "ID" if "ID" in annot.columns else annot.columns[0]
    gene_col = next((c for c in ["Gene symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol"] if c in annot.columns), None)
    if gene_col is None:
        return pd.DataFrame(columns=["probe_id", "gene_symbol"])
    out = annot[[id_col, gene_col]].rename(columns={id_col: "probe_id", gene_col: "gene_symbol"})
    out["probe_id"] = out["probe_id"].astype(str)
    out["gene_symbol"] = (
        out["gene_symbol"].astype(str)
        .str.replace("///", ";", regex=False)
        .str.replace("//", ";", regex=False)
        .str.split(";").str[0]
        .map(_clean_gene)
    )
    out = out[out["probe_id"].astype(bool) & out["gene_symbol"].astype(bool)]
    return out.drop_duplicates("probe_id")


def _load_geo_platform_maps(root: Path) -> dict[str, pd.DataFrame]:
    """Load available platform maps from downloaded GEO annotation files."""
    plat = root / "data" / "raw" / "external_geo" / "platforms"
    paths = {
        "GPL570": plat / "GPL570.annot.gz",
        "GPL6244": plat / "GPL6244.annot.gz",
        "GPL13667": plat / "GPL13667_family.soft.gz",
    }
    return {key: _read_geo_platform_mapping(path) for key, path in paths.items() if path.exists()}


def _parse_geo_series_matrix(path: Path) -> tuple[str, pd.DataFrame, pd.Series]:
    """Parse a GEO series matrix and infer tumor/normal sample labels."""
    meta: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line.replace('"', ""))
            elif line.startswith(("!Sample_title", "!Sample_source_name_ch1", "!Sample_characteristics_ch1", "!Sample_geo_accession", "!Series_platform_id")):
                parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
                meta.setdefault(parts[0], []).append(parts[1:])
    platform = meta.get("!Series_platform_id", [[""]])[0][0]
    if not table_lines:
        return platform, pd.DataFrame(), pd.Series(dtype=str)
    matrix = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", low_memory=False)
    matrix = matrix.rename(columns={matrix.columns[0]: "probe_id"})
    matrix["probe_id"] = matrix["probe_id"].astype(str)
    accessions = meta.get("!Sample_geo_accession", [[]])[0]
    sample_cols = [c for c in matrix.columns if c != "probe_id"]
    labels = {}
    title = meta.get("!Sample_title", [[]])[0]
    source = meta.get("!Sample_source_name_ch1", [[]])[0]
    chars = meta.get("!Sample_characteristics_ch1", [])
    for i, col in enumerate(sample_cols):
        bits = []
        for values in [title, source, *chars]:
            if i < len(values):
                bits.append(values[i])
        text = " ".join(bits).lower()
        if re.search(r"non[- ]?tumou?r|normal|adjacent|sample:\s*n\b|tissue:\s*n\b|nontumor", text):
            labels[col] = "normal"
        elif re.search(r"tumou?r|cancer|sample:\s*t\b|tissue:\s*t\b|pdac", text):
            labels[col] = "tumor"
        else:
            labels[col] = "other"
    return platform, matrix, pd.Series(labels)


def build_geo_validation_features(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build gene-level external validation from GEO tumor-normal cohorts."""
    maps = _load_geo_platform_maps(root)
    raw = root / "data" / "raw" / "external_geo"
    processed = root / "data" / "processed" / "external_geo"
    processed.mkdir(parents=True, exist_ok=True)
    cohort_tables: list[pd.DataFrame] = []
    cohort_summary: list[dict[str, Any]] = []
    for path in sorted(raw.glob("GSE*_series_matrix.txt.gz")):
        series = path.name.split("_")[0]
        platform, matrix, labels = _parse_geo_series_matrix(path)
        if matrix.empty:
            continue
        tumor_cols = [c for c in labels.index if labels[c] == "tumor" and c in matrix.columns]
        normal_cols = [c for c in labels.index if labels[c] == "normal" and c in matrix.columns]
        mapping = maps.get(platform)
        if platform == "GPL20769" and mapping is None:
            mapping = pd.DataFrame({"probe_id": matrix["probe_id"].astype(str), "gene_symbol": matrix["probe_id"].map(_clean_gene)})
        usable = mapping is not None and len(mapping) > 0 and len(tumor_cols) >= 3 and len(normal_cols) >= 3
        cohort_summary.append({
            "series": series,
            "platform": platform,
            "n_tumor": len(tumor_cols),
            "n_normal": len(normal_cols),
            "mapping_rows": 0 if mapping is None else len(mapping),
            "used_for_external_score": bool(usable),
        })
        if not usable:
            continue
        expr = matrix.merge(mapping, on="probe_id", how="inner")
        value_cols = tumor_cols + normal_cols
        expr[value_cols] = expr[value_cols].apply(pd.to_numeric, errors="coerce")
        if np.nanpercentile(expr[value_cols].to_numpy(float), 95) > 100:
            expr[value_cols] = np.log2(expr[value_cols].clip(lower=0) + 1)
        gene_expr = expr.groupby("gene_symbol")[value_cols].median()
        t = gene_expr[tumor_cols]
        n = gene_expr[normal_cols]
        logfc = t.mean(axis=1) - n.mean(axis=1)
        pvals = []
        for gene in gene_expr.index:
            tv = t.loc[gene].dropna().to_numpy(float)
            nv = n.loc[gene].dropna().to_numpy(float)
            if len(tv) >= 3 and len(nv) >= 3:
                pvals.append(float(stats.ttest_ind(tv, nv, equal_var=False, nan_policy="omit").pvalue))
            else:
                pvals.append(np.nan)
        out = pd.DataFrame({
            "gene_symbol": gene_expr.index,
            "series": series,
            "platform": platform,
            "geo_log2fc": logfc.to_numpy(),
            "geo_pvalue": pvals,
            "geo_direction": np.sign(logfc.to_numpy()),
            "geo_significant": ((logfc.abs() >= 0.5) & (pd.Series(pvals, index=logfc.index) < 0.05)).astype(int).to_numpy(),
            "n_tumor": len(tumor_cols),
            "n_normal": len(normal_cols),
        })
        out.to_csv(processed / f"{series}_gene_level_de.tsv", sep="\t", index=False)
        cohort_tables.append(out)
    summary = pd.DataFrame(cohort_summary)
    summary.to_csv(processed / "geo_external_validation_cohort_summary.tsv", sep="\t", index=False)
    if not cohort_tables:
        return pd.DataFrame(), summary
    all_de = pd.concat(cohort_tables, ignore_index=True)
    all_de.to_csv(processed / "geo_external_validation_gene_level.tsv", sep="\t", index=False)
    agg = all_de.groupby("gene_symbol").agg(
        geo_validation_cohorts=("series", "nunique"),
        geo_up_count=("geo_direction", lambda x: int((x > 0).sum())),
        geo_down_count=("geo_direction", lambda x: int((x < 0).sum())),
        geo_significant_count=("geo_significant", "sum"),
        geo_median_log2fc=("geo_log2fc", "median"),
    )
    agg["geo_direction_consistency"] = agg[["geo_up_count", "geo_down_count"]].max(axis=1) / agg["geo_validation_cohorts"].clip(lower=1)
    agg["geo_external_reproducibility_score"] = (
        (agg["geo_significant_count"].clip(0, 3) / 3) * agg["geo_direction_consistency"].fillna(0)
    ).clip(0, 1)
    return agg, summary


def load_known_targets(root: Path) -> dict[str, list[str]]:
    """Load known target categories from config, falling back to defaults."""
    path = root / "config" / "known_targets.yaml"
    if not path.exists():
        return KNOWN_TARGET_DEFAULTS
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    categories = payload.get("known_target_categories", {})
    return {key: [_clean_gene(g) for g in value.get("genes", [])] for key, value in categories.items()}


def assign_triage(df: pd.DataFrame) -> pd.Series:
    """Assign target triage classes from the two scores and limiting factors."""
    classes = []
    for _, row in df.iterrows():
        cas = row.get("clinical_actionability_score", 0)
        dps = row.get("discovery_potential_score", 0)
        cas_pct = row.get("clinical_actionability_percentile", 0)
        dps_pct = row.get("discovery_potential_percentile", 0)
        dep = row.get("dependency_score", 0)
        drug = row.get("druggability_score", 0)
        safety = row.get("safety_risk_score", 0)
        disease = row.get("disease_relevance_score", 0)
        if safety >= 0.75:
            cls = "Class IV: high-risk target"
        elif (cas >= 0.70 or cas_pct >= 0.985) and drug >= 0.25 and dep >= 0.15 and disease >= 0.10 and safety <= 0.65:
            cls = "Class I: clinically actionable candidate"
        elif (dps >= 0.70 or dps_pct >= 0.985) and dep >= 0.30 and disease >= 0.08 and safety <= 0.70:
            cls = "Class II: discovery-priority candidate"
        elif dep >= 0.45 and drug < 0.15:
            cls = "Class III: strong dependency, weak tractability"
        elif disease >= 0.35 and dep < 0.15:
            cls = "Class V: biomarker-like / low dependency"
        else:
            cls = "Class VI: insufficient or mixed evidence"
        classes.append(cls)
    return pd.Series(classes, index=df.index)


def build_scores(matrix: pd.DataFrame) -> pd.DataFrame:
    """Calculate disease relevance, dual scores and final triage classes."""
    df = matrix.copy()
    for col in [
        "patient_alteration_score", "proteomic_support_score", "external_reproducibility_score",
        "dependency_score", "druggability_score", "therapeutic_window_score", "target_tractability_score",
        "safety_risk_score",
    ]:
        if col not in df.columns:
            df[col] = 0.0
    patient = df["patient_alteration_score"].fillna(0).clip(0, 1)
    protein = df["proteomic_support_score"].fillna(0).clip(0, 1)
    external = df["external_reproducibility_score"].fillna(0).clip(0, 1)
    dependency = df["dependency_score"].fillna(0).clip(0, 1)
    druggability = df["druggability_score"].fillna(0).clip(0, 1)
    therapeutic_window = df["therapeutic_window_score"].fillna(0).clip(0, 1)
    tractability = df["target_tractability_score"].fillna(0).clip(0, 1)
    safety = df["safety_risk_score"].fillna(0).clip(0, 1)

    disease_components = pd.DataFrame({
        "patient_alteration_score": patient,
        "proteomic_support_score": protein,
        "external_reproducibility_score": external,
    }, index=df.index)
    df["disease_relevance_score"] = _mean_frame(disease_components, list(disease_components.columns))
    safety_factor = (1.0 - safety).clip(0, 1)
    df["safety_factor"] = safety_factor
    df["safety_modifier"] = safety_factor
    clinical_components = pd.DataFrame({
        "disease_relevance_score": df["disease_relevance_score"],
        "dependency_score": dependency,
        "druggability_score": druggability,
    }, index=df.index)
    discovery_components = pd.DataFrame({
        "disease_relevance_score": df["disease_relevance_score"],
        "dependency_score": dependency,
        "therapeutic_window_score": therapeutic_window,
    }, index=df.index)
    df["clinical_actionability_core"] = _geomean_frame(clinical_components, list(clinical_components.columns))
    df["clinical_actionability_score"] = df["clinical_actionability_core"].mul(safety_factor).clip(0, 1)
    df["discovery_potential_core"] = _geomean_frame(discovery_components, list(discovery_components.columns))
    df["discovery_potential_score"] = df["discovery_potential_core"].mul(np.sqrt(tractability)).clip(0, 1)
    df["clinical_actionability_percentile"] = df["clinical_actionability_score"].rank(pct=True, ascending=True)
    df["discovery_potential_percentile"] = df["discovery_potential_score"].rank(pct=True, ascending=True)
    df["clinical_actionability_rank"] = df["clinical_actionability_score"].rank(ascending=False, method="min").astype(int)
    df["discovery_potential_rank"] = df["discovery_potential_score"].rank(ascending=False, method="min").astype(int)
    df["triage_class"] = assign_triage(df)
    return df.sort_values("clinical_actionability_rank")


def _save_indexed_table(df: pd.DataFrame, path: Path) -> None:
    """Save a gene-indexed table with an explicit gene_symbol column."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.reset_index().rename(columns={"index": "gene_symbol"}).to_csv(path, sep="\t", index=False)
