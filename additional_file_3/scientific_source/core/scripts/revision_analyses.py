"""Reviewer-requested sensitivity and audit analyses, 2026-10-08.

Frozen target-level TSV inputs; no patient records. Random seeds are explicit.
Environment: results/environment.json. Primary model: explicit elastic-net.
Descriptive perturbation ranges are not biological confidence intervals.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_config
from scipy.stats import spearmanr, binom
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from run_revision import ROOT, prepare, save
from framework import pipeline as fp

np.random.seed(42)

FAMILIES = {
    "drug_knowledge": {"chemical_modulation_score", "clinical_precedence_score", "druggability_score"},
    "dependency": {"dependency_score", "safety_risk_score", "safety_factor"},
    "patient_omics": {"patient_alteration_score", "tcga_gtex_expression_score", "tcga_mutation_score", "tcga_cnv_score", "tcga_survival_score"},
    "protein_omics": {"proteomic_support_score", "cptac_mrna_score", "cptac_protein_score", "phosphoproteomic_support_score"},
    "external_knowledge": {"opentargets_overall_score", "external_reproducibility_score", "geo_external_reproducibility_score"},
}


def axis_data(cancer):
    matrix, labels, benchmarks, roles, features, cfg, out = prepare(cancer, write=False)
    X = features.set_index("gene_symbol")
    lab = labels.set_index("gene_symbol")
    scores = pd.read_csv(out / "framework_priority_scores.tsv", sep="\t").set_index("gene_symbol")
    reserved_bench = set().union(*benchmarks.values())
    for axis in fp.AXES:
        pos = set(lab.index[lab[axis.label_column].eq(1)])
        yield axis, X[fp.axis_feature_columns(features, cfg, axis.name)], pos, pos | reserved_bench, benchmarks[axis.name], scores[axis.score_column], matrix, cfg, out, roles


def fit_one(X, pos, reserved, seeds, cfg, strata, matched=True, ratio=8):
    import fast_model
    fast_model.install(fp)
    cfg = dict(cfg, compute_repeat_rank_frequency=False)
    with threadpool_limits(limits=1):
        score, coef, audit = fp.fit_sampled_background_model(X, pos, reserved, seeds, cfg, strata, matched, ratio)
    return score.set_index("gene_symbol")["score_mean"].reindex(X.index), len(audit), int(audit.used_in_ensemble.sum())


def null_fit(X, pos, reserved_bench, benchmark, cfg, strata, iteration):
    seeds = [int(s) + iteration * 1000 for s in cfg["random_seeds"]]
    score, attempts, accepted = fit_one(X, pos, pos | reserved_bench, seeds, cfg, strata)
    top = set(score.nlargest(50).index)
    return {"permutation": iteration, "recovery": len(top & benchmark) / len(benchmark),
            "attempted_fits": attempts, "accepted_fits": accepted}


def permutation(cancer, jobs):
    rows = []
    for axis, X, pos, reserved, benchmark, primary, matrix, cfg, out, roles in axis_data(cancer):
        reserved_bench = set(roles.loc[roles.record_role.eq("V"), "gene_symbol"])
        # Match the observed axis-specific label embargo: a development label
        # may occur in the other axis's validation list. Exclude only this
        # axis's held-out targets from randomized positive labels. Background
        # sampling continues to exclude the union of both validation lists.
        eligible = np.array(sorted(set(X.index) - set(benchmark)))
        strata = fp.build_background_strata(matrix)
        stratum_pools = {s: members.index.to_numpy() for s, members in strata.reindex(eligible).groupby(strata.reindex(eligible))}
        positive_stratum_counts = strata.reindex(sorted(pos)).value_counts()
        observed = len(set(primary.nlargest(50).index) & benchmark) / len(benchmark)
        for null, B in [("uniform_gene_labels", 500), ("coverage_conditioned_labels", 500)]:
            rng = np.random.default_rng(42 + (0 if axis.name == "translational" else 10000))
            sets = []
            label_audit = []
            for i in range(B):
                fallback_count = 0
                if null == "uniform_gene_labels":
                    pp = set(rng.choice(eligible, size=len(pos), replace=False))
                else:
                    # Same annotation/measurement/expression strata as observed
                    # positives, with the same documented fallback as sampling.
                    pp = set()
                    for stratum, count in positive_stratum_counts.items():
                        pool = stratum_pools.get(stratum, np.array([], dtype=object))
                        draw = min(int(count), len(pool))
                        if draw:
                            pp.update(rng.choice(pool, draw, replace=False))
                    if len(pp) < len(pos):
                        fallback_count = len(pos)-len(pp)
                        pool = [g for g in eligible if g not in pp]
                        pp.update(rng.choice(pool, len(pos)-len(pp), replace=False))
                    if len(pp) != len(pos):
                        raise RuntimeError("Conditional null changed positive count")
                sets.append(pp)
                assert len(pp)==len(pos) and not (pp & set(benchmark))
                label_audit.append({"permutation":i,"positive_genes":";".join(sorted(pp)),
                    "positive_genes_json":json.dumps(sorted(pp)),
                    "n_positive":len(pp),"n_eligible_positive_candidates":len(eligible),
                    "conditional_label_fallback_count":fallback_count,
                    "other_axis_validation_genes_as_positive":len(pp & (reserved_bench-set(benchmark))),
                    "positive_embargo":"current-axis validation list", "background_embargo":"both-axis validation union plus sampled positives"})
            save(pd.DataFrame(label_audit), out / f"revision_permutation_labels_{axis.name}_{null}.tsv")
            print(f"{cancer} {axis.name}: {null}, {B} permutations", flush=True)
            with parallel_config(backend="loky", inner_max_num_threads=1):
                detail = Parallel(n_jobs=jobs)(delayed(null_fit)(X, pp, reserved_bench, benchmark, cfg, strata, i) for i, pp in enumerate(sets))
            d = pd.DataFrame(detail)
            d["axis"] = axis.name
            d["null"] = null
            save(d, out / f"revision_permutation_{axis.name}_{null}.tsv")
            p = (int(d.recovery.ge(observed).sum()) + 1) / (B + 1)
            rows.append({"axis": axis.name, "null": null, "metric": "top50_recovery",
                         "evaluation_set": "heldout_benchmark", "n_development_positives": len(pos),
                         "n_benchmark_targets": len(benchmark), "n_model_repeats_per_permutation": 30,
                         "n_permutations": B, "n_converged_fits": int(d.accepted_fits.sum()),
                         "n_fit_attempts": int(d.attempted_fits.sum()),
                         "n_numerical_retry_attempts": int((d.attempted_fits-d.accepted_fits).sum()),
                         "fit_convergence_rate": d.accepted_fits.sum()/d.attempted_fits.sum(),
                         "observed_value": observed, "permutation_mean": d.recovery.mean(),
                         "permutation_sd": d.recovery.std(ddof=0), "empirical_p_value": p})
    result = pd.DataFrame(rows)
    save(result, out / "revision_permutation_summary.tsv")
    save(result[result.null.eq("uniform_gene_labels")].drop(columns=["null"]), out / "framework_label_permutation_test.tsv")


def perturb_fit(X, pos, reserved, benchmark, primary, cfg, strata, iteration, kind):
    rng = np.random.default_rng(42000 + iteration)
    if kind == "label_thinning":
        n_remove = min(len(pos)-2, max(1, int(np.ceil(.20*len(pos)))))
        removed = set(rng.choice(sorted(pos), n_remove, replace=False))
        training = pos-removed
    else:
        removed, training = set(), pos
    seeds = list(range(100000+iteration*100, 100030+iteration*100))
    score, attempts, accepted = fit_one(X, training, reserved, seeds, cfg, strata)
    ranks = score.rank(ascending=False, method="first").to_numpy()
    top = set(score.nlargest(50).index)
    primary_top = set(primary.nlargest(50).index)
    result = {"replicate": iteration, "analysis": kind, "n_training_positives": len(training),
              "removed_genes": ";".join(sorted(removed)),
              "spearman_with_primary": spearmanr(primary.reindex(X.index), score).statistic,
              "top50_overlap_fraction": len(top & primary_top)/50,
              "top50_recovery": len(top & benchmark)/len(benchmark),
              "attempted_fits": attempts, "accepted_fits": accepted}
    return result, ranks


def perturbation(cancer, jobs):
    summaries = []
    for axis, X, pos, reserved, benchmark, primary, matrix, cfg, out, roles in axis_data(cancer):
        strata = fp.build_background_strata(matrix)
        for kind, R in [("background_ensemble", 20), ("label_thinning", 50)]:
            print(f"{cancer} {axis.name}: {kind}, {R} ensembles", flush=True)
            with parallel_config(backend="loky", inner_max_num_threads=1):
                result = Parallel(n_jobs=jobs)(delayed(perturb_fit)(X,pos,reserved,benchmark,primary,cfg,strata,i,kind) for i in range(R))
            detail = pd.DataFrame([r[0] for r in result]); detail["axis"] = axis.name
            save(detail, out / f"revision_{kind}_{axis.name}_replicates.tsv")
            ranks = np.stack([r[1] for r in result], axis=1)
            frequency = (ranks <= 50).mean(axis=1)
            candidates = pd.DataFrame({"gene_symbol": X.index, "axis": axis.name,
                                       "analysis": kind, "n_replicates": R,
                                       "top50_frequency": frequency,
                                       "median_rank": np.median(ranks,axis=1),
                                       "rank_p025": np.quantile(ranks,.025,axis=1),
                                       "rank_p975": np.quantile(ranks,.975,axis=1)})
            save(candidates, out / f"revision_{kind}_{axis.name}_candidate_stability.tsv")
            for metric in ["top50_overlap_fraction","top50_recovery","spearman_with_primary"]:
                values=detail[metric]
                summaries.append({"cancer":cancer,"axis":axis.name,"analysis":kind,"metric":metric,
                                  "n_replicates":R,"estimate":values.median(),
                                  "range_p025":values.quantile(.025),"range_p975":values.quantile(.975),
                                  "interval_type":"empirical perturbation range; not a biological CI"})
    save(pd.DataFrame(summaries), out / "revision_stability_summary.tsv")


def descendant_exclusions(cancer, family):
    remove=set(FAMILIES[family])
    if cancer=="LUAD" and family=="dependency":
        remove.add("external_reproducibility_score")
    if cancer=="LUAD" and family=="patient_omics":
        # TCGA-derived RPPA fallback is part of the historical LUAD inputs.
        remove.update(["cptac_mrna_score","cptac_protein_score","proteomic_support_score","external_reproducibility_score"])
    return remove | {f"{c}_missing" for c in remove}


def strata_without_family(matrix, family):
    remove=[]
    for c in matrix:
        if family=="drug_knowledge" and any(k in c for k in ["drugcentral","dgidb","ttd_","iuphar","clinical_precedence"]): remove.append(c)
        if family=="dependency" and c.startswith("depmap_"): remove.append(c)
        if family=="patient_omics" and c.startswith("tcga_"): remove.append(c)
        if family=="external_knowledge" and c.startswith("opentargets_"): remove.append(c)
    return fp.build_background_strata(matrix.drop(columns=remove))


def fixed_ablation(cancer):
    rows=[]
    for axis,X,pos,reserved,benchmark,primary,matrix,cfg,out,roles in axis_data(cancer):
        for family in FAMILIES:
            excluded=descendant_exclusions(cancer,family)
            columns=[c for c in X if c not in excluded]
            if not columns: continue
            score,attempts,accepted=fit_one(X[columns],pos,reserved,cfg["random_seeds"],cfg,strata_without_family(matrix,family))
            top=set(score.nlargest(50).index)
            rows.append({"cancer":cancer,"axis":axis.name,"removed_family":family,
                         "n_training_positives":len(pos),"n_benchmark_targets":len(benchmark),
                         "excluded_model_columns":";".join(sorted(set(X.columns)&excluded)),
                         "n_features_remaining":len(columns),"top50_hits":len(top&benchmark),
                         "top50_recovery":len(top&benchmark)/len(benchmark),
                         "spearman_with_primary":spearmanr(primary,score.reindex(primary.index)).statistic,
                         "top50_overlap_fraction":len(top&set(primary.nlargest(50).index))/50,
                         "attempted_fits":attempts,"accepted_fits":accepted})
            save(score.rename("score").reset_index(),out/f"revision_ablation_{axis.name}_{family}_scores.tsv")
    save(pd.DataFrame(rows),out/"revision_fixed_benchmark_source_ablation.tsv")


def common_review(matrix, cfg, axis):
    m=matrix.set_index("gene_symbol")
    dep=m.dependency_score.fillna(0)
    disease=m.disease_relevance_score.fillna(0)
    safety=m.safety_factor.fillna(0).ge(cfg["review"]["safety_factor_min"]) & m.safety_risk_score.fillna(1).lt(cfg["review"]["safety_risk_high_threshold"])
    dep_gate=dep.ge(.5*dep.quantile(.9))
    if axis=="translational":
        override=m.clinical_precedence_score.fillna(0).ge(.6) | m.opentargets_overall_score.fillna(0).ge(.5)
        return safety & m.druggability_score.fillna(0).ge(.10) & (dep_gate | override)
    return safety & dep_gate & disease.ge(.5*disease.quantile(.9)) & m.druggability_score.fillna(0).ge(.15)


def fair_comparisons(cancer):
    rows=[]
    for axis,X,pos,reserved,benchmark,primary,matrix,cfg,out,roles in axis_data(cancer):
        strata=fp.build_background_strata(matrix)
        en,_,_=fit_one(X,pos,reserved,cfg["random_seeds"],cfg,strata)
        # Match information, solver, scaling, tolerance, seeds and training data.
        original=fp.make_elastic_net_model
        def same_solver_l2(C,l1_ratio,max_iter,solver="saga"):
            return Pipeline([("scale",StandardScaler()),("clf",LogisticRegression(penalty="elasticnet",l1_ratio=0.0,solver=solver,C=C,max_iter=max_iter,tol=1e-2,class_weight="balanced",random_state=0))])
        fp.make_elastic_net_model=same_solver_l2
        try:
            l2,_,_=fit_one(X,pos,reserved,cfg["random_seeds"],cfg,strata)
        finally:
            fp.make_elastic_net_model=original
        methods={"Elastic-net (same information)":en,
                 "L2 (same information)":l2,
                 "Mean percentile rank (same information)":pd.concat([1-X[c].rank(pct=True,method="average") if c=="safety_risk_score" else X[c].rank(pct=True,method="average") for c in X if not c.endswith("_missing")],axis=1).mean(axis=1)}
        if axis.name=="translational":
            cols=["clinical_precedence_score","chemical_modulation_score","druggability_score","opentargets_overall_score"]
            methods["Rank aggregation (native information)"]=fp._rank_aggregation(matrix.set_index("gene_symbol"),cols)
        else:
            methods["Dependency-only (native information)"]=matrix.set_index("gene_symbol").dependency_score.fillna(0)
        gate=common_review(matrix,cfg,axis.name).reindex(X.index)
        rng=np.random.default_rng(42)
        genes=np.array(sorted(benchmark))
        for name,s in methods.items():
            top=set(s.nlargest(50).index)
            top_reviewed=set(s[gate].nlargest(50).index)
            hits=np.array([g in top for g in genes],dtype=float)
            boot=binom.ppf([.025,.975],len(hits),hits.mean())/len(hits)
            rows.append({"cancer":cancer,"axis":axis.name,"method":name,
                         "n_benchmark":len(benchmark),"top50_hits":int(hits.sum()),
                         "top50_recovery":hits.mean(),"benchmark_bootstrap_ci_low":boot[0],
                         "benchmark_bootstrap_ci_high":boot[1],
                         "top50_common_review_retained":len(top&set(gate.index[gate])),
                         "top50_after_common_review_hits":len(top_reviewed&benchmark),
                         "eligible_universe_count":int(gate.sum()),
                         "review_metrics_are_descriptive":True})
    save(pd.DataFrame(rows),out/"revision_information_matched_comparators.tsv")


def overlap_and_class_vi(cancer):
    matrix,labels,benchmarks,roles,features,cfg,out=prepare(cancer,write=False)
    idx=matrix.set_index("gene_symbol")
    feature_rows=[]
    for axis in fp.AXES:
        bm=idx.loc[sorted(benchmarks[axis.name])]
        selected=set(fp.axis_feature_columns(features,cfg,axis.name))
        for family,columns in FAMILIES.items():
            columns=sorted(columns&set(idx.columns))
            positive=bm[columns].fillna(0).gt(0).any(axis=1)
            feature_rows.append({"cancer":cancer,"axis":axis.name,"feature_family":family,
                                 "benchmark_count":len(bm),"benchmark_with_positive_source_signal":int(positive.sum()),
                                 "fraction":positive.mean(),"retained_model_columns":";".join(sorted(set(columns)&selected)),
                                 "interpretation":"target-level annotation overlap proxy; not proof of shared facts"})
    save(pd.DataFrame(feature_rows),out/"revision_benchmark_feature_overlap.tsv")
    classes=pd.read_csv(out/"framework_review_classes.tsv",sep="\t").set_index("gene_symbol")
    blocks={"patient":["patient_alteration_score"],"protein":["proteomic_support_score"],
            "functional":["dependency_score"],"drug":["chemical_modulation_score","clinical_precedence_score"],
            "external":["external_reproducibility_score","opentargets_overall_score","geo_external_reproducibility_score"]}
    count=sum(idx[cols].fillna(0).gt(0).any(axis=1).astype(int) for cols in blocks.values())
    high=classes[fp.CLINICAL_SCORE].ge(classes[fp.CLINICAL_SCORE].quantile(.98)) | classes[fp.DISCOVERY_SCORE].ge(classes[fp.DISCOVERY_SCORE].quantile(.98))
    vi=classes.review_class.str.startswith("Class VI:")
    reason=pd.Series("multiple positive blocks below both score gates",index=idx.index)
    reason.loc[count.eq(0)]="no positive support block"
    reason.loc[count.eq(1)]="one positive support block"
    reason.loc[count.ge(2)&high]="multiple positive blocks with high score but review constraints"
    partition=reason[vi].value_counts().rename_axis("reason").reset_index(name="count")
    assert int(partition['count'].sum())==int(vi.sum())
    partition["fraction_of_class_vi"]=partition['count']/vi.sum()
    partition["cancer"]=cancer
    partition["measurement_note"]="positive-support proxy; unavailable and observed zero may coincide"
    save(partition,out/"revision_class_vi_partition.tsv")
    save(pd.DataFrame({"gene_symbol":idx.index,"positive_support_block_count":count,
                       "class_vi":vi,"class_vi_reason":reason.where(vi,"")}),out/"revision_class_vi_gene_audit.tsv")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",choices=["permutation","perturbation","extended"],required=True)
    parser.add_argument("--cancer",choices=["LUAD","PDAC","both"],default="both")
    parser.add_argument("--jobs",type=int,default=8)
    args=parser.parse_args()
    for cancer in (["LUAD","PDAC"] if args.cancer=="both" else [args.cancer]):
        if args.stage=="permutation": permutation(cancer,args.jobs)
        elif args.stage=="perturbation": perturbation(cancer,args.jobs)
        else:
            fixed_ablation(cancer)
            fair_comparisons(cancer)
            overlap_and_class_vi(cancer)
        print(f"{cancer}: {args.stage} completed",flush=True)
