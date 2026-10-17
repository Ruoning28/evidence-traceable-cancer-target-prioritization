"""Numerically equivalent ensemble engine with cached sampling pools.

Scaling is fitted only on each sampled training set. Full-universe scoring uses
an algebraically equivalent affine map. Convergence retries retain the same data
and seed, rather than selecting new backgrounds conditional on convergence.
"""
from __future__ import annotations
import warnings
import json,os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits


def install(fp):
    fp.fit_sampled_background_model = fit_sampled_background_model


def fit_sampled_background_model(X, positive_genes, all_positive_genes, seeds,
                                 cfg, strata, matched_background, background_ratio,
                                 penalty="elasticnet"):
    from framework import pipeline as fp
    genes = X.index.astype(str)
    positive = [g for g in genes if g in positive_genes]
    unlabeled = [g for g in genes if g not in all_positive_genes]
    if len(positive)<2 or len(unlabeled)<10:
        raise ValueError("Not enough positives or unlabeled genes")
    values=np.ascontiguousarray(X.to_numpy(dtype=float))
    position={g:i for i,g in enumerate(genes)}
    n_bg=min(len(unlabeled),max(len(positive)*background_ratio,len(positive)+1))
    if matched_background and strata is not None:
        eligible=set(unlabeled)
        pools={s:np.array([g for g in members.index.astype(str) if g in eligible],dtype=object)
               for s,members in strata.groupby(strata)}
        quotas=strata.reindex(positive).fillna("missing").astype(str).value_counts()
    predictions=[]; coefficients=[]; audits=[]
    for seed in seeds:
        rng=np.random.default_rng(int(seed))
        if matched_background and strata is not None:
            background=[]; used=set()
            for stratum,count in quotas.items():
                pool=pools.get(stratum,np.array([],dtype=object))
                draw=min(int(count)*background_ratio,len(pool))
                if draw:
                    selected=rng.choice(pool,draw,replace=False).tolist()
                    background.extend(selected);used.update(selected)
            fallback_n=n_bg-len(background)
            if fallback_n:
                remaining=[g for g in unlabeled if g not in used]
                background.extend(rng.choice(remaining,min(fallback_n,len(remaining)),replace=False).tolist())
            background=background[:n_bg]
        else:
            background=rng.choice(unlabeled,n_bg,replace=False).tolist()
            fallback_n=0
        assert not set(background)&all_positive_genes
        train_rows=[position[g] for g in positive+background]
        y=np.array([1]*len(positive)+[0]*len(background))
        model_cfg=cfg["model"]
        requested=int(model_cfg.get("max_iter",10000))
        converged=False
        initial_normalization=bool(model_cfg.get("normalize_class_weights_initially",False))
        budgets=[1,2,5,10] if initial_normalization else [1,1,2,5]
        for retry,multiplier in enumerate(budgets):
            limit=requested*multiplier
            if penalty=="elasticnet":
                model=fp.make_elastic_net_model(float(model_cfg.get("C",1)),float(model_cfg.get("l1_ratio",.5)),limit,str(model_cfg.get("solver","saga")))
            else:
                model=fp.make_l2_model(float(model_cfg.get("C",1)),limit)
            if "tol" in model_cfg:
                model.set_params(clf__tol=float(model_cfg["tol"]))
            # If standard class weighting causes unstable SAGA steps, use an
            # exactly equivalent objective with weights divided by their maximum
            # and C multiplied by that maximum. Preserve observations and seed.
            normalized = (initial_normalization or retry > 0) and penalty == "elasticnet"
            if normalized:
                weights=np.where(y==1,len(y)/(2*len(positive)),len(y)/(2*len(background)))
                weight_scale=float(weights.max())
                model.set_params(clf__class_weight=None,clf__C=float(model_cfg.get("C",1))*weight_scale)
            else:
                weights=None;weight_scale=1.0
            with warnings.catch_warnings(record=True) as caught, threadpool_limits(limits=1):
                warnings.simplefilter("always",ConvergenceWarning)
                if normalized:
                    model.fit(values[train_rows],y,clf__sample_weight=weights/weight_scale)
                else:
                    model.fit(values[train_rows],y)
            iterations=int(model.named_steps["clf"].n_iter_[0])
            failed=any(issubclass(w.category,ConvergenceWarning) for w in caught)
            converged=not failed and iterations<limit
            audits.append({"seed":int(seed),"attempt":retry+1,"n_positive":len(positive),
                "n_background":len(background),"background_ratio":background_ratio,
                "matched_background":int(matched_background and strata is not None),
                "penalty":penalty,"C":float(model_cfg.get("C",1)),
                "l1_ratio":float(model_cfg.get("l1_ratio",.5)) if penalty=="elasticnet" else np.nan,
                "iterations":iterations,"max_iter_used":limit,"converged":int(converged),
                "used_in_ensemble":int(converged),"fallback_background_count":fallback_n,
                "weight_normalization_factor":weight_scale,"convergence_retry_policy":"same observations/seed/objective; normalized weights and scaled C, then increased iteration budget"})
            if converged:break
        if not converged:raise RuntimeError(f"Seed {seed} failed to converge on unchanged sample")
        scaler=model.named_steps["scale"];clf=model.named_steps["clf"]
        weight=clf.coef_[0]/scaler.scale_
        intercept=clf.intercept_[0]-scaler.mean_@weight
        with threadpool_limits(limits=1):
            predictions.append(expit(values@weight+intercept))
        coefficients.append(pd.DataFrame({"feature":X.columns,"coefficient":clf.coef_[0],"seed":int(seed)}))
    pred=np.column_stack(predictions)
    score=pred.mean(axis=1)
    scores=pd.DataFrame({"gene_symbol":genes,"score_mean":score,
        "score_sd":pred.std(axis=1,ddof=1) if len(seeds)>1 else np.zeros(len(genes))})
    if cfg.get("compute_repeat_rank_frequency",True):
        ranks=pd.DataFrame(pred).rank(ascending=False,method="min")
        for k in [20,50,100]:scores[f"top{k}_frequency"]=(ranks<=k).mean(axis=1).to_numpy()
    else:
        for k in [20,50,100]:scores[f"top{k}_frequency"]=np.nan
    scores["rank"]=scores.score_mean.rank(ascending=False,method="min").astype(int)
    audit_frame=pd.DataFrame(audits)
    audit_root=Path(cfg.get("_run_output", Path(__file__).resolve().parents[2]))/"logs/convergence"
    audit_root.mkdir(parents=True,exist_ok=True)
    record={"stage":os.environ.get("REVISION_STAGE","unscoped_earlier_run"),
        "run_id":os.environ.get("REVISION_RUN_ID","unscoped_earlier_run"),
        "cancer":cfg.get("cancer_code"),"positive_genes":sorted(positive_genes),
        "model_columns":list(X.columns),"C":float(model_cfg.get("C",1)),
        "l1_ratio":float(model_cfg.get("l1_ratio",.5)),"background_ratio":background_ratio,
        "tol":float(model_cfg.get("tol",1e-2)),"initial_weight_normalization":initial_normalization,
        "matched_background":bool(matched_background),"n_requested_fits":len(seeds),
        "n_attempts":len(audits),"n_accepted_fits":int(audit_frame.used_in_ensemble.sum()),
        "seeds":list(map(int,seeds)),"numerical_retries":[a for a in audits if a["attempt"]>1 or not a["converged"]]}
    with (audit_root/f"ensembles_{os.getpid()}.jsonl").open("a",encoding="utf8") as stream:
        stream.write(json.dumps(record)+"\n")
    return scores,pd.concat(coefficients,ignore_index=True),audit_frame
