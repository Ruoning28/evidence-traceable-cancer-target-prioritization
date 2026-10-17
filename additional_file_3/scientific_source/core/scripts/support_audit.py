"""Descriptive gate, coverage, threshold and benchmark provenance audits."""
from __future__ import annotations
import json,itertools
import numpy as np
import pandas as pd
from run_revision import ROOT,prepare,fp,save

BLOCKS={"patient":["patient_alteration_score"],"protein":["proteomic_support_score"],
"dependency":["dependency_score"],"drug":["chemical_modulation_score","clinical_precedence_score"],
"external":["external_reproducibility_score","opentargets_overall_score","geo_external_reproducibility_score"]}


def review(matrix,scores,cfg):
 m=matrix.set_index("gene_symbol").join(scores.set_index("gene_symbol"))
 r=cfg["review"]
 tc=m[fp.CLINICAL_SCORE].quantile(1-r["high_score_percentile"])
 fc=m[fp.DISCOVERY_SCORE].quantile(1-r["high_score_percentile"])
 dc=m.dependency_score.quantile(1-r["dependency_top_percentile"])
 rc=m.disease_relevance_score.quantile(1-r["disease_top_percentile"])
 override=m.clinical_precedence_score.ge(r["clinical_precedence_override_min"])|m.opentargets_overall_score.ge(r["opentargets_override_min"])
 conditions=pd.DataFrame({"translational_score_pass":m[fp.CLINICAL_SCORE].ge(tc),
 "functional_score_pass":m[fp.DISCOVERY_SCORE].ge(fc),"dependency_half_reference_pass":m.dependency_score.ge(dc*.5),
 "dependency_full_reference_pass":m.dependency_score.ge(dc),"disease_half_reference_pass":m.disease_relevance_score.ge(rc*.5),
 "safety_acceptability_pass":m.safety_factor.ge(r["safety_factor_min"]),
 "high_safety_risk":m.safety_risk_score.ge(r["safety_risk_high_threshold"]),
 "translational_actionability_pass":m.druggability_score.ge(.1),"limited_actionability":m.druggability_score.lt(r["limited_actionability_druggability_max"]),
 "translational_override":override},index=m.index)
 t=conditions.translational_score_pass&conditions.translational_actionability_pass&conditions.safety_acceptability_pass&(conditions.dependency_half_reference_pass|override)
 limited=conditions.dependency_full_reference_pass&conditions.limited_actionability
 f=conditions.functional_score_pass&conditions.dependency_half_reference_pass&conditions.disease_half_reference_pass&conditions.safety_acceptability_pass
 biomarker=m.disease_relevance_score.ge(rc)&m.dependency_score.lt(dc*.5)
 cls=pd.Series(np.select([conditions.high_safety_risk,t,limited,f,biomarker],["IV","I","III","II","V"],default="VI"),index=m.index,name="class")
 return cls,conditions


def run(cancer):
 matrix,labels,benchmarks,roles,features,cfg,out=prepare(cancer,write=False)
 scores=pd.read_csv(out/"framework_priority_scores.tsv",sep="\t")
 m=matrix.set_index("gene_symbol")
 cls,conditions=review(matrix,scores,cfg)
 stored=pd.read_csv(out/"framework_review_classes.tsv",sep="\t").set_index("gene_symbol").review_class.str.extract(r"Class ([IVX]+):",expand=False)
 assert cls.equals(stored.rename("class")),"Vectorized gate audit differs from executable review"
 support=sum(m[cols].fillna(0).gt(0).any(axis=1).astype(int) for cols in BLOCKS.values())
 vi=cls.eq("VI");high=conditions.translational_score_pass|conditions.functional_score_pass
 reason=pd.Series("multiple support blocks below both score gates",index=m.index)
 reason.loc[support.eq(0)]="no positive support block"
 reason.loc[support.eq(1)]="one positive support block"
 reason.loc[support.ge(2)&high]="multiple support blocks and high score; review constraints"
 audit=conditions.copy();audit["review_class"]=cls;audit["positive_support_block_count"]=support
 audit["class_vi_reason"]=reason.where(vi,"")
 save(audit.reset_index(),out/"revision_gate_failure_gene_audit.tsv")
 partition=reason[vi].value_counts().rename_axis("reason").reset_index(name="count")
 partition["fraction_of_class_vi"]=partition["count"]/int(vi.sum())
 partition["note"]="Positive-support proxy; zero may represent observed zero or absent annotation."
 save(partition,out/"revision_class_vi_gate_partition.tsv")
 failure=[]
 for g in conditions.columns:
  failure.append({"cancer":cancer,"gate":g,"class_vi_count":int(vi.sum()),"class_vi_gate_true_count":int(conditions.loc[vi,g].sum()),"note":"Overlapping descriptions, not an additive partition"})
 save(pd.DataFrame(failure),out/"revision_class_vi_overlapping_gate_failures.tsv")
 moves=[];sizes=[];primary_ii=set(cls.index[cls.isin(["I","II"])])
 for pct,safe in itertools.product([.01,.02,.05],[.3,.4,.5]):
  cfg2=json.loads(json.dumps(cfg));cfg2["review"].update(high_score_percentile=pct,safety_factor_min=safe)
  alternate,_=review(matrix,scores,cfg2)
  members=set(alternate.index[alternate.isin(["I","II"])])
  counts=alternate.value_counts()
  sizes.append({"cancer":cancer,"high_score_percentile":pct,"safety_factor_min":safe,
    "class_I_II_count":len(members),"primary_I_II_retained":len(members&primary_ii),
    "I_II_jaccard_with_primary":len(members&primary_ii)/max(1,len(members|primary_ii)),
    "total_class_changed":int((alternate!=cls).sum()),**{f"class_{x}_count":int(counts.get(x,0)) for x in ["I","II","III","IV","V","VI"]}})
  transition=pd.crosstab(cls,alternate)
  for a in transition.index:
   for b in transition.columns:
    moves.append({"cancer":cancer,"high_score_percentile":pct,"safety_factor_min":safe,"primary_class":a,"alternate_class":b,"count":int(transition.loc[a,b])})
 save(pd.DataFrame(sizes),out/"revision_threshold_candidate_retention.tsv")
 save(pd.DataFrame(moves),out/"revision_threshold_class_transitions.tsv")
 stats=[]
 for feature in cfg["features"]:
  x=pd.to_numeric(m[feature],errors="coerce")
  stats.append({"cancer":cancer,"feature":feature,"n_candidates":len(x),"input_missing":int(x.isna().sum()),"input_zero":int(x.eq(0).sum()),"input_positive":int(x.gt(0).sum()),"constant_nonmissing":x.dropna().nunique()<=1,"minimum":x.min(),"maximum":x.max(),"availability_note":"Input missingness; upstream zero imputation is not reversible from this column alone"})
 save(pd.DataFrame(stats),out/"revision_input_feature_coverage_audit.tsv")
 corr=features.set_index("gene_symbol").corr(method="spearman")
 save(corr.rename_axis("feature").reset_index(),out/"revision_feature_spearman_matrix.tsv")
 rows=[];s=scores.set_index("gene_symbol")
 for axis in fp.AXES:
  values=s[axis.score_column];rank=values.rank(ascending=False,method="first")
  bgenes=benchmarks[axis.name]
  for gene in sorted(bgenes):
   rr=roles[(roles.axis==axis.name)&roles.gene_symbol.eq(gene)&roles.record_role.eq("V")]
   row={"cancer":cancer,"axis":axis.name,"gene_symbol":gene,"rank":int(rank[gene]),"score":float(values[gene]),"top50":bool(rank[gene]<=50),"benchmark_sources":";".join(sorted(set(rr.source))),"fact_independence_guaranteed":False}
   for block,columns in BLOCKS.items():row[block+"_positive_support"]=bool(pd.to_numeric(m.loc[gene,columns],errors="coerce").fillna(0).gt(0).any())
   rows.append(row)
 save(pd.DataFrame(rows),out/"revision_benchmark_gene_ranks_and_source_overlap.tsv")
 print(cancer,"support audits PASS",flush=True)

if __name__=="__main__":
 for cancer in ["LUAD","PDAC"]:run(cancer)
