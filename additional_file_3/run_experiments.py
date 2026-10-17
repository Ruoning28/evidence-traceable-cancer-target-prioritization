"""Run the manuscript's primary models and reviewer-requested diagnostics."""
from pathlib import Path
import argparse
import copy
import importlib.metadata
import json
import os
import sys
import datetime

for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
 os.environ.setdefault(key,'1')
DLL_HANDLES=[]
if os.name=='nt':
 lib=Path(sys.prefix)/'Library/bin'
 if lib.exists():DLL_HANDLES.append(os.add_dll_directory(str(lib)))
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'scientific_source/core/scripts'))
sys.path.insert(0,str(HERE/'scientific_source/core/src'))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import run_revision as rr
from framework import pipeline as fp

STAGES=['prepare','primary','enrichment','oof','source-out','comparators','legacy-comparators',
        'ablation','perturbation','permutation','background','regularization','thresholds','support']

def audit_inputs(cancer):
 m,l,b,r,f,c,o=rr.prepare(cancer)
 rebuilt,manifest=fp.build_feature_matrix(m,c)
 np.testing.assert_allclose(rebuilt[f.columns[1:]].to_numpy(float),f.iloc[:,1:].to_numpy(float),rtol=0,atol=2e-15)
 saved=pd.read_csv(rr.DATA/'frozen_inputs'/cancer/'framework_background_strata.tsv',sep='\t',keep_default_na=False)
 strata=fp.build_background_strata(m)
 assert strata.reindex(saved.gene_symbol).tolist()==saved.background_stratum.tolist()
 rows=[]
 for a in fp.AXES:
  pos=set(l.loc[pd.to_numeric(l[a.label_column]).eq(1),'gene_symbol'])
  assert not pos&b[a.name]
  cols=fp.axis_feature_columns(f,c,a.name)
  excluded=set(c['evaluation']['axes'][a.name]['feature_exclusions'])
  assert not set(cols)&(excluded|{x+'_missing' for x in excluded})
  rows.append(dict(cancer=cancer,axis=a.name,n_candidates=len(m),n_positive=len(pos),n_benchmark=len(b[a.name]),n_features=len(cols),n_repeats=c['background_sampling']['n_repeats'],tol=c['model']['tol']))
 rr.save(manifest,o/'framework_feature_manifest.tsv')
 rr.save(pd.DataFrame(rows),o/'input_checks.tsv')
 print(json.dumps(rows),flush=True)

def parameter_diagnostic(cancer,kind):
 m,l,b,r,f,c,o,s=rr.initialize_diagnostics(cancer)
 X=f.set_index('gene_symbol');lab=l.set_index('gene_symbol');base=s.set_index('gene_symbol')
 all_bench=set().union(*b.values());rows=[]
 variants=[]
 if kind=='background':
  variants=[('background_ratio_'+str(q),q,True,{}) for q in c['sensitivity']['background_ratio_values']]
  variants += [('matched_background_'+str(int(v)),c['background_sampling']['background_ratio'],bool(v),{}) for v in c['sensitivity']['matched_background_values']]
 else:
  for key in ['C','l1_ratio']:
   for value in c['sensitivity'][key+'_values']:variants.append((key+'_'+str(value),8,True,{key:value}))
 for name,ratio,matched,changes in variants:
  cc=copy.deepcopy(c);cc['model'].update(changes)
  for a in fp.AXES:
   pos=set(lab.index[pd.to_numeric(lab[a.label_column]).eq(1)])
   values,_,audit=fp.fit_sampled_background_model(X[fp.axis_feature_columns(f,cc,a.name)],pos,pos|all_bench,cc['random_seeds'][:int(cc['sensitivity']['n_repeats'])],cc,fp.build_background_strata(m) if matched else None,matched,int(ratio))
   score=values.set_index('gene_symbol').score_mean;reference=base[a.score_column].reindex(score.index)
   top=set(score.nlargest(50).index)
   rows.append(dict(cancer=cancer,axis=a.name,variant=name,profile=rr.PROFILE,
    C=cc['model']['C'],l1_ratio=cc['model']['l1_ratio'],background_ratio=ratio,matched_background=int(matched),
    spearman_with_primary=spearmanr(reference,score).statistic,top50_overlap_with_primary=len(top&set(reference.nlargest(50).index))/50,
    top50_benchmark_recovery=len(top&b[a.name])/len(b[a.name]),n_accepted_fits=int(audit.used_in_ensemble.sum())))
 rr.save(pd.DataFrame(rows),o/('reproduced_'+kind+'_sensitivity.tsv'))

def run_stage(cancer,stage,jobs):
 if stage=='prepare':return audit_inputs(cancer)
 if stage=='primary':
  m,l,b,r,f,c,o=rr.prepare(cancer)
  s=fp.train_dual_axis_models(f,l,b,m,c,o)
  fp.run_independent_benchmark_validation(s,b,c,o)
  cl=fp.assign_review_classes(m,s,c,o);fp.build_evidence_cards(m,s,cl,o)
  return
 if stage in ['background','regularization']:return parameter_diagnostic(cancer,stage)
 m,l,b,r,f,c,o,s=rr.initialize_diagnostics(cancer)
 if stage=='enrichment':return fp.run_independent_benchmark_validation(s,b,c,o)
 if stage=='oof':return fp.run_out_of_fold_validation(f,l,b,r,m,c,o)
 if stage=='source-out':return fp.run_leave_source_out(f,r,l,b,m,c,o)
 if stage=='legacy-comparators':
  cl=pd.read_csv(o/'framework_review_classes.tsv',sep='\t')
  return fp.run_comparators(f,m,l,b,s,cl,c,o)
 if stage in ['support','thresholds']:
  import support_audit
  return support_audit.run(cancer)
 import revision_analyses as ra
 if stage=='comparators':return ra.fair_comparisons(cancer)
 if stage=='ablation':return ra.fixed_ablation(cancer)
 if stage=='perturbation':return ra.perturbation(cancer,jobs)
 if stage=='permutation':return ra.permutation(cancer,jobs)
 raise ValueError(stage)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--stage',choices=STAGES+['all'],required=True)
 p.add_argument('--cancer',choices=['PDAC','LUAD','both'],default='both')
 p.add_argument('--profile',choices=['strict','reference'],default='strict');p.add_argument('--jobs',type=int,default=1)
 a=p.parse_args();rr.configure(a.data,a.output,a.profile)
 record=dict(stage=a.stage,profile=a.profile,data=str(a.data.resolve()),started_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
  python=sys.version,packages={x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','scikit-learn','joblib','threadpoolctl']})
 name='run_'+a.stage.replace('-','_')+'_'+a.profile+'.json'
 (rr.OUTPUT/name).write_text(json.dumps(record,indent=2),encoding='utf8')
 # Primary already writes recovery; support already includes threshold tables.
 stages=[s for s in STAGES if s not in {'enrichment','thresholds'}] if a.stage=='all' else [a.stage]
 for cancer in ['PDAC','LUAD'] if a.cancer=='both' else [a.cancer]:
  for stage in stages:
   os.environ['REVISION_STAGE']=stage+'_'+cancer;os.environ['REVISION_RUN_ID']=record['started_at']
   print('START',cancer,stage,flush=True);run_stage(cancer,stage,a.jobs);print('DONE',cancer,stage,flush=True)
 record['completed_at']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 (rr.OUTPUT/name).write_text(json.dumps(record,indent=2),encoding='utf8')

if __name__=='__main__':main()
