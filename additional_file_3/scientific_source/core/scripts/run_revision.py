"""Frozen-input loading and execution context for the study analyses."""
from pathlib import Path
import copy
import json
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from framework import pipeline as fp
import fast_model
fast_model.install(fp)

DATA=None
OUTPUT=None
PROFILE='strict'

def configure(data,output,profile='strict'):
 global DATA,OUTPUT,PROFILE
 DATA=Path(data).resolve();OUTPUT=Path(output).resolve();PROFILE=profile
 OUTPUT.mkdir(parents=True,exist_ok=True)

def save(df,path):
 path.parent.mkdir(parents=True,exist_ok=True)
 df.to_csv(path,sep='\t',index=False)

def prepare(cancer,write=True):
 if DATA is None:raise RuntimeError('Configure the data and output directories before preparing inputs')
 folder=DATA/'frozen_inputs'/cancer
 read=lambda name:pd.read_csv(folder/name,sep='\t',keep_default_na=False,low_memory=False)
 matrix=read('source_target_evidence_matrix.tsv')
 # Convert numeric source fields without treating literal identifier keys as NA.
 for col in matrix:
  if col!='gene_symbol':
   parsed=pd.to_numeric(matrix[col],errors='coerce')
   if parsed.notna().any() or matrix[col].eq('').all():matrix[col]=parsed
 labels=read('framework_development_positive_sets.tsv')
 roles=read('framework_record_role_manifest.tsv')
 features=read('framework_feature_matrix.tsv')
 for col in features:
  if col!='gene_symbol':features[col]=pd.to_numeric(features[col],errors='raise')
 cfg=copy.deepcopy(json.loads((folder/'analysis_config.json').read_text(encoding='utf-8-sig')))
 cfg['_run_output']=str(OUTPUT)
 if PROFILE=='reference':
  cfg['model'].update(tol=1e-2,normalize_class_weights_initially=False)
 elif PROFILE!='strict':raise ValueError(PROFILE)
 cfg['compute_repeat_rank_frequency']=True
 out=OUTPUT/cancer/'tables';out.mkdir(parents=True,exist_ok=True)
 benchmark={axis.name:set(roles.loc[(roles.record_role=='V')&(roles.axis==axis.name),'gene_symbol']) for axis in fp.AXES}
 if write:
  for name,table in [('framework_development_positive_sets.tsv',labels),('framework_record_role_manifest.tsv',roles),('framework_feature_matrix.tsv',features)]:save(table,out/name)
  save(fp.build_background_strata(matrix).reset_index(),out/'framework_background_strata.tsv')
  save(roles[roles.record_role.eq('V')],out/'framework_validation_benchmarks.tsv')
  (out.parent/'analysis_config.json').write_text(json.dumps(cfg,indent=2),encoding='utf8')
 return matrix,labels,benchmark,roles,features,cfg,out

def scores_for(cancer,out):
 p=out/'framework_priority_scores.tsv'
 if not p.exists():
  if PROFILE!='strict':raise FileNotFoundError('Run primary with --profile reference before reference-profile diagnostics')
  return pd.read_csv(DATA/'review_inputs'/cancer/'archived_strict_priority_scores.tsv',sep='\t',float_precision='round_trip')
 return pd.read_csv(p,sep='\t',float_precision='round_trip')

def initialize_diagnostics(cancer):
 m,l,b,r,f,c,o=prepare(cancer)
 s=scores_for(cancer,o)
 save(s,o/'framework_priority_scores.tsv')
 fp.assign_review_classes(m,s,c,o)
 return m,l,b,r,f,c,o,s
