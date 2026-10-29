"""Reproduce original strict-tolerance Top-50 recovery from frozen matrices."""
from pathlib import Path
import os,sys
if os.name=='nt':
    lib=Path(sys.prefix)/'Library/bin'
    if lib.exists():
        os.environ['PATH']=str(lib)+os.pathsep+os.environ.get('PATH','')
        DLL_HANDLE=os.add_dll_directory(str(lib))
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'scientific_source/core/src'))
sys.path.insert(0,str(HERE/'scientific_source/core/scripts'))
import argparse,json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from framework import pipeline as fp
import fast_model

def explicit_elastic_net(C,l1_ratio,max_iter,solver='saga'):
    # In the specified 1.7.2 runtime penalty is explicit; this preserves the
    # elastic-net objective implemented by l1_ratio in the original 1.9 runtime.
    return Pipeline([('scale',StandardScaler()),('clf',LogisticRegression(
      penalty='elasticnet',solver=solver,C=C,l1_ratio=l1_ratio,max_iter=max_iter,
      tol=1e-6,class_weight='balanced',random_state=0))])
fp.make_elastic_net_model=explicit_elastic_net

def run(data,out):
    out.mkdir(parents=True,exist_ok=True);rows=[]
    (Path(fast_model.__file__).resolve().parents[2]/'logs/convergence').mkdir(parents=True,exist_ok=True)
    for cancer in ['PDAC','LUAD']:
        f=data/'frozen_inputs'/cancer
        read=lambda n:pd.read_csv(f/n,sep='\t',keep_default_na=False)
        x=read('framework_feature_matrix.tsv').set_index('gene_symbol')
        labels=read('framework_development_positive_sets.tsv').set_index('gene_symbol')
        roles=read('framework_record_role_manifest.tsv')
        strata=read('framework_background_strata.tsv').set_index('gene_symbol').background_stratum
        cfg=json.loads((f/'analysis_config.json').read_text(encoding='utf-8-sig'))
        cfg['compute_repeat_rank_frequency']=False
        all_bench=set(roles.loc[roles.record_role=='V','gene_symbol'])
        for axis,label in [('translational','clinical_positive'),('functional','discovery_positive')]:
            pos=set(labels.index[pd.to_numeric(labels[label]).eq(1)])
            columns=fp.axis_feature_columns(x,cfg,axis)
            values=x[columns].apply(pd.to_numeric)
            scores,_,audit=fast_model.fit_sampled_background_model(values,pos,pos|all_bench,
              cfg['random_seeds'][:30],cfg,strata,True,8)
            benchmark=set(roles.loc[(roles.record_role=='V')&(roles.axis==axis),'gene_symbol'])
            # Scores are RAM-only. Only aggregate benchmark recovery is exported.
            top=set(scores.nlargest(50,'score_mean').gene_symbol)
            hits=len(top&benchmark)
            rows.append(dict(cancer=cancer,axis=axis,n_benchmark=len(benchmark),top50_hits=hits,
              top50_recovery=hits/len(benchmark),n_accepted_fits=int(audit.used_in_ensemble.sum())))
            print('Primary recovery',cancer,axis,hits,len(benchmark),flush=True)
    expected={('PDAC','translational'):4,('PDAC','functional'):4,('LUAD','translational'):4,('LUAD','functional'):2}
    assert all(r['top50_hits']==expected[(r['cancer'],r['axis'])] and r['n_accepted_fits']==30 for r in rows)
    pd.DataFrame(rows).to_csv(out/'primary_recovery.tsv',sep='\t',index=False)
    (out/'primary_checks.json').write_text(json.dumps(dict(all_matched=True,n_settings=4,n_fits=120,score_vectors_exported=False),indent=2),encoding='utf8')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.data,a.output)
