"""Reproduce aggregate measurement-conditional review for three backends."""
from pathlib import Path
import argparse,json,sys,os
if os.name=='nt':
    lib=Path(sys.prefix)/'Library/bin'
    if lib.exists():DLL_HANDLE=os.add_dll_directory(str(lib));os.environ['PATH']=str(lib)+os.pathsep+os.environ.get('PATH','')
import numpy as np
import pandas as pd
import gate_engine as gate
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'scientific_source/core/src'))
from framework import pipeline as fp
BACKENDS=['archived_strict_elastic_net','direction_corrected_mean_same_information','native_rank_aggregation_trans_plus_dependency_func']
PROTECTED={'IL1A','SLC5A2','ACOT8','SETD1A','PRDX1','ERO1A','UQCRC1','GOT1'}

def index(path):return pd.read_csv(path,sep='\t',float_precision='round_trip',low_memory=False).set_index('gene_symbol')

def run(data,out):
    out.mkdir(parents=True,exist_ok=True);rows=[]
    for cancer in ['PDAC','LUAD']:
        folder=data/'frozen_inputs'/cancer;m=index(folder/'source_target_evidence_matrix.tsv')
        f=index(folder/'framework_feature_matrix.tsv').reindex(m.index)
        cfg=json.loads((folder/'analysis_config.json').read_text(encoding='utf-8-sig'))
        q=data/'review_inputs'/cancer;primary=index(q/'archived_strict_priority_scores.tsv').reindex(m.index)
        lo=np.zeros(len(m));hi=np.ones(len(m));observed=np.zeros(len(m),bool)
        with np.load(q/'measurement_intervals.npz',allow_pickle=False) as a:
            ix=m.index.get_indexer(a['genes']);assert (ix>=0).all()
            lo[ix]=a['risk_lower'];hi[ix]=a['risk_upper'];observed[ix]=a['constraint_observed']
        cols=['translational_precedence_score','functional_discovery_score'];mean={};native={}
        for axis,col in zip(['translational','functional'],cols):
            excluded=set(cfg['evaluation']['axes'][axis]['feature_exclusions'])
            fields=[c for c in f if c not in excluded and not c.endswith('_missing')]
            votes=[1-f[c].rank(pct=True) if c=='safety_risk_score' else f[c].rank(pct=True) for c in fields]
            mean[col]=pd.concat(votes,axis=1).mean(axis=1)
            native[col]=fp._rank_aggregation(m,['clinical_precedence_score','chemical_modulation_score','druggability_score','opentargets_overall_score']) if axis=='translational' else m.dependency_score.fillna(0)
        policies=gate.policies(cfg);assert len(policies)==31
        eligible=~m.index.isin(PROTECTED)
        for backend,scores in zip(BACKENDS,[primary[cols],pd.DataFrame(mean,index=m.index),pd.DataFrame(native,index=m.index)]):
            reference=gate.classify(m,scores,policies[0]['parameters'])
            bits=[gate.completion_bits(m,scores,p['parameters'],lo,hi) for p in policies]
            union=np.bitwise_or.reduce(bits);ref=bits[0];nclasses=np.array([int(b).bit_count() for b in union])
            priority=np.isin(reference,[1,2])&eligible
            robust=(ref>0)&((ref&np.uint8(0b111100))==0)&eligible
            rows.append(dict(cancer=cancer,backend=backend,n_candidates=int(eligible.sum()),
              n_archived_reference_priority=int(priority.sum()),n_reference_priority_all_completions=int(robust.sum()),
              n_archived_priority_not_completion_robust=int((priority&~robust).sum()),
              n_archived_priority_with_unobserved_constraint=int((priority&~observed).sum()),
              n_single_class_all31policies=int(((nclasses==1)&eligible).sum()),
              n_multi_class_all31policies=int(((nclasses>1)&eligible).sum())))
        print('Conditional review',cancer,flush=True)
    result=pd.DataFrame(rows)
    expected=pd.read_csv(data/'studies/measurement_review_20261012_v1/review_summary.tsv',sep='\t')
    compare=result.merge(expected,on=['cancer','backend'],suffixes=('_actual','_expected'))
    for k in result.columns:
        if k not in ['cancer','backend']:assert (compare[k+'_actual']==compare[k+'_expected']).all(),k
    result.to_csv(out/'review_summary.tsv',sep='\t',index=False)
    (out/'review_checks.json').write_text(json.dumps(dict(all_matched=True,n_rows=6,n_policies=31,full_class_vectors_exported=False),indent=2),encoding='utf8')
    print('Review reproduction passed',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.data,a.output)
