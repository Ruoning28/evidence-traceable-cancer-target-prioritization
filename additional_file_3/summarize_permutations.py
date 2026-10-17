"""Apply within-null and joint Holm corrections to the eight reported tests."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

def holm(values):
 p=np.asarray(values,float);order=np.argsort(p,kind='stable');out=np.empty(len(p))
 out[order]=np.minimum(1,np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p)))))
 return out

def summarize(frame):
 out=frame.copy().reset_index(drop=True)
 if len(out)!=8 or out[['cancer','axis','null']].duplicated().any():
  raise ValueError('Expected both nulls for both axes in both cancers (eight tests)')
 out['holm_p_across_all_eight_tests']=holm(out.empirical_p_value)
 out['holm_p_within_null_family']=np.nan
 for _,idx in out.groupby('null').groups.items():
  out.loc[idx,'holm_p_within_null_family']=holm(out.loc[idx,'empirical_p_value'])
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--data',type=Path);p.add_argument('--results',type=Path);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 if bool(a.data)==bool(a.results):p.error('Supply exactly one of --data or --results')
 if a.data:
  frame=pd.read_csv(a.data/'historical_tables/strict_tolerance_revision_permutation_summary.tsv',sep='\t')
 else:
  frames=[]
  for cancer in ['PDAC','LUAD']:
   x=pd.read_csv(a.results/cancer/'tables/revision_permutation_summary.tsv',sep='\t');x['cancer']=cancer;frames.append(x)
  frame=pd.concat(frames,ignore_index=True)
 out=summarize(frame);a.output.parent.mkdir(parents=True,exist_ok=True);out.to_csv(a.output,sep='\t',index=False)
 print('Holm correction completed for eight tests')
