"""Construct the model feature matrices and manifests from supplied source matrices."""
from pathlib import Path
import argparse
import sys
import os
if os.name=='nt':
 lib=Path(sys.prefix)/'Library/bin'
 DLL_HANDLE=os.add_dll_directory(str(lib)) if lib.exists() else None
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'scientific_source/core/scripts'))
import run_revision as rr

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();rr.configure(a.data,a.output)
 for cancer in ['PDAC','LUAD']:
  m,l,b,r,f,c,o=rr.prepare(cancer)
  features,manifest=rr.fp.build_feature_matrix(m,c)
  rr.save(features,o/'framework_feature_matrix.tsv');rr.save(manifest,o/'framework_feature_manifest.tsv')
  print(cancer,len(features),'records')
