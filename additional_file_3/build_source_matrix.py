"""Construct target-level evidence matrices from the retained public-source layout."""
from pathlib import Path
import argparse
import sys
import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'scientific_source/reference'))

def pdac_matrix(root,output):
 import v4_end_to_end as source
 patient=source.build_patient_and_protein_features(root)
 tcga=source.build_tcga_features(root)
 geo,_=source.build_geo_validation_features(root)
 dependency=source.build_dependency_features(root)
 drug,_=source.build_druggability_features(root)
 common=dependency['is_common_essential'] if 'is_common_essential' in dependency else None
 safety=source.build_safety_features(root,common)
 known=source.load_known_targets(root)
 genes=set().union(*(set(x.index) for x in [patient,tcga,geo,dependency,drug] if not x.empty))
 genes|={g for group in known.values() for g in group}
 literature=source.extract_literature_benchmark(root,genes)
 if not literature.empty:genes|=set(literature.gene_symbol)
 literature_counts=pd.DataFrame(index=sorted(genes))
 if not literature.empty:
  literature_counts['literature_benchmark_source_count']=literature.groupby('gene_symbol').source.nunique()
  literature_counts['literature_benchmark_annotation_score']=literature_counts.literature_benchmark_source_count.fillna(0).clip(0,2)/2
 else:
  literature_counts['literature_benchmark_source_count']=0
  literature_counts['literature_benchmark_annotation_score']=0
 matrix=pd.concat([x.reindex(sorted(genes)) for x in [patient,tcga,geo,dependency,drug,safety,literature_counts]],axis=1).fillna(0)
 matrix['external_reproducibility_score']=matrix.geo_external_reproducibility_score.fillna(0).clip(0,1) if 'geo_external_reproducibility_score' in matrix else 0.
 matrix=matrix[matrix.index.astype(str).map(bool)]
 columns=[c for c in ['patient_alteration_score','tcga_patient_score'] if c in matrix]
 if columns:matrix['patient_alteration_score']=source._noisy_or_frame(matrix,columns)
 scored=source.build_scores(matrix)
 output.mkdir(parents=True,exist_ok=True)
 source._save_indexed_table(scored,output/'pdac_target_evidence_matrix.tsv')

def luad_matrix(root,shared,output):
 import build_luad_inputs as source
 source.LUAD_ROOT=root
 source.RAW=root/'data/raw'
 source.PROCESSED=output
 source.LEGACY_ROOT=shared
 source.main()

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--cancer',choices=['PDAC','LUAD'],required=True)
 p.add_argument('--source-root',type=Path,required=True)
 p.add_argument('--shared-pdac-root',type=Path)
 p.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 if a.cancer=='PDAC':pdac_matrix(a.source_root.resolve(),a.output.resolve())
 else:
  if a.shared_pdac_root is None:p.error('LUAD uses --shared-pdac-root for shared annotations and DepMap inputs')
  luad_matrix(a.source_root.resolve(),a.shared_pdac_root.resolve(),a.output.resolve())
