from pathlib import Path
import ast
import json
import numpy as np
import pandas as pd

import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--reference',type=Path,required=True)
parser.add_argument('--code',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
HERE=args.output.resolve()
HERE.mkdir(parents=True,exist_ok=True)
ROOT=HERE
DATA=args.data.resolve()
CODE=args.code.resolve()/'gate_engine.py'
tree=ast.parse(CODE.read_text(encoding='utf-8-sig'))
nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in {'base_policy','cutoffs','classify'}]
namespace={'np':np,'pd':pd}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(CODE),'exec'),namespace)
get_cutoffs=namespace['cutoffs']
expected={'PDAC':{1:267,2:79,3:1883,4:1184,5:2974,6:34204},
          'LUAD':{1:289,2:15,3:574,4:1184,5:1390,6:32026}}
rows=[]
for cancer in ('PDAC','LUAD'):
    p=DATA/'frozen_inputs'/cancer
    m=pd.read_csv(p/'source_target_evidence_matrix.tsv',sep='\t',low_memory=False).set_index('gene_symbol')
    s=pd.read_csv(DATA/'review_inputs'/cancer/'archived_strict_priority_scores.tsv',sep='\t').set_index('gene_symbol').reindex(m.index)
    cfg=json.loads((p/'analysis_config.json').read_text(encoding='utf-8-sig'))
    policy=namespace['base_policy'](cfg)
    full=namespace['classify'](m,s,policy)
    assert pd.Series(full).value_counts().to_dict()==expected[cancer]
    full_cuts=get_cutoffs(m,s,policy)
    mapping=pd.read_csv(HERE/'identifier_results'/f'{cancer}_record_mapping.tsv',sep='\t',keep_default_na=False)
    mapping['key']=np.where(mapping.hgnc_id.ne(''),mapping.hgnc_id,'UNRESOLVED:'+mapping.gene_symbol)
    mapping['prefer_approved']=mapping.mapping_status.ne('approved_symbol').astype(int)
    # Common representative for both axes, selected without using any score:
    # an approved-symbol record if present, otherwise first frozen input row.
    preferred=mapping.sort_values('prefer_approved',kind='stable').drop_duplicates('key',keep='first')
    for name,subset in [('full_record_universe',mapping),
                        ('collapse_unique_mappings_keep_unresolved',preferred),
                        ('mapped_human_genes_only',preferred[preferred.hgnc_id.ne('')])]:
        ids=subset.gene_symbol.tolist(); mm=m.loc[ids];ss=s.loc[ids]
        original=pd.Series(full,index=m.index).loc[ids].to_numpy()
        for mode in ('historical_full_universe_cutoffs','recomputed_subset_cutoffs'):
            namespace['cutoffs']=(lambda _m,_s,_p:full_cuts) if mode=='historical_full_universe_cutoffs' else get_cutoffs
            classes=namespace['classify'](mm,ss,policy)
            row=dict(cancer=cancer,universe=name,cutoff_mode=mode,n_records=len(ids),
                     representative_selection='Approved-symbol row preferred, otherwise first frozen input row; no score-based selection.',
                     changed_classes_among_retained_records=int((classes!=original).sum()),
                     I_II=int(np.isin(classes,[1,2]).sum()),VI_fraction=float((classes==6).mean()))
            priority_ids=np.asarray(ids)[np.isin(classes,[1,2])]
            priority_map=mapping.set_index('gene_symbol').loc[priority_ids]
            row['I_II_unique_mapped_genes']=int(priority_map.loc[priority_map.hgnc_id.ne(''),'hgnc_id'].nunique())
            row['I_II_unresolved_records']=int(priority_map.hgnc_id.eq('').sum())
            row.update({f'class_{i}':int((classes==i).sum()) for i in range(1,7)})
            row.update({f'cutoff_{k}':v for k,v in namespace['cutoffs'](mm,ss,policy).items()})
            rows.append(row)
    namespace['cutoffs']=get_cutoffs
pd.DataFrame(rows).to_csv(HERE/'identifier_results/review_universe_sensitivity.tsv',sep='\t',index=False)
print(pd.DataFrame(rows)[['cancer','universe','cutoff_mode','n_records','I_II','VI_fraction','changed_classes_among_retained_records']].to_string(index=False))
