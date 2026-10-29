from pathlib import Path
import pandas as pd,numpy as np,json
import sys
A=Path(__file__).resolve().parent
D=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else A.parent.parent/'additional_file_2'
O=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else A/'outputs'
O.mkdir(parents=True,exist_ok=True)
gene_tables=[];summaries=[];sensitivity=[];background=[]
for cancer in ['PDAC','LUAD']:
    f=pd.read_csv(D/f'frozen_inputs/{cancer}/source_target_evidence_matrix.tsv',sep='\t',low_memory=False).set_index('gene_symbol')
    score=pd.read_csv(D/f'review_inputs/{cancer}/archived_strict_priority_scores.tsv',sep='\t').set_index('gene_symbol')
    b=pd.read_csv(D/f'frozen_inputs/{cancer}/framework_validation_benchmarks.tsv',sep='\t')
    g=pd.read_csv(D/f'historical_tables/{cancer}/revision_benchmark_gene_ranks_and_source_overlap.tsv',sep='\t')
    # Existing gene-level support flags describe transformed evidence, not identity of underlying facts.
    for axis,grp in g.groupby('axis',sort=False):
        x={'cancer':cancer,'axis':axis,'benchmark_genes':len(grp)}
        for col in [c for c in g if c.endswith('_positive_support')]:x[col]=int(grp[col].sum())
        summaries.append(x)
    gene_tables.append(g)
    s=pd.read_csv(D/f'historical_tables/{cancer}/framework_model_sensitivity.tsv',sep='\t')
    s=s[s.parameter.isin(['background_ratio','matched_background'])].copy();s.insert(0,'cancer',cancer)
    s['benchmark_genes']=s.axis.map({'translational':11 if cancer=='PDAC' else 7,'functional':30 if cancer=='PDAC' else 36})
    s['top50_hits']=(s.top50_benchmark_recovery*s.benchmark_genes).round().astype(int)
    background.append(s)
    if cancer=='LUAD':
        full=set(b.loc[b.axis=='functional','gene_symbol']);dep=set(b.loc[(b.axis=='functional')&b.source.eq('depmap_luad_high_dependency'),'gene_symbol'])
        assert len(full)==36 and len(dep)==15
        methods={'elastic_net':score.functional_discovery_rank,'dependency_only':f.dependency_score.fillna(0).rank(ascending=False,method='min')}
        for method,ranks in methods.items():
            for name,genes in [('full',full),('DepMap_origin',dep),('other_origins',full-dep)]:
                for k in [20,50,100,200]:
                    hits=sorted(x for x in genes if ranks.loc[x]<=k)
                    sensitivity.append(dict(cancer=cancer,axis='functional',method=method,benchmark_subset=name,k=k,benchmark_n=len(genes),hits=len(hits),recovery=len(hits)/len(genes),precision=len(hits)/k,hit_genes=';'.join(hits),ranking_scope='unchanged full candidate universe'))
pd.concat(gene_tables,ignore_index=True).to_csv(O/'benchmark_gene_source_support.tsv',sep='\t',index=False)
pd.DataFrame(summaries).to_csv(O/'benchmark_source_support_summary.tsv',sep='\t',index=False)
pd.DataFrame(sensitivity).to_csv(O/'LUAD_functional_source_subset_sensitivity.tsv',sep='\t',index=False)
pd.concat(background,ignore_index=True).to_csv(O/'background_design_sensitivity.tsv',sep='\t',index=False)
print(pd.DataFrame(summaries).to_string(index=False))
print(pd.DataFrame(sensitivity).query('k==50').to_string(index=False))
print(pd.concat(background)[['cancer','variant','axis','top50_hits','benchmark_genes','top50_overlap_with_primary']].to_string(index=False))
