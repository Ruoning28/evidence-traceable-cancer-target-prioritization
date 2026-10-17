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
CODE=args.code.resolve()/'scientific_source/core/src/framework/pipeline.py'
tree = ast.parse(CODE.read_text(encoding='utf-8-sig'))
names = {'_numeric', 'build_background_strata', '_matched_background_sample'}
functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
assert len(functions) == 3
namespace = {'np': np, 'pd': pd}
exec(compile(ast.Module(body=functions, type_ignores=[]), str(CODE), 'exec'), namespace)
sample = namespace['_matched_background_sample']
summaries, selected_aliases, selection_rows = [], [], []
for cancer in ('PDAC','LUAD'):
    p = DATA / 'frozen_inputs' / cancer
    matrix = pd.read_csv(p/'source_target_evidence_matrix.tsv', sep='\t', low_memory=False)
    features = pd.read_csv(p/'framework_feature_matrix.tsv', sep='\t')
    genes = features.gene_symbol.astype(str).tolist()
    roles = pd.read_csv(p/'framework_record_role_manifest.tsv', sep='\t', keep_default_na=False)
    strata = pd.read_csv(p/'framework_background_strata.tsv', sep='\t').set_index('gene_symbol').background_stratum
    recreated = namespace['build_background_strata'](matrix)
    assert strata.sort_index().equals(recreated.sort_index())
    cfg = json.loads((p/'analysis_config.json').read_text(encoding='utf-8-sig'))
    benchmarks = set(roles.loc[roles.record_role.eq('V'),'gene_symbol'])
    for axis in ('translational','functional'):
        pos = set(roles.loc[roles.record_role.eq('L') & roles.axis.eq(axis),'gene_symbol'])
        ordered_pos = [g for g in genes if g in pos]
        reserved = pos | benchmarks
        unlabeled = [g for g in genes if g not in reserved]
        alias_path = HERE/'identifier_results'/f'{cancer}_{axis}_reserved_gene_aliases_in_eligible_background.tsv'
        aliases = pd.read_csv(alias_path, sep='\t', keep_default_na=False)
        alias_set = set(aliases.gene_symbol)
        n_seeds = 0
        n_selected = 0
        for seed in cfg['random_seeds'][:30]:
            background = sample(ordered_pos, unlabeled, strata, cfg['background_sampling']['background_ratio'], np.random.default_rng(seed))
            hits = sorted(set(background) & alias_set)
            n_seeds += bool(hits)
            n_selected += len(hits)
            selection_rows.append(dict(cancer=cancer,axis=axis,seed=seed,n_background=len(background),
                                       n_reserved_aliases_drawn=len(hits)))
            for g in hits:
                selected_aliases.append(dict(cancer=cancer,axis=axis,seed=seed,gene_symbol=g,
                                             approved_symbol=aliases.set_index('gene_symbol').at[g,'approved_symbol']))
        summaries.append(dict(cancer=cancer,axis=axis,n_seeds=30,
                              n_eligible_reserved_aliases=len(alias_set),n_seeds_with_reserved_aliases=n_seeds,
                              n_reserved_alias_draws=n_selected, frozen_strata_reconstructed=True))
pd.DataFrame(selection_rows).to_csv(HERE/'identifier_results/primary_background_sampling_audit.tsv',sep='\t',index=False)
pd.DataFrame(selected_aliases,columns=['cancer','axis','seed','gene_symbol','approved_symbol']).to_csv(HERE/'identifier_results/primary_background_selected_reserved_aliases.tsv',sep='\t',index=False)
result=dict(scope='Reconstruct the 120 primary sampling draws from frozen row order, strata, labels and seeds; no model fitting and no claim to recover historical execution logs.',
            summaries=summaries,selected_aliases=selected_aliases)
(HERE/'identifier_results/background_sampling_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
