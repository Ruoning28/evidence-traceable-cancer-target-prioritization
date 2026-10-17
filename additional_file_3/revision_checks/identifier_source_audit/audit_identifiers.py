from pathlib import Path
from collections import defaultdict
import hashlib
import json
import re
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
HGNC=args.reference.resolve()/'hgnc_complete_set.txt'
OUT = HERE / 'identifier_results'
OUT.mkdir(exist_ok=True)

h = pd.read_csv(HGNC, sep='\t', dtype=str, keep_default_na=False)
assert h.hgnc_id.is_unique
approved = {r.symbol.upper(): r.hgnc_id for r in h.itertuples()}
info = h.set_index('hgnc_id')
aliases = defaultdict(set)
identifiers = defaultdict(set)
uniprot_to_hgnc = defaultdict(set)
for r in h.itertuples():
    for field in ('alias_symbol', 'prev_symbol'):
        for token in getattr(r, field).split('|'):
            if token:
                aliases[token.upper()].add(r.hgnc_id)
    for field in ('hgnc_id', 'ensembl_gene_id', 'entrez_id'):
        for token in getattr(r, field).split('|'):
            if token:
                identifiers[token.upper()].add(r.hgnc_id)
    for token in r.uniprot_ids.split('|'):
        if token:
            identifiers[token.upper()].add(r.hgnc_id)
            uniprot_to_hgnc[token.upper()].add(r.hgnc_id)

uniprot_path = args.reference.resolve()/'uniprot_gene_names.tsv.gz'
u = pd.read_csv(uniprot_path, sep='\t', dtype=str, keep_default_na=False, usecols=['Entry', 'Gene Names'])
frozen_names = defaultdict(set)
for r in u.itertuples(index=False, name=None):
    entry, names = r
    if names and len(uniprot_to_hgnc[entry.upper()]) == 1:
        frozen_names[names.upper().strip()].update(uniprot_to_hgnc[entry.upper()])

def simple(key):
    if key in approved:
        return {approved[key]}, 'approved_symbol'
    key2 = re.sub(r'^(ENSG\d+)\.\d+$', r'\1', key)
    if key2 in identifiers:
        return identifiers[key2], 'stable_identifier'
    if key in aliases:
        return aliases[key], 'previous_or_alias_symbol'
    return set(), 'unmapped'

def resolve(value):
    key = str(value).strip().upper()
    targets, how = simple(key)
    if targets:
        if len(targets) == 1:
            return next(iter(targets)), how, ';'.join(sorted(targets))
        return '', 'ambiguous_alias_or_identifier', ';'.join(sorted(targets))
    if key in frozen_names:
        targets = frozen_names[key]
        if len(targets) == 1:
            return next(iter(targets)), 'frozen_uniprot_gene_names', ';'.join(sorted(targets))
        return '', 'ambiguous_frozen_uniprot', ';'.join(sorted(targets))
    if re.search(r'[\s;,/|]', key):
        tokens = [t for t in re.split(r'[\s;,/|]+', key) if t]
        resolved = [simple(t)[0] for t in tokens]
        if resolved and all(len(s) == 1 for s in resolved) and len(set.union(*resolved)) == 1:
            targets = set.union(*resolved)
            return next(iter(targets)), 'same_gene_multi_token', ';'.join(sorted(targets))
        return '', 'unresolved_multi_token', ';'.join(sorted(set().union(*resolved)))
    return '', 'unmapped', ''

all_summary, metrics, top_rows, collision_rows, overlap_rows = [], [], [], [], []
checks = []
for cancer in ('PDAC', 'LUAD'):
    srcdir = DATA / 'frozen_inputs' / cancer
    source = pd.read_csv(srcdir / 'source_target_evidence_matrix.tsv', sep='\t', usecols=['gene_symbol'], keep_default_na=False)
    genes = source.gene_symbol.astype(str)
    rows = []
    for g in genes:
        hg, how, candidates = resolve(g)
        rows.append(dict(gene_symbol=g, hgnc_id=hg, mapping_status=how, candidate_hgnc_ids=candidates,
                         approved_symbol=info.at[hg, 'symbol'] if hg else '',
                         locus_group=info.at[hg, 'locus_group'] if hg else ''))
    mapping = pd.DataFrame(rows)
    mapping.to_csv(OUT / f'{cancer}_record_mapping.tsv', sep='\t', index=False)
    resolved = mapping[mapping.hgnc_id.ne('')]
    collision = resolved[resolved.hgnc_id.duplicated(keep=False)]
    collision.to_csv(OUT / f'{cancer}_duplicate_gene_records.tsv', sep='\t', index=False)
    summary = dict(cancer=cancer, candidate_records=len(genes), unique_record_keys=genes.nunique(),
                   uniquely_mapped_records=len(resolved), uniquely_mapped_human_genes=resolved.hgnc_id.nunique(),
                   duplicate_gene_groups=collision.hgnc_id.nunique(), duplicate_gene_records=len(collision),
                   redundant_mapped_records=len(resolved)-resolved.hgnc_id.nunique(),
                   mapping_status_counts=mapping.mapping_status.value_counts().to_dict(),
                   locus_group_counts_unique_genes=resolved.drop_duplicates('hgnc_id').locus_group.value_counts().to_dict())
    all_summary.append(summary)
    m = mapping.set_index('gene_symbol')
    scores = pd.read_csv(DATA / 'review_inputs' / cancer / 'archived_strict_priority_scores.tsv', sep='\t')
    assert scores.gene_symbol.is_unique
    roles = pd.read_csv(srcdir / 'framework_record_role_manifest.tsv', sep='\t', keep_default_na=False)
    benchmarks_all = pd.read_csv(srcdir / 'framework_validation_benchmarks.tsv', sep='\t', keep_default_na=False)
    benchmarks_union = set(benchmarks_all.gene_symbol)
    for axis, prefix in [('translational', 'translational_precedence'), ('functional', 'functional_discovery')]:
        scorecol = prefix + '_score'
        ordered = scores.sort_values(scorecol, ascending=False).copy()
        ordered['original_position'] = np.arange(1, len(ordered)+1)
        for col in ('hgnc_id', 'approved_symbol', 'mapping_status'):
            ordered[col] = ordered.gene_symbol.map(m[col])
        ordered['identity_key'] = np.where(ordered.hgnc_id.ne(''), ordered.hgnc_id, 'UNRESOLVED:' + ordered.gene_symbol)
        collapsed = ordered.drop_duplicates('identity_key', keep='first')
        mapped_only = ordered[ordered.hgnc_id.ne('')].drop_duplicates('hgnc_id', keep='first')
        benchmark = set(benchmarks_all.loc[benchmarks_all.axis.eq(axis), 'gene_symbol'])
        benchmark_hg = {resolve(g)[0] for g in benchmark if resolve(g)[0]}
        assert len(benchmark_hg) == len(benchmark), 'Benchmark mapping must be complete and unique for denominator comparison'
        positives = set(roles.loc[roles.axis.eq(axis) & roles.record_role.eq('L'), 'gene_symbol'])
        pos_hg = {resolve(g)[0] for g in positives if resolve(g)[0]}
        background_excluded = positives | benchmarks_union
        reserved_hg = {resolve(g)[0] for g in background_excluded if resolve(g)[0]}
        alias_bg = resolved[~resolved.gene_symbol.isin(background_excluded) & resolved.hgnc_id.isin(reserved_hg)]
        alias_bg.assign(cancer=cancer, axis=axis).to_csv(OUT / f'{cancer}_{axis}_reserved_gene_aliases_in_eligible_background.tsv', sep='\t', index=False)
        overlap_rows.append(dict(cancer=cancer, axis=axis, positives=len(positives), benchmark=len(benchmark),
                                 unmapped_positives=sum(not resolve(g)[0] for g in positives),
                                 normalized_L_V_overlap=len(pos_hg & benchmark_hg),
                                 reserved_gene_aliases_in_eligible_background=len(alias_bg),
                                 reserved_genes_with_alias_background=alias_bg.hgnc_id.nunique()))
        for k in (20, 50, 100, 200):
            original_top = ordered.head(k)
            mapped_counts = original_top.loc[original_top.hgnc_id.ne(''), 'hgnc_id'].value_counts()
            for policy, table in [('original_record_ranking', ordered),
                                  ('collapse_unique_mappings_keep_unresolved', collapsed),
                                  ('mapped_human_genes_only', mapped_only)]:
                top = table.head(k)
                exact_hits = len(set(top.gene_symbol) & benchmark)
                normalized_hits = len(set(top.hgnc_id[top.hgnc_id.ne('')]) & benchmark_hg)
                original_exact = len(set(original_top.gene_symbol) & benchmark)
                metrics.append(dict(cancer=cancer, axis=axis, k=k, policy=policy, universe_size=len(table),
                                    benchmark_denominator=len(benchmark), exact_key_hits=exact_hits,
                                    normalized_gene_hits=normalized_hits, original_exact_hits=original_exact,
                                    unique_mapped_genes_in_original_topk=original_top.loc[original_top.hgnc_id.ne(''), 'hgnc_id'].nunique(),
                                    unresolved_records_in_original_topk=int(original_top.hgnc_id.eq('').sum()),
                                    duplicate_gene_groups_in_original_topk=int(mapped_counts.gt(1).sum()),
                                    redundant_mapped_records_in_original_topk=int((mapped_counts-1).clip(lower=0).sum()),
                                    original_record_overlap=len(set(top.gene_symbol) & set(original_top.gene_symbol))/k))
                if k == 50:
                    out = top.assign(cancer=cancer, axis=axis, policy=policy, position=np.arange(1,len(top)+1),
                                     exact_benchmark_hit=top.gene_symbol.isin(benchmark),
                                     normalized_benchmark_hit=top.hgnc_id.isin(benchmark_hg))
                    top_rows.extend(out[['cancer','axis','policy','position','gene_symbol','hgnc_id','approved_symbol',
                                         'mapping_status','original_position','exact_benchmark_hit','normalized_benchmark_hit']].to_dict('records'))
            if k == 50:
                dup_hg = set(mapped_counts[mapped_counts.gt(1)].index)
                collision_rows.extend(original_top[original_top.hgnc_id.isin(dup_hg)].assign(cancer=cancer,axis=axis)[
                    ['cancer','axis','gene_symbol','hgnc_id','approved_symbol','mapping_status','original_position']].to_dict('records'))
        checks.append(dict(cancer=cancer, axis=axis,
                           original_top50_hits=len(set(ordered.head(50).gene_symbol)&benchmark),
                           input_records_equal_scores=set(genes)==set(scores.gene_symbol)))

pd.DataFrame(metrics).to_csv(OUT / 'ranking_identity_sensitivity.tsv', sep='\t', index=False)
pd.DataFrame(top_rows).to_csv(OUT / 'top50_identity_details.tsv', sep='\t', index=False)
pd.DataFrame(collision_rows, columns=['cancer','axis','gene_symbol','hgnc_id','approved_symbol','mapping_status','original_position']).to_csv(OUT / 'top50_duplicate_gene_groups.tsv', sep='\t', index=False)
pd.DataFrame(overlap_rows).to_csv(OUT / 'label_and_background_identity_audit.tsv', sep='\t', index=False)
result = dict(reference_rows=len(h), hgnc_sha256=hashlib.sha256(HGNC.read_bytes()).hexdigest(),
              uniprot_sha256=hashlib.sha256(uniprot_path.read_bytes()).hexdigest(),
              mapping_rules='Exact approved symbol takes precedence; then stable identifier; then unique union of previous/alias symbols. Frozen UniProt full Gene Names fields may map via unique HGNC accession. Multi-token strings map only if every token maps uniquely to the same gene. Ambiguous and unmapped records are retained unresolved.',
              interpretation='Frozen-ranking post hoc identity sensitivity, not a refitted model or new independent benchmark.',
              summaries=all_summary, label_audit=overlap_rows, checks=checks)
(OUT / 'identifier_audit_summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=True, indent=2))
print(pd.DataFrame(metrics).query('k == 50').to_string(index=False))
