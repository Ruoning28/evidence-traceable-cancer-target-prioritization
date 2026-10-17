from pathlib import Path
from collections import defaultdict
import ast
import hashlib
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
CODE=args.code.resolve()/'scientific_source'
OUT = HERE / 'source_results'
OUT.mkdir(exist_ok=True)
revision_path = CODE/'core/scripts/revision_analyses.py'
tree = ast.parse(revision_path.read_text(encoding='utf-8-sig'))
families = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
                and any(isinstance(t,ast.Name) and t.id=='FAMILIES' for t in n.targets))
fn = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='descendant_exclusions')
ns = {'FAMILIES': families}
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(revision_path),'exec'),ns)

def descendants(roots, graph):
    seen=set(roots)
    while True:
        new={child for node in seen for child in graph.get(node,[])}-seen
        if not new:
            return seen
        seen |= new

def graph_for(cancer):
    graph=defaultdict(set)
    edges=[]
    def add(parent,children,evidence):
        for child in children:
            graph[parent].add(child)
            edges.append(dict(cancer=cancer,parent=parent,child=child,evidence=evidence))
    add('chemical_modulation_score',['druggability_score'],'reference/v4_end_to_end.py:682-686')
    add('clinical_precedence_score',['druggability_score'],'reference/v4_end_to_end.py:682-686')
    add('target_tractability_score',['druggability_score'],'reference/v4_end_to_end.py:682-686')
    add('_shared_OpenTargets_annotation',['chemical_modulation_score','clinical_precedence_score','opentargets_overall_score'],
        'reference/v4_end_to_end.py:561-572,666-680; shared drug annotations retained in LUAD UNIVERSAL_COLUMNS')
    add('_DepMap_resource',['dependency_score','common_essential_risk_score'],'reference/v4_end_to_end.py:build_dependency_features/build_safety_features')
    add('common_essential_risk_score',['safety_risk_score'],'reference/v4_end_to_end.py:721-728')
    add('normal_expression_risk_score',['safety_risk_score'],'reference/v4_end_to_end.py:695-705,726-728')
    add('lof_constraint_risk_score',['safety_risk_score'],'reference/v4_end_to_end.py:706-728')
    add('safety_risk_score',['safety_factor'],'core/src/framework/pipeline.py:build_feature_matrix')
    if cancer=='PDAC':
        add('opentargets_overall_score',['chemical_modulation_score'],
            'reference/v4_end_to_end.py:561-572,666-668; the overall field is one counted channel')
        add('opentargets_clinical_score',['chemical_modulation_score','clinical_precedence_score'],
            'reference/v4_end_to_end.py:666-680')
        for c in ['cptac_mrna_score','cnv_score','mutation_score','tcga_patient_score']:
            add(c,['patient_alteration_score'],'reference/v4_end_to_end.py:497-499,1675-1677')
        tcga=['tcga_gtex_expression_score','tcga_mutation_score','tcga_cnv_score','tcga_survival_score']
        add('_TCGA_resource',tcga,'reference/v4_end_to_end.py:build_tcga_features')
        for c in tcga:
            add(c,['tcga_patient_score'],'reference/v4_end_to_end.py:424-425')
        for c in ['cptac_protein_score','phosphoproteomic_support_score']:
            add(c,['proteomic_support_score'],'reference/v4_end_to_end.py:500-501')
        add('geo_external_reproducibility_score',['external_reproducibility_score'],'reference/v4_end_to_end.py:1670-1673')
        add('external_reproducibility_score',['geo_external_reproducibility_score'],'Verified numeric identity / equivalent information, reverse edge for deletion audit only')
    else:
        add('_TCGA_resource',['tcga_gtex_expression_score','tcga_cnv_score','tcga_survival_score','_TCGA_RPPA'],
            'reference/build_luad_inputs.py:239-291; mutation is zero at patient-composite construction')
        add('tcga_gtex_expression_score',['patient_alteration_score','cptac_mrna_score','external_reproducibility_score'],
            'reference/build_luad_inputs.py:283-288,386-389; cptac_mrna path is fallback')
        add('tcga_cnv_score',['patient_alteration_score'],'reference/build_luad_inputs.py:283-287')
        add('_TCGA_RPPA',['cptac_protein_score','proteomic_support_score'],
            'reference/build_luad_inputs.py:289-290,374-378; per-field fallback')
        add('dependency_score',['external_reproducibility_score'],'reference/build_luad_inputs.py:386-389')
        for c in ['cptac_protein_score','phosphoproteomic_support_score','cptac_acetylprotein_score']:
            add(c,['proteomic_support_score'],'reference/build_luad_inputs.py:207-209; CPTAC-preferred branch')
    for c in ['cptac_protein_score','phosphoproteomic_support_score','geo_external_reproducibility_score','opentargets_overall_score']:
        add(c,[c+'_missing'],'core/src/framework/pipeline.py:291-320; missingness may carry source availability')
    return graph,edges

fixed_rows, source_rows, edge_rows, checks = [],[],[],[]
archived = pd.read_csv(DATA/'historical_tables/strict_tolerance_revision_fixed_benchmark_source_ablation.tsv',sep='\t')
lso = pd.read_csv(DATA/'historical_tables/strict_tolerance_framework_leave_source_out_validation.tsv',sep='\t',keep_default_na=False)
for cancer in ['PDAC','LUAD']:
    cfg=json.loads((DATA/'frozen_inputs'/cancer/'analysis_config.json').read_text(encoding='utf-8-sig'))
    f=pd.read_csv(DATA/'frozen_inputs'/cancer/'framework_feature_matrix.tsv',sep='\t')
    m=pd.read_csv(DATA/'frozen_inputs'/cancer/'source_target_evidence_matrix.tsv',sep='\t',low_memory=False).set_index('gene_symbol')
    graph,edges=graph_for(cancer);edge_rows+=edges
    for axis,spec in cfg['evaluation']['axes'].items():
        excluded=set(spec['feature_exclusions'])
        excluded|={c+'_missing' for c in list(excluded)}
        allowed=set(f.columns)-{'gene_symbol'}-excluded
        for family, roots in families.items():
            actual=allowed & ns['descendant_exclusions'](cancer,family)
            # Trace descendants of the declared feature block. This is not a
            # claim that those blocks partition raw databases into disjoint sets.
            complete=allowed & descendants(roots,graph)
            retained=complete-actual
            same_resource=set()
            if family=='external_knowledge':
                # In LUAD the exposed disease-specific OT field was replaced,
                # while shared chemical/clinical annotations retain the prior
                # PDAC/general OT input. They share a resource, not a value.
                same_resource=(allowed & descendants({'_shared_OpenTargets_annotation'},graph))-actual
            row=archived[(archived.cancer==cancer)&(archived.axis==axis)&(archived.removed_family==family)].iloc[0]
            stored=set(str(row.excluded_model_columns).split(';')) if pd.notna(row.excluded_model_columns) else set()
            assert stored==actual,(cancer,axis,family,stored,actual)
            fixed_rows.append(dict(cancer=cancer,axis=axis,family=family,
                                   actual_deleted_columns=';'.join(sorted(actual)),
                                   retained_descendant_columns=';'.join(sorted(retained)),
                                   extra_deletions_for_declared_block_closure=';'.join(sorted(retained)),
                                   retained_same_resource_path_columns=';'.join(sorted(same_resource)),
                                   extra_deletions_for_audited_source_closure=';'.join(sorted(retained|same_resource)),
                                   complete_for_declared_feature_dependencies=not retained,
                                   complete_for_audited_resource_paths=not (retained|same_resource),
                                   archived_top50_hits=int(row.top50_hits),benchmark_n=int(row.n_benchmark_targets),
                                   stored_column_list_matches_code=True))
    for r in lso[lso.cancer.eq(cancer)].itertuples():
        spec=cfg['evaluation']['axes'][r.axis]
        base=set(spec['feature_exclusions']);base|={c+'_missing' for c in list(base)}
        allowed=set(f.columns)-{'gene_symbol'}-base
        explicit=set(cfg['evaluation']['source_feature_exclusions'].get(r.heldout_source,[]))
        actual=allowed & (explicit|{c+'_missing' for c in explicit})
        roots=set(explicit)
        if r.heldout_source in ['known_dependency','depmap_luad_high_dependency']:
            roots.add('_DepMap_resource')
        if r.heldout_source in ['known_pdac_driver','tcga_luad_nature_2014']:
            roots.add('_TCGA_resource')
        if r.heldout_source=='curated_literature:OpenTargets2025':
            roots.add('_shared_OpenTargets_annotation')
        complete=allowed & descendants(roots,graph)
        remaining=complete-actual
        active=[c for c in remaining if c in f and f[c].nunique(dropna=False)>1]
        source_rows.append(dict(cancer=cancer,axis=r.axis,heldout_source=r.heldout_source,status=r.status,
                                configured_exclusions=';'.join(sorted(explicit)),
                                actual_additional_deleted_model_columns=';'.join(sorted(actual)),
                                retained_potential_path_columns=';'.join(sorted(remaining)),
                                nonconstant_retained_path_columns=';'.join(sorted(active)),
                                background_matching='Full frozen strata retained, including database coverage, dependency availability and expression coverage.',
                                diagnosis='Resource/derivation-path audit; not proof of identical label-defining facts in these columns.'))
    # Numeric verification of the main shared composites, without fitting.
    drug=m[['chemical_modulation_score','clinical_precedence_score','target_tractability_score']].fillna(0).mean(axis=1)
    risk=m[['normal_expression_risk_score','lof_constraint_risk_score','common_essential_risk_score']].fillna(0).mean(axis=1)
    checks.append(dict(cancer=cancer,formula='druggability_mean',max_abs_error=float((drug-m.druggability_score).abs().max())))
    checks.append(dict(cancer=cancer,formula='risk_mean',max_abs_error=float((risk-m.safety_risk_score).abs().max())))
    if cancer=='PDAC':
        cols=['cptac_mrna_score','cnv_score','mutation_score','tcga_patient_score']
        value=1-(1-m[cols].fillna(0).clip(0,1)).prod(axis=1)
        checks.append(dict(cancer=cancer,formula='patient_noisy_OR_CPTAC_TCGA',max_abs_error=float((value-m.patient_alteration_score).abs().max())))
        checks.append(dict(cancer=cancer,formula='GEO_external_numeric_alias',max_abs_error=float((m.external_reproducibility_score-m.geo_external_reproducibility_score.fillna(0)).abs().max())))
    else:
        value=.5*m.tcga_gtex_expression_score.fillna(0)+.5*m.dependency_score.fillna(0)
        checks.append(dict(cancer=cancer,formula='TCGA_DepMap_external_composite',max_abs_error=float((value-m.external_reproducibility_score).abs().max())))
        # Zero mutation component at construction: later CPTAC mutation update
        # must not be added retroactively to this patient-composite dependency.
        value=.45*m.tcga_cnv_score.fillna(0)+.45*m.tcga_gtex_expression_score.fillna(0)
        checks.append(dict(cancer=cancer,formula='patient_without_later_mutation_update',max_abs_error=float((value-m.patient_alteration_score).abs().max())))
for c in checks:
    assert c['max_abs_error']<1e-12,c
pd.DataFrame(edge_rows).to_csv(OUT/'verified_derivation_edges.tsv',sep='\t',index=False)
pd.DataFrame(fixed_rows).to_csv(OUT/'fixed_benchmark_ablation_scope.tsv',sep='\t',index=False)
pd.DataFrame(source_rows).to_csv(OUT/'source_role_out_scope.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(OUT/'composite_numeric_checks.tsv',sep='\t',index=False)
result=dict(fixed_benchmark_rows=len(fixed_rows),fixed_rows_with_residual_descendants=[r for r in fixed_rows if not r['complete_for_audited_resource_paths']],
            source_role_out=source_rows,numeric_checks=checks,
            code_sha256={str(p.relative_to(CODE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [revision_path,CODE/'core/src/framework/pipeline.py',CODE/'reference/v4_end_to_end.py',CODE/'reference/build_luad_inputs.py']},
            boundaries='Paths verified against delivered code and selected frozen numeric identities. Source-wide deletion is distinct from deletion of a named feature block; no new fitting, transfer evaluation, or causal attribution is performed.')
(OUT/'source_audit_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(fixed_rows)[['cancer','axis','family','retained_descendant_columns','retained_same_resource_path_columns']].to_string(index=False))
print(pd.DataFrame(source_rows)[['cancer','axis','heldout_source','status','actual_additional_deleted_model_columns','retained_potential_path_columns']].to_string(index=False))
print('Numeric checks:',json.dumps(checks))
