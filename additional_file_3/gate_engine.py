"""Parameterized operational review; frozen score inputs, no efficacy labels."""
import numpy as np
import pandas as pd

def base_policy(cfg):
    r=cfg['review']
    return dict(high_score_fraction=float(r['high_score_percentile']),safety_min=float(r['safety_factor_min']),
        dependency_top_fraction=float(r['dependency_top_percentile']),disease_top_fraction=float(r['disease_top_percentile']),
        dependency_reference_multiplier=.5,disease_reference_multiplier=.5,safety_risk_high=float(r['safety_risk_high_threshold']),
        clinical_druggability_min=.1,limited_druggability_max=float(r.get('limited_actionability_druggability_max',.15)),
        clinical_override_min=float(r.get('clinical_precedence_override_min',.6)),opentargets_override_min=float(r.get('opentargets_override_min',.5)),
        override_enabled=bool(r.get('allow_strong_translational_evidence_without_dependency',False)),
        limited_before_functional=bool(r.get('limited_actionability_before_functional_class',False)))

GRID=dict(high_score_fraction=[.01,.02,.05],safety_min=[.3,.4,.5],
    dependency_top_fraction=[.05,.1,.2],disease_top_fraction=[.05,.1,.2],
    dependency_reference_multiplier=[.25,.5,1.],disease_reference_multiplier=[.25,.5,1.],
    safety_risk_high=[.65,.75,.85],clinical_druggability_min=[.05,.1,.2],
    limited_druggability_max=[.1,.15,.2],clinical_override_min=[.4,.6,.8],opentargets_override_min=[.3,.5,.7])

def policies(cfg):
    base=base_policy(cfg);out=[];seen=set()
    def add(name,kind,changes):
        p=dict(base,**changes);key=tuple(p.items())
        if key in seen:return
        seen.add(key);out.append(dict(policy_id=name,policy_type=kind,changed_parameters=','.join(changes),parameters=p))
    add('reference','reference',{})
    for field,values in GRID.items():
        assert base[field] in values,(field,base[field])
        for val in values:
            if val!=base[field]:add(field+'_'+str(val),'one_factor',{field:val})
    for field in ('override_enabled','limited_before_functional'):
        add(field+'_toggle','rule_stress',{field:not base[field]})
    for high in (.01,.02,.05):
        for safety in (.3,.4,.5):add('legacy_grid_'+str(high)+'_'+str(safety),'legacy_grid',dict(high_score_fraction=high,safety_min=safety))
    add('joint_more_selective','joint_stress',dict(high_score_fraction=.01,safety_min=.5,
        dependency_top_fraction=.05,disease_top_fraction=.05,dependency_reference_multiplier=1.,disease_reference_multiplier=1.,
        safety_risk_high=.65,clinical_druggability_min=.2,limited_druggability_max=.1,clinical_override_min=.8,opentargets_override_min=.7,override_enabled=False))
    add('joint_less_selective','joint_stress',dict(high_score_fraction=.05,safety_min=.3,
        dependency_top_fraction=.2,disease_top_fraction=.2,dependency_reference_multiplier=.25,disease_reference_multiplier=.25,
        safety_risk_high=.85,clinical_druggability_min=.05,limited_druggability_max=.2,clinical_override_min=.4,opentargets_override_min=.3))
    assert len(out)==31,len(out)
    return out

def cutoffs(m,scores,p):
    disease=m['disease_relevance_score'] if 'disease_relevance_score' in m else m.patient_alteration_score
    return dict(trans=float(scores.translational_precedence_score.quantile(1-p['high_score_fraction'])),
        func=float(scores.functional_discovery_score.quantile(1-p['high_score_fraction'])),
        dep=float(m.dependency_score.quantile(1-p['dependency_top_fraction'])),
        disease=float(disease.quantile(1-p['disease_top_fraction'])))

def classify(m,scores,p,risk=None):
    c=cutoffs(m,scores,p);get=lambda n:m[n].to_numpy(dtype=float)
    dep=get('dependency_score');drug=get('druggability_score')
    disease=get('disease_relevance_score' if 'disease_relevance_score' in m else 'patient_alteration_score')
    old=get('safety_risk_score');risk=old if risk is None else np.asarray(risk,dtype=float)
    safety=get('safety_factor') if np.array_equal(risk,old,equal_nan=True) else 1-risk
    strong=(get('clinical_precedence_score')>=p['clinical_override_min'])|(get('opentargets_overall_score')>=p['opentargets_override_min'])
    clinical=(scores.translational_precedence_score.to_numpy()>=c['trans'])&(drug>=p['clinical_druggability_min'])&(safety>=p['safety_min'])&((dep>=c['dep']*p['dependency_reference_multiplier'])|(p['override_enabled']&strong))
    limited=(dep>=c['dep'])&(drug<p['limited_druggability_max'])
    functional=(scores.functional_discovery_score.to_numpy()>=c['func'])&(dep>=c['dep']*p['dependency_reference_multiplier'])&(disease>=c['disease']*p['disease_reference_multiplier'])&(safety>=p['safety_min'])
    cond=[risk>=p['safety_risk_high'],clinical];value=[4,1]
    if p['limited_before_functional']:cond.append(limited);value.append(3)
    cond += [functional,limited,(disease>=c['disease'])&(dep<c['dep']*p['dependency_reference_multiplier'])];value += [2,3,5]
    return np.select(cond,value,default=6).astype(np.uint8)

def scalar_reference(row,score,p,c,risk=None):
    """Explicit ordered decision reference, independently from vector selection."""
    r=float(row.safety_risk_score) if risk is None else float(risk)
    sf=float(row.safety_factor) if risk is None else 1-r
    if r>=p['safety_risk_high']:return 4
    strong=(row.clinical_precedence_score>=p['clinical_override_min'] or row.opentargets_overall_score>=p['opentargets_override_min'])
    if score.translational_precedence_score>=c['trans'] and row.druggability_score>=p['clinical_druggability_min'] and sf>=p['safety_min'] and (row.dependency_score>=c['dep']*p['dependency_reference_multiplier'] or (p['override_enabled'] and strong)):return 1
    limited=row.dependency_score>=c['dep'] and row.druggability_score<p['limited_druggability_max']
    if p['limited_before_functional'] and limited:return 3
    disease=getattr(row,'disease_relevance_score',row.patient_alteration_score)
    if score.functional_discovery_score>=c['func'] and row.dependency_score>=c['dep']*p['dependency_reference_multiplier'] and disease>=c['disease']*p['disease_reference_multiplier'] and sf>=p['safety_min']:return 2
    if limited:return 3
    if disease>=c['disease'] and row.dependency_score<c['dep']*p['dependency_reference_multiplier']:return 5
    return 6

def completion_bits(m,scores,p,lo,hi):
    bits=np.zeros(len(m),dtype=np.uint8)
    points=[lo,hi]
    for b in (1-p['safety_min'],p['safety_risk_high']):
        for step in (-1e-12,0.,1e-12):points.append(np.clip(b+step,lo,hi))
    for risk in points:
        cl=classify(m,scores,p,risk=risk);bits|=np.left_shift(np.uint8(1),cl-1)
    return bits
