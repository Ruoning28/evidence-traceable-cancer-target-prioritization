"""Canonical gnomAD observation states and within-namespace constraint intervals."""
import numpy as np
import pandas as pd
FIELDS=['gene','gene_id','transcript','canonical','gene_flags','lof.oe_ci.upper','lof.pLI']

def normalize(s):return s.astype(str).str.strip().str.replace(r'\s+\(\d+\)$','',regex=True).str.split('|',regex=False).str[0].str.upper()

def read_constraints(raw_path):
    frames=[];n=0
    for chunk in pd.read_csv(raw_path,sep='\t',compression='gzip',usecols=lambda c:c in FIELDS,chunksize=50000,low_memory=False):
        n+=len(chunk);chunk=chunk[chunk.canonical.astype(str).str.lower().isin(['true','1'])].copy()
        chunk['symbol']=normalize(chunk.gene);frames.append(chunk)
    return pd.concat(frames,ignore_index=True),n

def constraint_bounds(g):
    g=g.copy();lo=pd.to_numeric(g['lof.oe_ci.upper'],errors='coerce');pli=pd.to_numeric(g['lof.pLI'],errors='coerce')
    gl=lo.notna()&np.isfinite(lo)&lo.ge(0);gp=pli.notna()&np.isfinite(pli)&pli.between(0,1)
    lower,upper=float(lo[gl].min()),float(lo[gl].max());assert upper>lower
    ls=pd.Series(np.nan,index=g.index);ls[gl]=(1-(lo[gl]-lower)/(upper-lower)).clip(0,1)
    ps=pli.where(gp)
    g['lower']=np.maximum(ls.fillna(0),ps.fillna(0));g['upper']=np.maximum(ls.fillna(1),ps.fillna(1))
    g['loeuf_observed']=gl;g['pli_observed']=gp;g['both_observed']=gl&gp
    g['quality_flagged']=~g.gene_flags.fillna('').astype(str).str.replace(' ','',regex=False).isin(['','[]','NA','nan','None'])
    ids=g.gene_id.fillna('').astype(str);g['validid']=ids.where(~ids.isin(['','NA','nan']),None)
    out=g.groupby('symbol').agg(lower=('lower','min'),upper=('upper','max'),loeuf_observed=('loeuf_observed','any'),
        pli_observed=('pli_observed','any'),both_observed=('both_observed','all'),canonical_rows=('symbol','size'),
        quality_flagged=('quality_flagged','any'),distinct_geneids=('validid','nunique'))
    g['id_namespace']=np.select([ids.str.match(r'^ENSG[0-9]+(?:\.[0-9]+)?$'),ids.str.match(r'^[0-9]+$')],['ENSEMBL','NCBI'],default='OTHER')
    namespace_counts=g.groupby(['symbol','id_namespace']).validid.nunique()
    ambiguous=namespace_counts.gt(1).groupby(level='symbol').any().reindex(out.index).fillna(False)
    out.loc[ambiguous,['lower','upper']]=[0.,1.]
    out['ambiguous']=ambiguous
    return out,dict(min_LOEUF=lower,max_LOEUF=upper,scaling='Original finitecanonicalmin/max, retained only forformula audit; not biologicallycalibrated risk.')
