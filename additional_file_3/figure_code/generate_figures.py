"""Generate Figures 3--6 from frozen scores and supplied figure tables."""
from pathlib import Path
import argparse, hashlib, importlib.util, json, os, sys
if os.name=='nt':
    lib=Path(sys.prefix)/'Library/bin'
    if lib.exists():
        DLL_HANDLE=os.add_dll_directory(str(lib))
        os.environ['PATH']=str(lib)+os.pathsep+os.environ.get('PATH','')
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
SOURCES=HERE/'figure_source'
spec=importlib.util.spec_from_file_location('figure_style',HERE/'plotting_style.py')
style=importlib.util.module_from_spec(spec);spec.loader.exec_module(style)
plt=style.plt;Rectangle=style.Rectangle
BLUE=style.AXIS_COLORS['translational'];ORANGE=style.AXIS_COLORS['functional'];GRAY=style.AXIS_COLORS['baseline']
PROTECTED={'IL1A','SLC5A2','ACOT8','SETD1A','PRDX1','ERO1A','UQCRC1','GOT1'}
EXAMPLES=['KRAS','ARID1A','KIF23','PALB2','ERBB2']
DATA=None
OUTPUT=None

def read(c,n):
    if n=='framework_priority_scores.tsv':
        return pd.read_csv(DATA/'review_inputs'/c/'archived_strict_priority_scores.tsv',sep='\t',float_precision='round_trip')
    num='3' if c=='LUAD' else '4'
    if n=='framework_review_classes.tsv':return pd.read_csv(SOURCES/('fig'+num+'_classes.tsv'),sep='\t')
    if n=='framework_evidence_cards.tsv':return pd.read_csv(SOURCES/('fig'+num+'_original_case_inputs.tsv'),sep='\t')
    cache={'revision_permutation_summary.tsv':'fig5_all_permutation_summaries.tsv',
           'framework_leave_source_out_validation.tsv':'fig5_source_role_out.tsv',
           'revision_stability_summary.tsv':'fig5_top50_stability.tsv',
           'framework_oof_validation_summary.tsv':c+'_fig5_oof.tsv',
           'framework_independent_benchmark_validation.tsv':c+'_benchmark_aggregate.tsv'}
    path=SOURCES/cache[n]
    if not path.exists():path=DATA/'historical_tables'/c/n
    frame=pd.read_csv(path,sep='\t')
    if 'cancer' in frame:frame=frame[frame.cancer.eq(c)].copy()
    return frame

def safe(c):
    scores=read(c,'framework_priority_scores.tsv');scores=scores[~scores.gene_symbol.isin(PROTECTED)].copy()
    cards=read(c,'framework_evidence_cards.tsv');cards=cards[cards.gene_symbol.isin(EXAMPLES)].copy()
    return scores,read(c,'framework_review_classes.tsv'),cards

def source(frame,name):
    folder=OUTPUT/'figure_source';folder.mkdir(exist_ok=True)
    frame.to_csv(folder/(name+'.tsv'),sep='\t',index=False,encoding='utf-8-sig')

def save(fig,name):
    for ext in ['pdf','svg','png']:
        fig.savefig(OUTPUT/(name+'.'+ext),bbox_inches='tight',**({'dpi':350} if ext=='png' else {}))
    plt.close(fig)

def panel(ax,letter,title):style.panel(ax,letter);ax.set_title(title,x=.08,ha='left',pad=12)

def curves(ax,cancer):
    x=read(cancer,'framework_independent_benchmark_validation.tsv');source(x,cancer+'_benchmark_aggregate')
    for r in x.itertuples():
        vals=[getattr(r,'top'+str(k)+'_recovery') for k in [20,50,100,200]]
        ax.plot([20,50,100,200],vals,'o-',color=BLUE if r.axis=='translational' else ORANGE,label=r.axis.title())
    ax.set_xlabel('Top-K cutoff');ax.set_ylabel('Held-out recovery');ax.set_ylim(-.03,1.03);style.add_grid(ax);ax.legend(frameon=False,fontsize=7)

def density(ax,score):
    x=score.translational_precedence_score.to_numpy();y=score.functional_discovery_score.to_numpy()
    hist,xe,ye=np.histogram2d(x,y,bins=40)
    ax.pcolormesh(xe,ye,np.log1p(hist.T),cmap='Greys',shading='auto',rasterized=True)
    # Export density bins, not individual score vectors.
    cells=pd.DataFrame({'x_bin':np.repeat(np.arange(40),40),'y_bin':np.tile(np.arange(40),40),'n_records':hist.ravel().astype(int)})
    for i,g in enumerate(EXAMPLES):
        r=score[score.gene_symbol==g].iloc[0]
        ax.scatter(r.translational_precedence_score,r.functional_discovery_score,color=BLUE,s=17,zorder=3)
        xx=float(r.translational_precedence_score);yy=float(r.functional_discovery_score)
        dx=-5 if xx>.65 else 4;dy=-12 if yy>.95 else 5
        if g=='ERBB2':dy=-22
        if g=='PALB2':dy=5
        ax.annotate(g,(xx,yy),xytext=(dx,dy),ha='right' if dx<0 else 'left',textcoords='offset points',fontsize=6.5)
    ax.set_xlim(-.02,1.08);ax.set_ylim(-.02,1.1)
    ax.set_xlabel('Translational similarity');ax.set_ylabel('Functional similarity');return cells

def classes_plot(ax,classes):
    order=['I','II','III','IV','V','VI']
    x=classes.set_index('class').n_records if {'class','n_records'} <= set(classes.columns) else classes.review_class.str.extract(r'Class ([IV]+)')[0].value_counts()
    counts=np.array([x.get(c,0) for c in order]);colors=[style.CLASS_COLORS['Class '+c] for c in order]
    ax.barh(order,counts,color=colors);ax.set_xscale('symlog',linthresh=20);ax.invert_yaxis()
    ax.set_xlabel('Records (log scale)')
    for i,v in enumerate(counts):ax.text(v*1.04+1,i,f'{v:,}',va='center',fontsize=7)
    ax.set_xlim(0,max(counts)*2);style.add_grid(ax,'x')
    return pd.DataFrame({'class':order,'n_records':counts})

def figure3():
    scores,cl,cards=safe('LUAD');fig,axes=plt.subplots(1,3,figsize=(7.4,2.8))
    curves(axes[0],'LUAD');panel(axes[0],'A','Benchmark recovery')
    source(density(axes[1],scores),'fig3_density_bins');panel(axes[1],'B','Dual-axis evidence')
    source(classes_plot(axes[2],cl),'fig3_classes');panel(axes[2],'C','Review classes')
    axes[2].text(.02,.02,'N=35,470',transform=axes[2].transAxes,fontsize=7)
    source(cards,'fig3_original_case_inputs');fig.tight_layout(w_pad=1.6);save(fig,'fig3_luad_positive_control')

def figure4():
    scores,cl,cards=safe('PDAC');fig,axes=plt.subplots(1,3,figsize=(7.4,3.0))
    fields=['patient_alteration','proteomic_support','dependency','druggability','safety_factor','safety_risk']
    x=cards.set_index('gene_symbol').loc[EXAMPLES,fields];hm=axes[0].imshow(x.to_numpy(float),vmin=0,vmax=1,cmap='Blues',aspect='auto')
    axes[0].set_yticks(range(5),EXAMPLES);axes[0].set_xticks(range(6),['Patient','Protein','Depend.','Drug.','Safety\naccept.','Risk'],rotation=55,ha='right')
    axes[0].tick_params(axis='both',labelsize=7);panel(axes[0],'A','Support / risk proxies')
    source(density(axes[1],scores),'fig4_density_bins');panel(axes[1],'B','Dual-axis evidence')
    source(classes_plot(axes[2],cl),'fig4_classes');panel(axes[2],'C','Review classes')
    axes[2].text(.02,.02,'N=40,583',transform=axes[2].transAxes,fontsize=7)
    source(cards,'fig4_original_case_inputs');fig.tight_layout(w_pad=1.5);save(fig,'fig4_pdac_triage')

def figure5():
    perm=pd.concat([read(c,'revision_permutation_summary.tsv').assign(cancer=c) for c in ['PDAC','LUAD']],ignore_index=True)
    order=np.argsort(perm.empirical_p_value.to_numpy(),kind='stable')
    perm['holm8']=0.;perm.loc[order,'holm8']=np.maximum.accumulate(np.minimum(1,perm.empirical_p_value.to_numpy()[order]*np.arange(8,0,-1)))
    source(perm,'fig5_all_permutation_summaries')
    fig,axes=plt.subplots(2,3,figsize=(7.4,6.0))
    for ax,c,letter in [(axes[0,0],'LUAD','A'),(axes[0,1],'PDAC','B')]:
        q=perm[(perm.cancer==c)&(perm['null']=='uniform_gene_labels')].set_index('axis').loc[['translational','functional']]
        y=np.arange(2)
        ax.errorbar(q.permutation_mean,y,xerr=1.96*q.permutation_sd,fmt='o',color=GRAY,capsize=2,ms=4)
        ax.scatter(q.observed_value,y,color=[BLUE,ORANGE],s=28,zorder=3)
        for i,(_,r) in enumerate(q.iterrows()):ax.text(.97,i+.20,f'Holm(8)={r.holm8:.3f}',ha='right',va='center',transform=ax.get_yaxis_transform(),fontsize=6.5)
        ax.set_yticks(y,['Translational','Functional']);ax.set_ylim(1.4,-.4);ax.set_xlim(0,1.05);ax.set_xlabel('Top-50 recovery')
        panel(ax,letter,c+' permutation');style.add_grid(ax,'x')
    q=pd.concat([read(c,'framework_leave_source_out_validation.tsv').assign(cancer=c) for c in ['PDAC','LUAD']])
    q=q[q.status=='ESTIMATED'].copy();source(q,'fig5_source_role_out')
    ax=axes[0,2];y=np.arange(len(q));ax.scatter(q.median_rank_percentile,y,color=[BLUE if a=='translational' else ORANGE for a in q.axis],s=20)
    short={'known_direct_drug':'drug','known_biomarker':'biomarker','known_dependency':'dependency','known_pdac_driver':'driver','luad_actionable_anchor':'anchors','tcga_luad_nature_2014':'TCGA','depmap_luad_high_dependency':'DepMap'}
    ax.set_yticks(y,[r.cancer+' '+short.get(r.heldout_source,str(r.heldout_source).replace('curated_literature:','')) for r in q.itertuples()],fontsize=6.6)
    ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('Median rank percentile');panel(ax,'C','Source-role-out');style.add_grid(ax,'x')
    handles=[];labels=[]
    for c in ['PDAC','LUAD']:
        q=read(c,'framework_independent_benchmark_validation.tsv')
        for r in q.itertuples():
            line,=axes[1,0].plot([20,50,100,200],[getattr(r,'top'+str(k)+'_recovery') for k in [20,50,100,200]],marker='o',ms=3,color=BLUE if r.axis=='translational' else ORANGE,ls='--' if c=='PDAC' else '-')
            handles.append(line);labels.append(c+' '+r.axis[:5]+'.')
        q=read(c,'framework_oof_validation_summary.tsv');source(q,c+'_fig5_oof')
        for r in q.itertuples():axes[1,1].plot([20,50,100,200],[getattr(r,'top'+str(k)+'_recovery') for k in [20,50,100,200]],marker='o',ms=3,color=BLUE if r.axis=='translational' else ORANGE,ls='--' if c=='PDAC' else '-')
    for ax,letter,title,ylabel in [(axes[1,0],'D','Held-out benchmarks','Benchmark recovery'),(axes[1,1],'E','Development positives','OOF recovery')]:
        ax.set_xlabel('Top-K cutoff');ax.set_ylabel(ylabel);ax.set_ylim(-.03,1.03);style.add_grid(ax);panel(ax,letter,title)
    q=pd.concat([read(c,'revision_stability_summary.tsv') for c in ['PDAC','LUAD']]);source(q,'fig5_top50_stability')
    q=q[q.metric=='top50_overlap_fraction'];ax=axes[1,2]
    settings=[('PDAC','translational'),('PDAC','functional'),('LUAD','translational'),('LUAD','functional')]
    for j,(c,a) in enumerate(settings):
        for condition,offset,color in [('background_ensemble',-.12,GRAY),('label_thinning',.12,BLUE if a=='translational' else ORANGE)]:
            r=q[(q.cancer==c)&(q.axis==a)&(q.analysis==condition)].iloc[0]
            ax.errorbar(r.estimate,j+offset,xerr=[[r.estimate-r.range_p025],[r.range_p975-r.estimate]],fmt='o',color=color,ms=4,capsize=2)
    ax.set_yticks(range(4),[c+' '+a[:5]+'.' for c,a in settings],fontsize=7);ax.invert_yaxis();ax.set_xlim(0,1.04);ax.set_xlabel('Top-50 overlap');panel(ax,'F','Background / label stress');style.add_grid(ax,'x')
    fig.legend(handles,labels,ncol=4,loc='lower center',frameon=False,fontsize=7,bbox_to_anchor=(.5,.005));fig.tight_layout(rect=(0,.05,1,1),h_pad=2,w_pad=1.2)
    save(fig,'fig5_ranking_validation')

def figure6():
    _,_,cards=safe('PDAC');cards=cards.set_index('gene_symbol').loc[EXAMPLES[:4]].reset_index()
    data=json.loads((DATA/'context/original_example_measurement_cards.json').read_text(encoding='utf8'))
    meta={r['gene_symbol']:r for r in data if r['cancer']=='PDAC' and r['backend']=='archived_strict_elastic_net'}
    fig,axes=plt.subplots(2,2,figsize=(7.4,5.7))
    fields=[('patient_alteration','Patient'),('proteomic_support','Protein'),('dependency','Depend.'),('druggability','Drug.'),('safety_factor','Safety accept.'),('safety_risk','Safety risk')]
    for ax,r in zip(axes.ravel(),cards.itertuples()):
        ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off');cls=r.review_class.split(':')[0];edge=style.CLASS_COLORS.get(cls,GRAY)
        ax.add_patch(Rectangle((.02,.03),.96,.94,fill=False,ec=edge,lw=1.3))
        ax.text(.06,.91,r.gene_symbol+' (PDAC)',fontsize=10,fontweight='bold')
        ax.text(.06,.82,f'{cls}; trans {r.translational_precedence_rank}; func {r.functional_discovery_rank}',fontsize=7.3)
        for i,(field,label) in enumerate(fields):
            yy=.71-i*.069;v=float(getattr(r,field));ax.text(.06,yy,label,va='center',fontsize=7)
            ax.add_patch(Rectangle((.36,yy-.016),.43,.027,fc='#efefef',ec='none'))
            ax.add_patch(Rectangle((.36,yy-.016),.43*np.clip(v,0,1),.027,fc='#59A14F' if field=='safety_factor' else BLUE,ec='none'))
            ax.text(.83,yy,f'{v:.2f}',va='center',fontsize=7)
        m=meta[r.gene_symbol];ax.text(.06,.225,'Constraint: '+m['constraint_measurement_state'].lower(),fontsize=7)
        ax.text(.06,.155,'Policy classes: '+', '.join(map(str,m['operational_classes_all_policies'])),fontsize=7)
        ax.text(.06,.085,'Source and intervention scope: see case text',fontsize=6.7,color='#606060')
    source(cards,'fig6_original_cases');fig.tight_layout(h_pad=1.5,w_pad=1);save(fig,'fig6_evidence_cards')
def main():
    global DATA,OUTPUT
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--figures',type=int,nargs='+',choices=[3,4,5,6],default=[3,4,5,6],help='Select figure numbers; default: 3 4 5 6')
    args=parser.parse_args();DATA=args.data.resolve();OUTPUT=args.output.resolve();OUTPUT.mkdir(parents=True,exist_ok=True)
    for number in dict.fromkeys(args.figures):
        {3:figure3,4:figure4,5:figure5,6:figure6}[number]()
    files=sorted(SOURCES.glob('*.tsv'))
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (OUTPUT/'figure_input_sha256.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

if __name__=='__main__':main()
