"""Run frozen-table figures in their separately locked Windows runtimes."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,tempfile,zipfile

ROOT=Path(__file__).resolve().parents[1]
PROBE=r'''
import sys,os,json,hashlib,platform
from pathlib import Path
if os.name=='nt':
 p=Path(sys.prefix)/'Library/bin'
 handle=os.add_dll_directory(str(p)) if p.exists() else None
import matplotlib,numpy,pandas
import matplotlib.ft2font as ft
import matplotlib.font_manager as fm
font=Path(fm.findfont(fm.FontProperties(family='Arial'),fallback_to_default=False))
print(json.dumps(dict(python=platform.python_version(),platform=sys.platform,
 matplotlib=matplotlib.__version__,numpy=numpy.__version__,pandas=pandas.__version__,
 freetype=ft.__freetype_version__,freetype_build=ft.__freetype_build_type__,
 arial_sha256=hashlib.sha256(font.read_bytes()).hexdigest(),matplotlib_file=matplotlib.__file__)))
'''

def environment(cache,figure2_runtime=None):
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env['PYTHONUTF8']='1';env['MPLCONFIGDIR']=str(cache)
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[name]='1'
    env.pop('PYTHONPATH',None)
    if figure2_runtime is not None:env['PYTHONPATH']=str(figure2_runtime)
    return env

def verify_profile(name,env,expected):
    result=subprocess.run([sys.executable,'-c',PROBE],env=env,text=True,encoding='utf8',capture_output=True,check=True)
    actual=json.loads(result.stdout)
    mismatches={key:{'expected':value,'actual':actual.get(key)} for key,value in expected.items() if actual.get(key)!=value}
    if mismatches:raise RuntimeError('Pixel-profile mismatch for '+name+': '+json.dumps(mismatches)+'. See figure_code/README.md; no figure was generated for this profile.')
    return actual

def input_archive(data,destination):
    rows=[]
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for cancer in ['PDAC','LUAD']:
            for name in ['framework_feature_matrix.tsv','framework_record_role_manifest.tsv','source_target_evidence_matrix.tsv']:
                contents=(data/'frozen_inputs'/cancer/name).read_bytes();entry='additional_file_2/frozen_inputs/'+cancer+'/'+name
                archive.writestr(entry,contents);rows.append(dict(entry=entry,sha256=hashlib.sha256(contents).hexdigest()))
    return rows

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True,help='Extracted Additional file 2 directory')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--figures',type=int,nargs='+',choices=[2,3,4,5,6],default=[2,3,4,5,6])
    parser.add_argument('--figure2-runtime',type=Path,default=ROOT/'.runtime/figure2')
    args=parser.parse_args();data=args.data.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    chosen=list(dict.fromkeys(args.figures));profiles=json.loads((ROOT/'environments/plotting-profiles.json').read_text('utf8'));checks={};environments={}
    if 2 in chosen:
        target=args.figure2_runtime.resolve()
        if not (target/'matplotlib/__init__.py').is_file():raise FileNotFoundError('Install the hashed Figure 2 wheel into '+str(target)+' as described in figure_code/README.md')
        environments['figure2']=environment(output/'.mpl-figure2',target)
        checks['figure2']=verify_profile('figure2',environments['figure2'],profiles['figure2'])
        if not Path(checks['figure2']['matplotlib_file']).resolve().is_relative_to(target):raise RuntimeError('Figure 2 loaded Matplotlib outside the specified runtime')
    if any(n>=3 for n in chosen):
        environments['figures3_6']=environment(output/'.mpl-figures3-6')
        checks['figures3_6']=verify_profile('figures3_6',environments['figures3_6'],profiles['figures3_6'])
    record=dict(figures=chosen,profiles=checks,inputs=dict(data=str(data)),validation_scope='Runtime/font preflight and frozen-table figure generation; no model refitting.')
    # Check every selected runtime before writing any figure.
    if 2 in chosen:
        with tempfile.TemporaryDirectory(prefix='figure2-input-',dir=output) as temporary:
            archive=Path(temporary)/'frozen_tables.zip';record['figure2_input_files']=input_archive(data,archive)
            subprocess.run([sys.executable,str(ROOT/'figure_code/generate_figure2.py'),'--data-archive',str(archive),'--output',str(output)],env=environments['figure2'],check=True)
    remaining=[n for n in chosen if n>=3]
    if remaining:subprocess.run([sys.executable,str(ROOT/'figure_code/generate_figures.py'),'--data',str(data),'--output',str(output),'--figures',*[str(n) for n in remaining]],env=environments['figures3_6'],check=True)
    (output/'figure_runtime_record.json').write_text(json.dumps(record,indent=2),encoding='utf8')
    print('Generated figures '+', '.join(map(str,chosen))+'; runtime record: '+str(output/'figure_runtime_record.json'))
if __name__=='__main__':main()
