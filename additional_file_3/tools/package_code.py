"""Refresh/verify the code manifest, build Additional file 3, and verify its contents."""
from pathlib import Path
import argparse,hashlib,json,tempfile,zipfile
EXCLUDED={'.git','__pycache__','.venv','.runtime','results','logs'}
def digest(data):return hashlib.sha256(data).hexdigest()
def included(root):return {p.relative_to(root).as_posix():p for p in root.rglob('*') if p.is_file() and not any(x in EXCLUDED for x in p.relative_to(root).parts) and p.suffix not in ['.pyc','.pyo']}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--refresh-manifest',action='store_true')
    parser.add_argument('--output',type=Path,help='New zip path outside the code directory; existing files are never overwritten')
    args=parser.parse_args();root=args.root.resolve();files=included(root)
    actual={name:digest(p.read_bytes()) for name,p in sorted(files.items()) if name!='sha256.json'}
    if args.refresh_manifest:(root/'sha256.json').write_text(json.dumps(actual,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    expected=json.loads((root/'sha256.json').read_text('utf-8-sig'))
    if expected!=actual:raise RuntimeError('Manifest differs from current files; review changes, then use --refresh-manifest')
    if args.output is None:print('Manifest verified: '+str(len(actual))+' files (sha256.json excludes itself)');return
    output=args.output.resolve()
    if output.is_relative_to(root):raise ValueError('Place the zip outside the code directory')
    output.parent.mkdir(parents=True,exist_ok=True);files=included(root)
    with zipfile.ZipFile(output,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name,path in sorted(files.items()):archive.write(path,'additional_file_3/'+name)
    with tempfile.TemporaryDirectory(prefix='af3-check-',dir=output.parent) as temporary:
        with zipfile.ZipFile(output) as archive:
            names=archive.namelist();wanted={'additional_file_3/'+name for name in files}
            if len(names)!=len(set(names)) or set(names)!=wanted:raise RuntimeError('Unexpected or missing zip entries')
            if archive.testzip() is not None:raise RuntimeError('Zip CRC verification failed')
            archive.extractall(temporary)
        extracted=Path(temporary)/'additional_file_3';ziphash={name:digest((extracted/name).read_bytes()) for name in files};sourcehash={name:digest(path.read_bytes()) for name,path in files.items()}
        if ziphash!=sourcehash:raise RuntimeError('Extracted package differs from source directory')
    sha=digest(output.read_bytes());output.with_suffix('.zip.sha256').write_text(sha+'  '+output.name+'\n',encoding='utf8')
    record=dict(zip_name=output.name,zip_sha256=sha,root_prefix='additional_file_3/',n_files=len(files),n_manifest_entries=len(expected),crc_ok=True,extracted_files_identical=True,files=ziphash)
    output.with_suffix('.verification.json').write_text(json.dumps(record,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in record.items() if k!='files'},indent=2))
if __name__=='__main__':main()
