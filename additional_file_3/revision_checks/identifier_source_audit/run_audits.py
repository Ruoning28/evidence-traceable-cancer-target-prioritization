from pathlib import Path
import argparse,subprocess,sys
parser=argparse.ArgumentParser(description='Read-only identifier and source-exclusion audits.')
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
here=Path(__file__).resolve().parent
code=here.parents[1]
data=args.data.resolve()
ref=data/'revision_checks/identifier_source_audit/reference'
for name in ['audit_identifiers.py','check_background_samples.py','check_review_denominators.py','audit_source_paths.py']:
    subprocess.run([sys.executable,str(here/name),'--data',str(data),'--reference',str(ref),
                    '--code',str(code),'--output',str(args.output.resolve())],check=True)
