from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from framework.verify import verify_results, write_verification_report

if __name__ == '__main__':
    path = write_verification_report()
    report = verify_results()
    print(f'Wrote verification report: {path}')
    print(report.to_string(index=False))
    raise SystemExit(0 if report[report['status'].ne('PASS')].empty else 1)


