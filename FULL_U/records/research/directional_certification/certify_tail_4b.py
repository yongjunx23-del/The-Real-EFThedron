"""Run the reusable tail family with the agreed positive-direction repair d=4b."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
import certify_directional_tail as base

base.CASES=[(b,F(0) if b<=0 else 4*b) for b,d in base.CASES]
if __name__=='__main__':
    if '--output' not in sys.argv:
        sys.argv.extend(['--output',str(base.ROOT/'output/directional_certification/tail_family_4b.json')])
    out=Path(sys.argv[sys.argv.index('--output')+1])
    base.main()
    data=json.loads(out.read_text());data['positive_repair']='delta_X4=4*b for b>0; zero otherwise'
    here=Path(__file__).resolve();data['sha256_manifest'][str(here.relative_to(base.ROOT))]=hashlib.sha256(here.read_bytes()).hexdigest()
    out.write_text(json.dumps(data,indent=2))
