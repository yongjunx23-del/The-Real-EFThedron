"""Check finite-width threshold caps for a new exact saved direction."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/all_direction_certification_20261007/head'))
from finite_fast import *

ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True)
ap.add_argument('--moments',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
a=ap.parse_args();ctx.prec=256;tic=time.time();rows=[]
for ell in [0,2,4,6,8,14,40]:
    ts=time.time();m=FiniteCentered(a.witness,ell,a.moments,N=128) if ell<4 else FiniteFast(a.witness,ell,a.moments,N=128)
    try:
        gaps=m.gaps(30);bounds=[m.gap_bound(lo,hi) for lo,hi in gaps]
        assert all(v.is_finite() and v>=0 for v in bounds)
        r=dict(spin=ell,status='PASS_FINITE_WIDTH_CAPS',gaps=[[str(x),str(y)]for x,y in gaps],
            bound_upper=[str(v.upper())for v in bounds],threshold=m.cap_status,seconds=time.time()-ts)
    except (ValueError,AssertionError,ZeroDivisionError) as e:
        r=dict(spin=ell,status='FAILED_PREFLIGHT',error=str(e),seconds=time.time()-ts)
    rows.append(r);print(json.dumps(r),flush=True)
out=dict(status='PASS_PREFLIGHT'if all(r['status'].startswith('PASS')for r in rows)else'FAILED_PREFLIGHT',
    witness_sha256=hashlib.sha256(a.witness.read_bytes()).hexdigest(),
    source_manifest_sha256=hashlib.sha256((a.moments/'manifest.json').read_bytes()).hexdigest(),
    rows=rows,seconds=time.time()-tic,scope='Finite-width caps for the listed spins only; not a complete support bound.')
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2))
