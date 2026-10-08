"""Finite-width threshold components for all even head spins, no head sum."""
from finite220 import *
ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--moments',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--lmax',type=int,default=40)
a=ap.parse_args();ctx.prec=256;tic=time.time();records=[]
for ell in range(2,a.lmax+1,2):
    m=FiniteSpin(a.witness,ell,a.moments,N=128);gaps=m.gaps()
    records.append(dict(ell=ell,cap=m.cap_status));print(ell,m.cap_status,time.time()-tic,flush=True)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(dict(witness_sha256=m.sha,seconds=time.time()-tic,records=records,complete=ell==a.lmax,scope='Finite-width beta threshold caps only; remaining energy and other spins not covered'),indent=2))
