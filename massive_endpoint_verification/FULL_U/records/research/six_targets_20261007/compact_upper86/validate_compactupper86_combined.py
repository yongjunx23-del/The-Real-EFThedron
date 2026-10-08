from pathlib import Path
import sys
_ADAPT_ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(_ADAPT_ROOT/'research/six_targets_20261007/tail_upper86'),str(_ADAPT_ROOT/'research/boundary_length_20261007/compact')]
"""Independent floating kernel comparison for the new compact evaluator."""
from compact_finiteupper86_combined import *
from angle_audit import fields
import numpy as np,json,time,argparse
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--spins',default='42');ap.add_argument('--moments',type=Path,default=ROOT/'output/six_targets_20261007/upper86/compact/moments');a=ap.parse_args()
    ctx.prec=512;started=time.time();spins=list(map(int,a.spins.split(',')));bank=CompactFinite220(a.moments,spins);print('BUILT',time.time()-started,flush=True)
    ss=np.array([4.0001,4.02,4.5,5.,7.,10.,11.8,12.001,18.,24.001,36.,48.001,49.]);st=np.load(WITNESS)
    if max(spins)>60:ss=ss[ss>=4.5]
    uu,vv,pp=fields(st['names'],st['multipliers'],st['direction'],ss,spins,kernel_order=128,angular_order=128);records=[]
    for i,mu in enumerate(ss):
        values=bank.fields(ar(F(str(mu))))
        for li,l in enumerate(spins):
            v=values[l];r=dict(mu=float(mu),ell=l,interval=str(v),numeric=float(vv[li,i]),relative=abs(float(v)-vv[li,i])/(1+abs(vv[li,i])))
            records.append(r);print(r,flush=True)
    out=ROOT/'output/six_targets_20261007/upper86/compact'/('point_comparison_combined_'+a.spins.replace(',','_')+'.json')
    assert max(r['relative']for r in records)<1e-8
    out.write_text(json.dumps(dict(records=records,seconds=time.time()-started,scope='Independent numerical comparison only, not an interval domain certificate'),indent=2))
