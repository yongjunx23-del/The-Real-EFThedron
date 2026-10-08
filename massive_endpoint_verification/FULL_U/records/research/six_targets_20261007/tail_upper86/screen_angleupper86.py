"""Feasibility screen only; the energy/spin continua require later proofs."""
from angleupper86_models import *
from finite_uv_bridge220 import ctx,upper
import time,argparse

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lmin',type=int,default=450);ap.add_argument('--eta-panels',type=int,default=16);ap.add_argument('--precision',type=int,default=320);a=ap.parse_args();ctx.prec=a.precision;tic=time.time();source_manifest();bank=Uniform(a.lmin);records=[]
 for tau in map(F,['0','1/20','1/10','1/5','1','3/2','2','9/4','5/2','3','4','8','16','21','24','32','64','128','256']):
  values=[bank.field(tau,tau,F(i,a.eta_panels*bank.T),F(i+1,a.eta_panels*bank.T))for i in range(a.eta_panels)]
  u=max(upper(v)for v in values)if all(v.is_finite()for v in values)else None
  records.append(dict(tau_exact=str(tau),uniform_upper_exact=str(u)if u is not None else None));print(str(tau),float(u)if u is not None else'nonfinite',flush=True)
 finite=Finite();rows=[]
 for ell in[42,80,200,448]:
  for tau in map(F,['2','9/4','5/2','3','21']):
   v=finite.field(ell,tau,tau);rows.append(dict(spin=ell,tau_exact=str(tau),upper_exact=str(upper(v))));print('finite',ell,str(tau),float(upper(v)),flush=True)
 report=dict(status='NUMERICAL_INTERVAL_SCREEN_NOT_CONTINUUM_CERTIFICATE',witness_sha256=SHA,objective_exact=list(map(str,OBJECTIVE)),lmin=a.lmin,eta_panels=a.eta_panels,uniform_point_records=records,finite_point_records=rows,seconds=time.time()-tic)
 (OUT/f'screen_L{a.lmin}_eta{a.eta_panels}_{a.precision}.json').write_text(json.dumps(report,indent=2))

if __name__=='__main__':main()
