"""Independent direct derivative-kernel comparison for the real220 coefficient field."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import sys,json,argparse,hashlib,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(Path(__file__).parent),str(ROOT/'research')]
from real220 import coefficient_fields,q
from forward_transfer_derivatives import build_derivative_bank
def main():
    p=argparse.ArgumentParser();p.add_argument('--witness',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();ww,groups,rows,alpha=coefficient_fields(a.witness)
    mu=np.array([4.05,4.5,5.,8.,12.,24.,48.,100.,1000.,10000.,1e6]);spins=np.array([0,2,14,42,100,400]);beta=np.sqrt(1-4/mu);omega=((2+1j*np.sqrt(mu*(mu-4)))/(mu-2))**2;ff=32/(mu-2)**3
    direct=np.broadcast_to(float(alpha)*ff*omega.imag,(len(spins),len(mu))).copy()
    for name,v in rows.items():
        if name.startswith('H')and name[1:].isdigit():direct+=float(v)*(ff*omega**int(name[1:])).imag
        if name.startswith('rational_'):
            _,kind,rr,cc=name.split('_');r,c=float(F(rr[1:])),float(F(cc[1:]));den=1-2*r*c*omega+r*r*omega**2
            mm=(1-r*c*omega)/den if kind=='cos'else r*omega/den
            at=(1+r*c)/(1+2*r*c+r*r)if kind=='cos'else-r/(1+2*r*c+r*r)
            direct+=float(v)*(ff*omega**2*(mm-at)).imag
    with np.errstate(over='ignore',invalid='ignore'):
        for name,u,v in build_derivative_bank(mu,spins,3,12):
            if name in rows:direct+=float(rows[name])*u
    errors=[]
    for i,l in enumerate(spins):
        T=int(l*(l+1))
        for j,ss in enumerate(mu):
            be=F(float(beta[j]));zz=be*be;e=1-zz;de=4*zz/e
            coeff=[F(str(w.n(q(be))/w.d(q(be))))for w in ww]
            B=[F(1)]
            for k in range(1,4):B.append(B[-1]*F(T-(k-1)*k,k*k))
            Uref=16*e**4/(be*(1+zz)**5)
            got=Uref*sum(coeff[k]*B[k]/de**k for k in range(4))
            err=abs(float(got)-direct[i,j])/(1+abs(float(got)));errors.append((err,int(l),float(ss),float(got),float(direct[i,j])))
    worst=max(errors);assert worst[0]<2e-8,worst
    result=dict(status='PASS_INDEPENDENT_REAL_KERNEL_COMPARISON',witness_sha256=hashlib.sha256(a.witness.read_bytes()).hexdigest(),points=len(errors),maximum_scaled_error=worst[0],worst_point=list(worst[1:]),scope='Forward/rational/derivative real kernel only; source-local strips bounded separately. This numerical comparison does not itself prove the infinite-spin norm formula.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
