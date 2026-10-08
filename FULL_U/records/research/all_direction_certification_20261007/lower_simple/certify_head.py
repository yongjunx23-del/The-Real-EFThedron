"""Exact continuum matrix certificate for the new small lower-g0 witness."""
from pathlib import Path
from fractions import Fraction as F
from math import ceil,sqrt
import sys,json,time,hashlib
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'research/larger_grams_20261006'),str(ROOT/'research/all_direction_certification_20261007/alternative')]
import certificate_general_head as gen
import certificate_full_w_head as base
from exact_bernstein import bernstein,q
from rational_head import build as fast_build
from flint import arb,ctx
OUT=ROOT/'output/all_direction_certification_20261007/lower_simple'
def build(ell,names,raw,a,b):
    beta,t,z=base.beta,base.t,base.z;E,delta,Q,e2=base.E,base.delta,base.Q,base.e2
    P=s.legendre(ell,z);Z=Q/delta;pp=P.subs(z,Z)/Q**3;lv=ell*(ell+1)
    CP=8/((3*E-8)*(3*E-4))*(1/(E-2)-4*P.subs(z,Z)/E)
    U=s.Rational(a)*gen.U0;V=s.Rational(a)*(gen.f0*gen.c0+CP)+s.Rational(b)*pp
    for name,rw in zip(names,raw):
        wt=-s.Rational(rw)
        if name=='X4':V+=wt*(2*z*z*s.diff(P,z,2)-3*z*s.diff(P,z)).subs(z,Z)/Q**5;continue
        deriv=name.startswith('Dt');sub=name=='Dt1G_regular';j=1 if sub else int(name[6:]if deriv else name[1:])
        tj=s.chebyshevt(j,z);uj=s.chebyshevu(j-1,z)
        if deriv:
            uu=gen.U0*(lv/delta*uj.subs(z,gen.c0)+((2-base.x-3*base.x**2)*uj.subs(z,gen.c0)-j*tj.subs(z,gen.c0))/(4*(1+base.x)))
            vv=(gen.fp+lv/delta*gen.f0)*tj.subs(z,gen.c0)+gen.f0*gen.cp*s.diff(tj,z).subs(z,gen.c0)
            if sub:vv-=gen.kap+lv/delta*base.ka.subs(t,0)+P.subs(z,Z)*gen.kbp
        else:uu=gen.U0*uj.subs(z,gen.c0);vv=gen.f0*tj.subs(z,gen.c0)
        U+=wt*uu;V+=wt*vv
    off=s.factor(-beta**ell*U/e2);bb=s.factor(-2*beta**(2*ell)*V/e2)
    den=s.lcm(s.denom(off),s.denom(bb))
    if den.subs(beta,s.Rational(1,2))<0:den=-den
    assert den.subs(beta,0)>0 and den.subs(beta,1)>0 and s.Poly(den,beta).count_roots(0,1)==0
    arr=[base.coeff(beta**(2*ell)*den),base.coeff(s.cancel(bb*den)),base.coeff(s.cancel(off*den)),base.coeff(den)]
    n=max(map(len,arr));arr=[base.pad(p,n)for p in arr]
    hf=s.lambdify(beta,[s.cancel(U/e2),s.cancel(V/e2)],'numpy')
    return arr,hf,str(s.factor(den))
def solve(arr,lo,hi):
    arrays=[bernstein(p,lo,hi)for p in arr];needed=0.
    for a,b,c,d in zip(*arrays):
        if a<0 or d<=0:return None
        af,bf,cf,df=map(arb,(a,b,c,d))
        if a==0:
            if b<0 or (b==0 and c):return None
            root=cf*cf/(bf*df)if b else arb(0)
        else:
            disc=(bf*bf+4*af*cf*cf/df).sqrt()
            root=(2*cf*cf/df)/(bf+disc)if b>=0 and bf+disc>0 else(-bf+disc)/(2*af)
        needed=max(needed,float(root.upper()))
    H=F(ceil((needed*(1+2e-12)+1e-13)*10**12),10**12)
    for _ in range(8):
        HH=q(H)
        if all(a*HH+b>=0 and d*HH>=0 and (a*HH+b)*(d*HH)>=c*c for a,b,c,d in zip(*arrays)):return H
        H+=F(1,10**10)+abs(H)*F(1,10**9)
    return None
def main():
    ctx.prec=192
    src=OUT/'candidate.json';data=json.loads(src.read_text());names=data['names'];raw=list(map(F,data['raw_weights_exact']));a,b=map(F,data['direction'])
    start=time.time();records=[];total=F(0);lastspin=40
    if (OUT/'head_certificate.json').exists():
        old=json.loads((OUT/'head_certificate.json').read_text());assert old['source_sha256']==hashlib.sha256(src.read_bytes()).hexdigest()
        records=old['finite_spins'];total=F(old['finite_head_upper_rational'])
    for ell in range(len(records)*2,lastspin+1,2):
        arr,den=fast_build(ell,names,[-v for v in raw],a,b);cells=[];stack=[(F(j,16),F(j+1,16),0)for j in reversed(range(16))]
        while stack:
            lo,hi,depth=stack.pop();H=solve(arr,lo,hi)
            mean=float(solve(arr,(lo+hi)/2,(lo+hi)/2))
            if H is not None and (float(H)<=1.025*mean+2e-7 or depth>=18):cells.append(dict(lo=str(lo),hi=str(hi),H=str(H)))
            else:
                if depth>=20:raise ArithmeticError(('unproved cell',ell,lo,hi,H))
                mid=(lo+hi)/2;stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
        assert cells[0]['lo']=='0'and cells[-1]['hi']=='1' and all(F(r['hi'])==F(t['lo'])for r,t in zip(cells,cells[1:]))
        cost=sum((F(r['hi'])-F(r['lo']))*F(r['H'])for r in cells)*F(4*(2*ell+1)*106,333);total+=cost
        records.append(dict(ell=ell,denominator=den,polynomials_A_B_C_D=[[str(v)for v in p]for p in arr],cells=cells,support_upper_rational=str(cost),support_upper_decimal=float(cost)))
        result=dict(status='PASS_FINITE_HEAD'if ell==lastspin else 'IN_PROGRESS',direction=[str(a),str(b)],names=names,raw_weights_exact=list(map(str,raw)),finite_spins=records,finite_head_upper_rational=str(total),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.time()-start)
        (OUT/'head_certificate.json').write_text(json.dumps(result,indent=2));print('spin',ell,'cells',len(cells),'upper',float(cost),'sum',float(total),'seconds',time.time()-start,flush=True)
if __name__=='__main__':main()
