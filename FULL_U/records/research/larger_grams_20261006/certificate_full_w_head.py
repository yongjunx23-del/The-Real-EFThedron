"""Continuum 2x2 polynomial-PSD certificates for the full unitarity disk.

H is constant on rational velocity cells; h=(1-beta^2)^2 H.
Every Bernstein coefficient matrix after a positive congruence is verified
PSD in exact rational arithmetic. Numerical work only selects H/partitions.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,ceil,sqrt
import json,sys,argparse,time
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/larger_grams_20261006';OUT.mkdir(parents=True,exist_ok=True)
beta,t,L,z=s.symbols('beta t L z')
E=4/(1-beta**2);delta=E-4;Q=E-s.Rational(4,3)
x=(1-beta**2)/(1+beta**2);e2=(1-beta**2)**2
A=2+t/2;qt=E-2+t/2;ff=8*A*A/qt**3;c=2*A*A/qt**2-1
ka=2/qt-1/Q-1/(E-s.Rational(8,3)+t)
kb=1/(E-t)+1/(E-s.Rational(8,3)+t)-2/Q
local=s.diff(ff*c-ka,t).subs(t,0)+L/delta*(ff*c-ka).subs(t,0)+2*ff.subs(t,0)*(2*c.subs(t,0)**2-1)
C=-Q**3*s.diff(kb,t).subs(t,0)
Uraw=16*beta*(1-beta**2)**4/(1+beta**2)**5*(L/delta+(3-x-5*x*x+16*(1+x)*(2*x*x-1))/(4*(1+x)))

def coeff(expr):
    p=s.Poly(s.expand(expr),beta);return [F(p.nth(j))for j in range(p.degree()+1)]

def pad(p,n):return p+[F(0)]*(n-len(p))

def bern_interval(p,lo,hi):
    n=len(p)-1;d=hi-lo
    powerslo=[lo**j for j in range(n+1)];powersd=[d**j for j in range(n+1)]
    mon=[sum(p[j]*comb(j,k)*powerslo[j-k] for j in range(k,n+1))*powersd[k]for k in range(n+1)]
    return [sum(mon[k]*F(comb(j,k),comb(n,k))for k in range(j+1))for j in range(n+1)]

def exact_ceiling(v,places=12):
    v=F(v);m=10**places;return F(-((-v.numerator*m)//v.denominator),m)

def build(ell,lam):
    P=s.legendre(ell,z);Z=Q/delta
    DD=(2*z*z*s.diff(P,z,2)-3*z*s.diff(P,z)).subs(z,Z)
    pp=P.subs(z,Z)/Q**3
    VV=pp-s.Rational(lam)*(local.subs(L,ell*(ell+1))+C*P.subs(z,Z)/Q**3+10*DD/Q**5)
    UU=-s.Rational(lam)*Uraw.subs(L,ell*(ell+1))
    off=s.factor(-beta**ell*UU/e2)
    bb=s.factor(-2*beta**(2*ell)*VV/e2)
    denominator=s.lcm(s.denom(off),s.denom(bb))
    if denominator.subs(beta,s.Rational(1,2))<0:denominator=-denominator
    # All denominator factors have no real roots in [0,1].
    assert denominator.subs(beta,0)>0 and denominator.subs(beta,1)>0
    assert s.Poly(denominator,beta).count_roots(0,1)==0
    arr=[coeff(beta**(2*ell)*denominator),coeff(s.cancel(bb*denominator)),
         coeff(s.cancel(off*denominator)),coeff(denominator)]
    n=max(map(len,arr));arr=[pad(p,n)for p in arr]
    hf=s.lambdify(beta,[s.cancel(UU/e2),s.cancel(VV/e2)],'numpy')
    return arr,hf,str(s.factor(denominator))

def solve_cell(arr,lo,hi):
    arrays=[bern_interval(p,lo,hi)for p in arr]
    needed=0.
    for a,b,c,d in zip(*arrays):
        if a<0 or d<=0:return None
        af,bf,cf,df=map(float,(a,b,c,d))
        if a==0:
            if b<0:return None
            if b==0:
                if c:return None
                root=0.
            else:root=cf*cf/(bf*df)
        else:
            disc=sqrt(bf*bf+4*af*cf*cf/df)
            root=(2*cf*cf/df)/(bf+disc)if bf>=0 and bf+disc>0 else(-bf+disc)/(2*af)
        needed=max(needed,root)
    H=F(ceil((needed*(1+2e-12)+1e-13)*10**12),10**12)
    for attempt in range(8):
        if all(a*H+b>=0 and d*H>=0 and (a*H+b)*(d*H)-c*c>=0 for a,b,c,d in zip(*arrays)):
            return H
        H=H+F(1,10**10)+abs(H)*F(1,10**9)
    return None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--lambda',dest='lam',default='1/4');ap.add_argument('--last-spin',type=int,default=18)
    ap.add_argument('--relative-slack',type=float,default=.035);ap.add_argument('--output',default='certificate_full_w_head.json');args=ap.parse_args()
    lam=F(args.lam);records=[];total=F(0);start=time.time()
    for ell in range(0,args.last_spin+1,2):
        arr,hf,den=build(ell,lam);cells=[];stack=[(F(j,16),F(j+1,16),0)for j in reversed(range(16))]
        while stack:
            lo,hi,depth=stack.pop();H=solve_cell(arr,lo,hi)
            pts=np.array([float(lo+(hi-lo)*F(j,8))for j in range(1,8)])
            uu,vv=hf(pts);rad=np.hypot(uu,vv)
            opt=np.where(vv<0,uu*uu/(rad-vv),rad+vv)
            mean=float(np.mean(opt));acceptable=(H is not None and float(H)<=(1+args.relative_slack)*mean+2e-7)
            if acceptable or (depth>=18 and H is not None):
                cells.append({'lo':str(lo),'hi':str(hi),'H':str(H)})
            else:
                if depth>=20:raise ArithmeticError(('unproved cell',ell,str(lo),str(hi)))
                mid=(lo+hi)/2;stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
        assert cells[0]['lo']=='0' and cells[-1]['hi']=='1'
        assert all(F(a['hi'])==F(b['lo'])for a,b in zip(cells,cells[1:]))
        # 333/106 < pi; using the lower endpoint rounds the objective outward.
        cost=sum((F(r['hi'])-F(r['lo']))*F(r['H'])for r in cells)*F(4*(2*ell+1)*106,333)
        total+=cost
        records.append({'ell':ell,'denominator':den,'polynomials_A_B_C_D':[[str(q)for q in p]for p in arr],
                        'cells':cells,'support_upper_rational':str(cost),'support_upper_decimal':float(exact_ceiling(cost,9))})
        print('spin',ell,'cells',len(cells),'upper',float(cost),'elapsed',time.time()-start,flush=True)
        report={'status':'IN_PROGRESS','direction':[0,1],'names':['Dt1G_regular','H2','X4'],
            'multipliers':[str(-lam),str(-2*lam),str(-10*lam)],'finite_spins':records,
            'finite_head_upper_rational':str(total),'finite_head_upper_decimal':float(exact_ceiling(total,9)),
            'warning':'Finite head only. An additional all-spin tail bound is required.',
            'proof_method':'Exact rational PSD of every 2x2 Bernstein coefficient matrix on gapless rational velocity cells.'}
        (OUT/args.output).write_text(json.dumps(report,indent=2),encoding='utf8')
    report['status']='PASS_FINITE_HEAD';(OUT/args.output).write_text(json.dumps(report,indent=2),encoding='utf8')
    print('TOTAL HEAD',float(total),flush=True)

if __name__=='__main__':main()
