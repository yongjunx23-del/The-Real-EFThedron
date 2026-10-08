"""Exact full-W Bernstein certificates for rational first-derivative witnesses."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,time
import sympy as s
import numpy as np
import certificate_full_w_head as base
beta,L,t,z=base.beta,base.L,base.t,base.z
E,delta,Q,e2,x=base.E,base.delta,base.Q,base.e2,base.x
U0=16*beta*(1-beta**2)**4/(1+beta**2)**5
f0=base.ff.subs(t,0);c0=base.c.subs(t,0)
fp=s.diff(base.ff,t).subs(t,0);cp=s.diff(base.c,t).subs(t,0)
kap=s.diff(base.ka,t).subs(t,0);kbp=s.diff(base.kb,t).subs(t,0)

def build_witness(ell,names,weights):
    lv=ell*(ell+1);P=s.legendre(ell,z);Z=Q/delta;pp=P.subs(z,Z)/Q**3
    U=0;V=pp
    for name,weight in zip(names,weights):
        weight=s.Rational(weight)
        if name=='X4':
            D=(2*z*z*s.diff(P,z,2)-3*z*s.diff(P,z)).subs(z,Z)
            V+=weight*D/Q**5;continue
        deriv=False;sub=False
        if name=='Dt1G_regular':j=1;deriv=True;sub=True
        elif name.startswith('Dt1H_j'):j=int(name[6:]);deriv=True
        elif name.startswith('H'):j=int(name[1:])
        else:raise ValueError('Unsupported exact polynomial row '+name)
        if not sub and j<2:raise ValueError('Only the proved harmonic nulls j>=2 are allowed')
        tj=s.chebyshevt(j,z);uj=s.chebyshevu(j-1,z)
        if deriv:
            uu=U0*(lv/delta*uj.subs(z,c0)+((2-x-3*x*x)*uj.subs(z,c0)-j*tj.subs(z,c0))/(4*(1+x)))
            vv=(fp+lv/delta*f0)*tj.subs(z,c0)+f0*cp*s.diff(tj,z).subs(z,c0)
            if sub:vv-=kap+lv/delta*base.ka.subs(t,0)+P.subs(z,Z)*kbp
        else:uu=U0*uj.subs(z,c0);vv=f0*tj.subs(z,c0)
        U+=weight*uu;V+=weight*vv
    off=s.factor(-beta**ell*U/e2);bb=s.factor(-2*beta**(2*ell)*V/e2)
    denominator=s.lcm(s.denom(off),s.denom(bb))
    if denominator.subs(beta,s.Rational(1,2))<0:denominator=-denominator
    assert denominator.subs(beta,0)>0 and denominator.subs(beta,1)>0
    assert s.Poly(denominator,beta).count_roots(0,1)==0
    arr=[base.coeff(beta**(2*ell)*denominator),base.coeff(s.cancel(bb*denominator)),
         base.coeff(s.cancel(off*denominator)),base.coeff(denominator)]
    n=max(map(len,arr));arr=[base.pad(p,n)for p in arr]
    hf=s.lambdify(beta,[s.cancel(U/e2),s.cancel(V/e2)],'numpy')
    return arr,hf,str(s.factor(denominator))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--first-spin',type=int,default=0);ap.add_argument('--last-spin',type=int,default=18)
    ap.add_argument('--relative-slack',type=float,default=.05);ap.add_argument('--output',required=True);args=ap.parse_args()
    raw=args.witness.read_bytes();witness=json.loads(raw);names=witness['names']
    weights=[F(str(q))for q in witness.get('multipliers_exact',witness.get('multipliers'))]
    if 'direction'in witness:assert witness['direction']==[0,1]
    records=[];total=F(0);start=time.time()
    for ell in range(args.first_spin,args.last_spin+1,2):
        arr,hf,den=build_witness(ell,names,weights);cells=[];stack=[(F(j,16),F(j+1,16),0)for j in reversed(range(16))]
        while stack:
            lo,hi,depth=stack.pop();H=base.solve_cell(arr,lo,hi)
            pts=np.array([float(lo+(hi-lo)*F(j,8))for j in range(1,8)])
            uu,vv=hf(pts);rad=np.hypot(uu,vv);opt=np.where(vv<0,uu*uu/(rad-vv),rad+vv)
            mean=float(np.mean(opt));acceptable=(H is not None and float(H)<=(1+args.relative_slack)*mean+2e-7)
            if acceptable or(depth>=18 and H is not None):cells.append({'lo':str(lo),'hi':str(hi),'H':str(H)})
            else:
                if depth>=20:raise ArithmeticError(('unproved cell',ell,str(lo),str(hi)))
                mid=(lo+hi)/2;stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
        assert cells[0]['lo']=='0' and cells[-1]['hi']=='1'
        assert all(F(a['hi'])==F(b['lo'])for a,b in zip(cells,cells[1:]))
        cost=sum((F(r['hi'])-F(r['lo']))*F(r['H'])for r in cells)*F(4*(2*ell+1)*106,333);total+=cost
        records.append({'ell':ell,'denominator':den,'polynomials_A_B_C_D':[[str(q)for q in p]for p in arr],
            'cells':cells,'support_upper_rational':str(cost),'support_upper_decimal':float(base.exact_ceiling(cost,9))})
        print('spin',ell,'cells',len(cells),'upper',float(cost),'elapsed',time.time()-start,flush=True)
        report={'status':'IN_PROGRESS','direction':[0,1],'names':names,'multipliers':list(map(str,weights)),
            'witness_sha256':hashlib.sha256(raw).hexdigest(),'witness_path':str(args.witness.resolve()),
            'finite_spins':records,'finite_head_upper_rational':str(total),'finite_head_upper_decimal':float(base.exact_ceiling(total,9)),
            'warning':'Finite head only. An additional all-spin tail bound is required.',
            'proof_method':'Exact rational PSD of every 2x2 Bernstein coefficient matrix on gapless rational velocity cells.'}
        (base.OUT/args.output).write_text(json.dumps(report,indent=2),encoding='utf8')
    report['status']='PASS_FINITE_HEAD';(base.OUT/args.output).write_text(json.dumps(report,indent=2),encoding='utf8')
    print('TOTAL HEAD',float(total),flush=True)

if __name__=='__main__':main()
