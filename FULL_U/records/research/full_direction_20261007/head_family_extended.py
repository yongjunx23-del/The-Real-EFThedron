"""Validated heads for g0+b*g2, with an exact increment d of the X4 null.

The frozen 150-row witness and source moments are unchanged.  This module
adds exactly b*p+d*X4 to V and nothing to U, including derivative bounds.
Every integration cell uses the existing validated disk-support formula.
"""
from __future__ import annotations
import argparse, functools, hashlib, json, math, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from flint import arb, arb_series, ctx
from certify_real_tail import exact
from head_certification_spin0 import SpinZero
from head_certification_lowspin import LowSpin, pe
C=F(4,3)
WITNESS=ROOT/'output/angles_resummed/angle_000p000.npz'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(q):return F(int(q[0]))*F(2)**int(q[1])
def upper(q):
    assert q.is_finite()
    return dyad(q.upper().man_exp())
def outward(q,digits=12):
    n=-((-q.numerator*10**digits)//q.denominator)
    return f'{n//10**digits}.{n%10**digits:0{digits}d}'

class DirectionalMixin:
    def setup(self,base,slope,crossing):
        self.__dict__=base.__dict__.copy()
        for name in ('fields','fields_series','normalized_interval','threshold_cap'):
            self.__dict__.pop(name,None)
        self.original=base;self.slope=F(slope);self.crossing=F(crossing)
        self.endpoint_cache={}

    def addition(self,a,b,series=False):
        ball=arb(exact((a+b)/2),exact((b-a)/2))
        beta=arb_series([ball,1],prec=3) if series else ball
        x=(1-beta*beta)/4
        ell=getattr(self,'ell',0)
        if ell:
            h=x/(1-4*x)
            pstar=sum(self.coeff[j]*(exact(C)*h)**j for j in range(ell+1))
            z=1+2*exact(C)*h
            cross=x**3/(16*(1-exact(C)*x)**5)*pe(self.crossZ,z)
        else:pstar=1;cross=0
        return exact(self.slope)*x*pstar/(16*(1-exact(C)*x)**3)+exact(self.crossing)*cross

    def fields(self,a,b,omit_sources=()):
        u,v,logs=self.original.fields(a,b,tuple(omit_sources))
        if a==b==1:return u,v,logs
        return u,v+self.addition(a,b),logs

    def fields_series(self,a,b):
        u,v,err=self.original.fields_series(a,b)
        return u,v+self.addition(a,b,True),err

class DirectionalZero(DirectionalMixin,SpinZero):
    def __init__(self,base,slope,crossing):self.setup(base,slope,crossing)

class DirectionalLow(DirectionalMixin,LowSpin):
    def __init__(self,base,slope,crossing):self.setup(base,slope,crossing)

    def ratio(self,a,b):
        beta=arb(exact((a+b)/2),exact((b-a)/2));e=1-beta*beta
        mu=4/e;delta=mu-4;Q=mu-exact(C)
        return pe(self.crossD,delta)/(Q*Q*pe(self.Rstar,delta))

    def normalized_interval(self,a,b):
        u,v,pinv,e=self.original.normalized_interval(a,b)
        return u,v+exact(self.slope)+exact(self.crossing)*self.ratio(a,b),pinv,e

    def threshold_cap(self,b):
        oldcap=self.original.threshold_cap(b)
        status=self.original.cap_status.copy()
        oldv=F(int(status['V_over_p_upper_dyadic'][0]))*F(2)**int(status['V_over_p_upper_dyadic'][1])
        newv=exact(oldv)+exact(self.slope)+exact(self.crossing)*self.ratio(F(0),b)
        if not newv.upper()<0:raise ValueError('Directional threshold sign failed: '+str(newv))
        # The old proof bounded numerator N=measure*U^2/p and used N/(-2Vn).
        # Multiplication by (-oldVn)/(-newVn) transfers that same numerator cap.
        bound=(oldcap*exact(-oldv)/(-newv.upper())).upper()
        self.cap_status=dict(old=status,slope=str(self.slope),crossing_increment=str(self.crossing),
            new_V_over_p_upper_dyadic=list(map(str,newv.upper().man_exp())),upper=str(bound))
        return bound

def verify_dependencies():
    path=ROOT/'output/certification_head/complete_upper_certificate.json'
    old=json.loads(path.read_text());assert old['status']=='PASS'
    for p,h in old['sha256_manifest'].items():
        assert sha(ROOT/p)==h,('Changed prerequisite',p)
    return {**{p:sha(ROOT/p) for p in old['sha256_manifest']},str(path.relative_to(ROOT)):sha(path)}

def run(spin,slopes,precision,out,subdivide=1,repair_ratio=F(18,5)):
    ctx.prec=precision;tic=time.time();deps=verify_dependencies()
    part=ROOT/('output/certification_head/spin0_component_tight.json' if spin==0 else f'output/certification_head/spin{spin}_tight_component.json')
    data=json.loads(part.read_text());deps[str(part.relative_to(ROOT))]=sha(part)
    base=SpinZero(WITNESS) if spin==0 else LowSpin(WITNESS,spin)
    base.fields=functools.lru_cache(maxsize=50000)(base.fields)
    base.fields_series=functools.lru_cache(maxsize=2)(base.fields_series)
    if spin:
        base.normalized_interval=functools.lru_cache(maxsize=2)(base.normalized_interval)
        threshold=base.threshold_cap
        def saved_threshold(b):
            if b not in thresholds:
                q=threshold(b);thresholds[b]=(q,base.cap_status.copy())
            q,status=thresholds[b];base.cap_status=status.copy();return q
        thresholds={};base.threshold_cap=saved_threshold
    specs=[(F(s),F(repair_ratio)*F(s) if F(s)>0 else F(0)) for s in slopes]
    models=[(DirectionalZero if spin==0 else DirectionalLow)(base,b,d) for b,d in specs]
    gaps=[tuple(map(F,x)) for x in data['gap_beta_intervals']]
    cells=[]
    for c in data['cells_data']:
        a,b=F(c['a']),F(c['b'])
        for i in range(subdivide):cells.append((a+(b-a)*i/subdivide,a+(b-a)*(i+1)/subdivide))
    partition=sorted(cells+gaps)
    assert partition[0][0]==0 and partition[-1][1]==1
    assert all(a<b for a,b in partition)
    assert all(b==c for (a,b),(c,d) in zip(partition,partition[1:]))
    totals=[F(0) for _ in specs];records=[[] for _ in specs];gap_records=[[] for _ in specs]
    for a,b in gaps:
        for j,m in enumerate(models):
            val=upper(m.gap_bound(a,b));assert val>=0;totals[j]+=val
            gap_records[j].append({'a':str(a),'b':str(b),'upper_exact':str(val)})
    print(json.dumps({'event':'built','spin':spin,'seconds':time.time()-tic,'directions':list(slopes)}),flush=True)
    for i,(a,b) in enumerate(cells):
        for j,m in enumerate(models):
            q,_=m.trapezoid(a,b);val=upper(q);assert val>=0;totals[j]+=val
            records[j].append({'a':str(a),'b':str(b),'upper_exact':str(val)})
        if (i+1)%1000==0:print(json.dumps({'event':'progress','spin':spin,'cells':i+1,'seconds':time.time()-tic}),flush=True)
    result={'status':'PASS','claim':'Directional single-spin head support upper bounds; finite block and infinite tail are separate.',
        'spin':spin,'precision_bits':precision,'witness_sha256':sha(WITNESS),'source_sha256':sha(__file__),
        'partition_exactly_covers_beta_0_1':True,'repair_ratio_exact':str(repair_ratio),'dependency_sha256':deps,'seconds':time.time()-tic,
        'directions':[{'slope_exact':str(b),'crossing_increment_exact':str(d),'head_upper_exact':str(totals[j]),
          'head_upper_outward_12dp':outward(totals[j]),'cells':records[j],'gaps':gap_records[j],
          'threshold':models[j].cap_status if spin else None} for j,(b,d) in enumerate(specs)]}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2))
    print(json.dumps({'event':'complete','spin':spin,'upper':[outward(q) for q in totals],'seconds':time.time()-tic}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--spin',type=int,choices=[0,2,4,6],required=True)
    ap.add_argument('--slopes',default='-1,-1/2,-1/4,-1/10,0,1/10,1/4,1/2')
    ap.add_argument('--precision',type=int,default=256);ap.add_argument('--subdivide',type=int,default=1)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--repair-ratio',default='18/5');a=ap.parse_args()
    run(a.spin,a.slopes.split(','),a.precision,a.out,a.subdivide,F(a.repair_ratio))
