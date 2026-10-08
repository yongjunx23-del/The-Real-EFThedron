"""Produce per-spin real masses and spectral V/X4 enclosures, read-only inputs.

The fixed finite-spin kernel evaluator and original partitions are reused.
Per-spin data permit sharper local coercivity and affine directional changes.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time,argparse
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research/tail_certification'),str(ROOT/'research/directional_certification'),str(ROOT/'tmp/certification_deps')]
from certified_finite_v import FiniteV,arb,ctx,ar,tight,rev,C,EXPECTED
from certify_finite_real import FiniteReal,exact,power,square
from certify_directional_tail import CASES

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dyad(v):
    assert v.is_finite();return list(map(str,v.upper().man_exp()))
def frac(q):return F(int(q[0]))*F(2)**int(q[1])

def dependencies():
    rp='output/certification_head/V8_40_kappa14_local_replay.json'
    rr=json.loads((ROOT/rp).read_text());assert rr['status']=='PASS' and rr['witness_sha256']==EXPECTED
    for p,h in rr['sha256_manifest'].items():assert sha(p)==h
    hashes=dict(rr['sha256_manifest']);hashes[rp]=sha(rp)
    for p in ['research/directional_certification/certify_finite_vectors.py','research/directional_certification/certify_directional_tail.py']:
        hashes[p]=sha(p)
    return hashes

class VectorReal(FiniteReal):
    def densities(self,a,b):
        beta=arb(exact((a+b)/2),exact((b-a)/2));x=beta*beta;one=1-x
        xb=exact(b)*exact(b);zb=(2+xb)/(3*xb)
        inv={l:(1/zb.legendre_p(l)).upper() for l in self.spins}
        sa=exact(4)/(1-exact(a)**2);sb=exact(4)/(1-exact(b)**2) if b<1 else None
        active=[(A,B,terms) for (A,B),terms in self.groups.items() if not(sa>=exact(B) or (sb is not None and sb<=exact(A)))]
        kval=self.K(beta)
        if not active:
            common=exact(65536)/(27*arb.pi())*x*(2+x)**3/(1+x)**10*power(one,3)*square(kval)
            return {l:(common*(2*l+1)*inv[l]).upper() for l in self.spins}
        assert b<1
        mu=4/one;Q=mu-exact(F(4,3));forward=16*beta*power(one,4)/(1+x)**5*kval
        us={l:forward for l in self.spins}
        for A,B,terms in active:
            clipped=mu.intersection(exact(A).union(exact(B)))
            en=(2*(clipped-exact(A))/exact(B-A)-1).intersection(arb(0,1));taper=1-square(en)
            aa=exact(F(4,3))/(clipped-4)
            pars={j:tight(self.energy[A,B,j],en) for j in {j for j,n,v in terms}}
            for l in self.spins:
                val=arb.pi()/exact(B-A)*taper*sum((pars[j]*tight(self.hp[l,j],aa) for j in pars),arb(0))
                if sa<exact(A) or (sb is not None and sb>exact(B)):val=val.union(arb(0))
                us[l]+=val
        common=Q**3*8/(2*arb.pi()*power(one,2))
        return {l:(common*(2*l+1)*square(us[l])*inv[l]).upper() for l in self.spins}

def run_real():
    hashes=dependencies();ap='output/certification_tail/finite_A8_40_tighter.json'
    aa=json.loads((ROOT/ap).read_text());assert aa['witness_sha256']==EXPECTED;hashes[ap]=sha(ap)
    bank=VectorReal(ROOT/'output/angles_resummed/angle_000p000.npz');records=[];totals={l:arb(0) for l in bank.spins};tic=time.time()
    for i,c in enumerate(aa['partition']):
        a,b=F(c['a']),F(c['b']);values=bank.densities(a,b);row={}
        for l,h in values.items():
            mass=(ar(b-a)*h).upper();assert mass.is_finite() and mass>=0
            row[str(l)]=dyad(mass);totals[l]+=mass
        records.append(dict(a=str(a),b=str(b),mass_upper_by_spin=row))
        if (i+1)%2000==0:print(json.dumps(dict(kind='real',cells=i+1,seconds=time.time()-tic)),flush=True)
    return dict(status='PASS',kind='real_norm',witness_sha256=EXPECTED,spins=bank.spins,precision_bits=ctx.prec,seconds=time.time()-tic,total_upper_by_spin={str(l):dyad(v) for l,v in totals.items()},partition=records,sha256_manifest=hashes)

def run_v():
    hashes=dependencies();vp='output/certification_head/V8_40_kappa14_exact.json'
    vv=json.loads((ROOT/vp).read_text());assert vv['witness_sha256']==EXPECTED;hashes[vp]=sha(vp)
    bank=FiniteV();records=[];tic=time.time();stack=[(c['variable'],F(c['a']),F(c['b']),0) for c in vv['leaves']]
    trials=0
    while stack:
        kind,a,b,depth=stack.pop();trials+=1;x=arb(ar((a+b)/2),ar((b-a)/2));row={};valid=True
        values=bank.fields(x) if kind=='mu' else bank.uv_fields(x)
        for l,v in values.items():
            if kind=='mu':
                delta=x-4;rs=tight(bank.roy.polys[l][1],delta)
                cross=tight(bank.cross[l],delta)/((x-ar(C))*rs)
            else:
                rs=tight(rev(bank.roy.polys[l][1]),x)
                cross=x*x*tight(rev(bank.cross[l],l-1),x)/((1+2*ar(C)*x)*rs)
            if not(v.is_finite() and cross.is_finite()):valid=False;break
            vu,xu=frac(dyad(v)),frac(dyad(cross))
            if not all(vu+bb+dd*xu<0 for bb,dd in CASES):valid=False;break
            row[str(l)]={'V_upper':dyad(v),'X4_upper':dyad(cross)}
        if not valid:
            assert depth<25,(kind,a,b,l,v,cross)
            mid=(a+b)/2;stack.extend([(kind,a,mid,depth+1),(kind,mid,b,depth+1)]);continue
        records.append(dict(variable=kind,a=str(a),b=str(b),upper_by_spin=row))
        if len(records)%1000==0:print(json.dumps(dict(kind='v',cells=len(records),remaining=len(stack),seconds=time.time()-tic)),flush=True)
    return dict(status='PASS',kind='V_and_crossing',witness_sha256=EXPECTED,spins=bank.spins,precision_bits=ctx.prec,seconds=time.time()-tic,certified_cases=[dict(b=str(b),delta_X4=str(d)) for b,d in CASES],acceptance='Every spin finite; exact saved V upper+b+d times exact saved X4 upper is strictly negative for every listed case.',partition=records,sha256_manifest=hashes)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('kind',choices=['real','v']);args=ap.parse_args();ctx.prec=280
    data=run_real() if args.kind=='real' else run_v()
    for p,h in data['sha256_manifest'].items():assert sha(p)==h
    out=ROOT/f'output/directional_certification/finite_{args.kind}_vectors.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(data,indent=2))
    print(json.dumps({k:v for k,v in data.items() if k not in ['partition','sha256_manifest']},indent=2))

if __name__=='__main__':main()
