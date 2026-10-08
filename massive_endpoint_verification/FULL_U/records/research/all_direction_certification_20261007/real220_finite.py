"""Finite-spin real norm retaining exact spin-polynomial cancellations."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import sys,json,time,argparse,heapq,hashlib
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent))
from real220 import Norm,coefficient_fields,R,fmpq_poly,q,arb,arb_poly,ctx,exact,square,tight,endpoint
def diff(r):return R(r.n.derivative()*r.d-r.n*r.d.derivative(),r.d*r.d)
def pair(r):return tuple(arb_poly([exact(F(str(v)))for v in p.coeffs()])for p in [r.n,r.d])
def val(p,x):return tight(p[0],x)/tight(p[1],x)
class Finite(Norm):
    def __init__(self,path,first=42,last=400,K=20):
        assert first % 2 == last % 2 == 0 and first >= max(K, 6) and last >= first
        super().__init__(path,first-2,16);self.first=first;self.last=last;self.K=K
        ww,_,_,_=coefficient_fields(path);b=R(fmpq_poly([0,1]));de=4*b*b/(1-b*b)
        zz=[w/de**k for k,w in enumerate(ww)];self.zfields=[(pair(v),pair(diff(v)),pair(diff(diff(v))))for v in zz]
        self.spins=[]
        for ell in range(first,last+1,2):
            T=ell*(ell+1);B=[F(1)]
            for k in range(1,K+1):B.append(B[-1]*F(T-(k-1)*k,k*k))
            self.spins.append((ell,[exact(x)for x in B]))
    def density(self,lo,hi):
        mid=(lo+hi)/2;h=(hi-lo)/2;beta=exact(lo).union(exact(hi));z=beta*beta;e=1-z;bm=exact(mid);dx=arb(0,exact(h))
        fields=[(val(p,bm),val(dp,bm),val(ddp,beta))for p,dp,ddp in self.zfields]
        if any(not v.is_finite()for r in fields for v in r):return arb.pos_inf()
        sa=4/(1-lo*lo);sb=4/(1-hi*hi)if hi<1 else None;strips=[]
        for (A,B),terms in self.groups.items():
            if sa>=B or(sb is not None and sb<=A):continue
            assert hi<1
            sc=(4/e).intersection(exact(A).union(exact(B)));en=(2*(sc-exact(A))/exact(B-A)-1).intersection(arb(0,1));taper=(1-square(en)).abs_upper()
            l2=sum((square(tight(self.energy[A,B,j],en)).upper()/(2*j+1)for j in {j for j,_,_ in terms}),arb(0))
            dd=sc-4;T=exact(F(4,3));fac=2*dd/(arb.pi()*(T*(dd-T)).sqrt())
            pref=arb.pi()/exact(B-A)*taper*l2.sqrt()*beta*(1+z)**5/(16*e**4)
            strips.append((pref.abs_upper(),fac.upper()))
        # Positive Taylor partial sum for P_l; no interval cancellation.
        powers=[(e/3)**k*z**(self.K-k)for k in range(self.K+1)];numerator=z**(self.K-1);ss=arb(0)
        for ell,B in self.spins:
            f0=sum((B[k]*fields[k][0]for k in range(4)),arb(0));f1=sum((B[k]*fields[k][1]for k in range(4)),arb(0));f2=sum((B[k]*fields[k][2]for k in range(4)),arb(0))
            value=(f0+f1*dx+f2*dx*dx/2).abs_upper()
            for pref,fac in strips:value+=pref*min(arb(1),(fac/exact(F(2*ell+1,2))).sqrt().upper(),key=float)
            den=sum((B[k]*powers[k]for k in range(self.K+1)),arb(0))
            ss+=(2*ell+1)*value*value*numerator/den
        return arb(65536)/(27*arb.pi())*e**3*(2+z)**3/(1+z)**10*ss
def main():
    p=argparse.ArgumentParser();p.add_argument('--witness',type=Path,required=True);p.add_argument('--first',type=int,default=42);p.add_argument('--last',type=int,default=400);p.add_argument('--cells',type=int,default=2048);p.add_argument('--bits',type=int,default=256);p.add_argument('--out',type=Path,required=True);a=p.parse_args();ctx.prec=a.bits;tic=time.time();model=Finite(a.witness,a.first,a.last);print('model ready',time.time()-tic,flush=True)
    heap=[];total=F(0);idx=0
    def insert(lo,hi):
        nonlocal idx,total
        v=model.density(lo,hi)
        if not v.is_finite():
            assert hi-lo>F(1,2**24);m=(lo+hi)/2;insert(lo,m);insert(m,hi);return
        mass=max(F(0),endpoint(v))*(hi-lo);total+=mass;heapq.heappush(heap,(-float(mass),idx,lo,hi,mass));idx+=1
    for i in range(64):insert(F(i,64),F(i+1,64))
    print('initial',len(heap),float(total),time.time()-tic,flush=True)
    while len(heap)<a.cells:
        _,_,lo,hi,mass=heapq.heappop(heap);total-=mass;m=(lo+hi)/2;insert(lo,m);insert(m,hi)
        if len(heap)%256==0:print('cells',len(heap),'Aupper',float(total),'secs',time.time()-tic,flush=True)
    out=dict(status='COMPUTED_FINITE_REAL_NORM_PENDING_INDEPENDENT_REVIEW',witness_sha256=hashlib.sha256(a.witness.read_bytes()).hexdigest(),first=a.first,last=a.last,K=model.K,A_upper_exact=str(total),A_upper_float=float(total),precision_bits=a.bits,seconds=time.time()-tic,cells=[dict(lo=str(r[2]),hi=str(r[3]),mass_upper=str(r[4]))for r in sorted(heap,key=lambda x:x[2])],scope='Finite spin real norm only; requires separate V sign and head proof.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='cells'},indent=2))
if __name__=='__main__':main()
