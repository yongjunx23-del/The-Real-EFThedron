"""Complete all-spin real norm for 220-row witnesses, derivative order <=3.

Exact rational coefficient discovery; directed reference integration. This is
only the real norm A, not a full coupling certificate without the V/head proof.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import sys,re,json,argparse,time,hashlib,heapq
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research/all_direction_certification_20261007/alternative'),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from rational_head import R,compose
from exact_bernstein import fmpq_poly,q
from flint import arb,arb_poly,ctx
from certify_real_tail import exact,square
from matrix_real_tail_gram import legendre_reference_parameter
def mul(a,b):return [sum((a[k]*b[j-k]for k in range(j+1)),R(0))for j in range(4)]
def scale(a,c):return [x*c for x in a]
def add(a,b):return [x+y for x,y in zip(a,b)]
def polyjet(pol,x):
    v=[R(0)]*4
    for co in reversed(pol.coeffs()):v=mul(v,x);v[0]+=R(co)
    return v
def cheb(n,second=True):
    x=fmpq_poly([0,1]);p=[fmpq_poly([1]),(2 if second else 1)*x]
    for j in range(2,n+1):p.append(2*x*p[-1]-p[-2])
    return p
def fractions(p):return [F(str(v))for v in p.coeffs()]
def coefficient_fields(path):
    st=np.load(path);rows={str(n):F(float(v))for n,v in zip(st['names'],st['multipliers'])};alpha=F(float(st['direction'][0]));groups={}
    for n,v in rows.items():
        m=re.fullmatch(r'Roy_B[^_]+_J(\d+)_n(\d+)',n)
        if m and int(m[1])==int(m[2])==0:alpha-=2*v
        m=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',n)
        if m:
            A,B=F(m[1]),F(m[2]);j,k=int(m[4]),int(m[5]);groups.setdefault((A,B),[]).append((j,k,v))
            if j==k==0:alpha-=F(2,3)*v
            if j==0 and k==2:alpha+=F(2,15)*v
    beta=R(fmpq_poly([0,1]));z=beta**2;e=1-z;x=(1-6*z+z*z)/(1+z)**2;us=cheb(44)
    forward=fmpq_poly([q(alpha)]);Dp={p:fmpq_poly([])for p in [1,2,3]}
    rats=[]
    for n,v in rows.items():
        m=re.fullmatch(r'H(\d+)',n)
        if m:forward+=q(v)*us[int(m[1])-1]
        m=re.fullmatch(r'Dt([123])(?:H_j(\d+)|G_regular)',n)
        if m:
            p=int(m[1]);j=int(m[2])if m[2]else 1
            Dp[p]+=q(v)*sum(((-1)**h*comb(2*p-2,h)*us[j+h-1]for h in range(2*p-1)),fmpq_poly([]))
        m=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',n)
        if m:rats.append((m[1],F(m[2]),F(m[3]),v))
    K=compose(forward,x)
    for kind,r,c,v in rats:
        y2=1-x*x;dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x
        if kind=='cos':nr=1-r*c*x;ni=-r*c;mm=(1+r*c)/(1+2*r*c+r*r)
        else:nr=r*x;ni=r;mm=-r/(1+2*r*c+r*r)
        den=(2*r*x-c*(1+r*r))**2+(1-c*c)*(1-r*r)**2
        mr=(nr*dr+y2*ni*di)/den;mi=(ni*dr-nr*di)/den
        K+=v*((2*x*x-1)*mi+2*x*(mr-mm))
    q0=2*(1+z)/e
    A=[R(2),R(F(1,2)),R(0),R(0)]
    invq2=[(-1)**k*(k+1)/(2*q0)**k/q0**2 for k in range(4)]
    xx=scale(mul(A,A),2);xx=mul(xx,invq2);xx[0]-=1
    sqrtjet=[];co=F(1)
    for k in range(4):
        if k:co*=F(1,2)-(k-1);co/=k
        sqrtjet.append(co*(e/4)**k)
    Z=mul([R(F(comb(3,k),4**k))for k in range(4)],[(-1)**k*comb(4+k,k)/(2*q0)**k for k in range(4)])
    Z=scale(mul(Z,sqrtjet),z)
    ww=[z*K,R(0),R(0),R(0)]
    for p in [1,2,3]:
        coeff=mul(Z,polyjet(Dp[p],xx))
        for k in range(p+1):ww[k]+=factorial(p)*coeff[p-k]
    return ww,groups,rows,alpha
def endpoint(v,upper=True):return F(str((v.upper()if upper else v.lower()).fmpq()))
def tight(p,x):return p(arb_poly([x.mid(),1]))(arb(0,x.rad()))
class Norm:
    def __init__(self,path,L=40,order=16):
        assert L >= order and L % 2 == 0 and order >= 8
        self.path=Path(path);self.L=L;ww,self.groups,self.rows,self.alpha=coefficient_fields(path)
        self.coefficients=[dict(numerator=list(map(str,fractions(v.n))),denominator=list(map(str,fractions(v.d))))for v in ww]
        self.ww=[(arb_poly([exact(x)for x in fractions(v.n)]),arb_poly([exact(x)for x in fractions(v.d)]))for v in ww]
        self.order=order;self.a,self.refs=legendre_reference_parameter(L,order);self.T=L*(L+1);self.ratio=F((L+2)*(L+3),self.T);self.factor=F(2*L+5,4*L+6)
        self.energy={}
        for (A,B),terms in self.groups.items():
            assert A>=12 and B>A
            for j in {j for j,_,_ in terms}:
                ps=[arb_poly([1]),arb_poly([0,1])]
                for n in range(1,max(n for _,n,_ in terms)):ps.append(((2*n+1)*arb_poly([0,1])*ps[-1]-n*ps[-2])*(arb(1)/(n+1)))
                self.energy[A,B,j]=sum((exact(v)*ps[n]for jj,n,v in terms if jj==j),arb_poly([]))
    def density(self,lo,hi):
        beta=exact(lo).union(exact(hi));z=beta*beta;e=1-z
        w=[(tight(n,beta)/tight(d,beta)).abs_upper()for n,d in self.ww]
        if any(not v.is_finite()for v in w):return arb.pos_inf()
        sa=4/(1-lo*lo);sb=4/(1-hi*hi)if hi<1 else None
        for (A,B),terms in self.groups.items():
            if sa>=B or(sb is not None and sb<=A):continue
            assert hi<1
            mu=4/e;sc=mu.intersection(exact(A).union(exact(B)));en=(2*(sc-exact(A))/exact(B-A)-1).intersection(arb(0,1));taper=(1-square(en)).abs_upper()
            l2=sum((square(tight(self.energy[A,B,j],en)).upper()/(2*j+1)for j in {j for j,_,_ in terms}),arb(0))
            dd=sc-4;T=exact(F(4,3));fac=2*dd/(arb.pi()*exact(F(2*self.L+5,2))*(T*(dd-T)).sqrt());fac=min(arb(1),fac.upper(),key=float)
            strip=arb.pi()/exact(B-A)*taper*(l2*fac).sqrt()
            invUref=beta*(1+z)**5/(16*e**4)
            w[0]+=(strip*invUref).abs_upper()
        def y_at(b):
            z=exact(b)*exact(b);return 4*z/(4*z+exact(self.a*self.T)*(1-z))
        y=y_at(lo).union(y_at(hi));s=arb(0)
        for i in range(4):
            for j in range(4):
                m=i+j
                integ=sum((exact(F((-1)**k*comb(m,k),self.order-1-m+k))*y**(self.order-1-m+k)for k in range(m+1)),arb(0))
                s+=w[i]*w[j]*exact(self.ratio**m/self.a**(m+1)/F(factorial(i)**2*factorial(j)**2))*integ
        common=arb(262144)/(27*arb.pi())*e*e*(2+z)**3/(1+z)**10*exact(self.factor)
        return common*s
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--L',type=int,default=40);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cells',type=int,default=4096);ap.add_argument('--bits',type=int,default=256);ap.add_argument('--probe',action='store_true');a=ap.parse_args();ctx.prec=a.bits;tic=time.time();model=Norm(a.witness,a.L);print('fields built',time.time()-tic,flush=True)
    if a.probe:
        for n in [64,256,1024]:
            vals=[model.density(F(i,n),F(i+1,n))for i in range(n)];print('uniform',n,'upper',str(sum((v.upper()for v in vals),arb(0))/n),'secs',time.time()-tic,flush=True)
        return
    heap=[];total=F(0);ident=0
    def insert(lo,hi):
        nonlocal total,ident
        v=model.density(lo,hi)
        if not v.is_finite():
            assert hi-lo>F(1,2**24);mid=(lo+hi)/2;insert(lo,mid);insert(mid,hi);return
        mass=max(F(0),endpoint(v))*(hi-lo);total+=mass;heapq.heappush(heap,(-float(mass),ident,lo,hi,mass));ident+=1
    for i in range(64):insert(F(i,64),F(i+1,64))
    while len(heap)<a.cells:
        _,_,lo,hi,mass=heapq.heappop(heap);total-=mass;mid=(lo+hi)/2;insert(lo,mid);insert(mid,hi)
        if len(heap)%512==0:print('cells',len(heap),'A upper',float(total),'secs',time.time()-tic,flush=True)
    cells=sorted(heap,key=lambda r:r[2]);out=dict(status='PASS_REAL_NORM_COMPONENT',witness=str(a.witness.resolve()),witness_sha256=hashlib.sha256(a.witness.read_bytes()).hexdigest(),L=a.L,precision_bits=a.bits,A_upper_exact=str(total),A_upper_float=float(total),coefficients=model.coefficients,a=str(model.a),reference_coefficients=list(map(str,model.refs)),ratio=str(model.ratio),factor=str(model.factor),alpha_exact=str(model.alpha),seconds=time.time()-tic,cells=[dict(lo=str(r[2]),hi=str(r[3]),mass_upper=str(r[4]))for r in cells],scope='Complete real norm sum over all even l>L and all mu>4, conditional on analytic positive Legendre/moment and strip Bernstein estimates. Not a complete disk support or coupling bound; independent formula review still required.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['cells','coefficients']},indent=2))
if __name__=='__main__':main()
