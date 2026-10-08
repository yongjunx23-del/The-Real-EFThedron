"""Directed Arb certificate for the real-kernel A-tail (not a full dual bound).

Supports forward harmonics, rational rows, crossing, finite-J Roy and power-1
tapered physical strips. Computes A=sum_even_l>L int U_l^2/p_l dnu_l.
All infinite spins and the unbounded spectral interval are included using an
analytic positive-coefficient Legendre inequality and the beta compactification.
The remaining V-tail condition and finite-head support are separate obligations.
"""
from __future__ import annotations
import argparse, hashlib, heapq, json, math, re, sys, time
from pathlib import Path
from fractions import Fraction
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'tmp/certification_deps'))
from flint import arb, fmpq, ctx


def exact(x):
    if isinstance(x,Fraction):return arb(fmpq(x.numerator,x.denominator))
    if isinstance(x,(float,np.floating)):
        a,b=float(x).as_integer_ratio();return arb(fmpq(a,b))
    return arb(x)


def dyadic(x):
    return [str(t) for t in x.man_exp()]


def power(x,n):
    """Integer power by multiplication, also for balls containing zero."""
    if n<0:return 1/power(x,-n)
    result=arb(1)
    while n:
        if n%2:result*=x
        x=x*x;n//=2
    return result


def square(x):
    lo,hi=x.lower(),x.upper()
    if x.contains(0):return arb(0).union(x.abs_upper()*x.abs_upper())
    return (lo*lo).union(hi*hi)


class RealTail:
    def __init__(self,path,L,max_k=30,combined_strips=False):
        self.path=Path(path);self.L=L
        self.combined_strips=combined_strips
        d=np.load(path);self.names=list(map(str,d['names']))
        self.lam=d['multipliers'];self.c0=exact(float(d['direction'][0]))
        self.h=[];self.r=[];self.groups={}
        for name,value in zip(self.names,self.lam):
            v=exact(value)
            m=re.fullmatch(r'H(\d+)',name)
            if m:self.h.append((int(m[1]),v));continue
            m=re.fullmatch(r'Rc_([^_]+)_([^_]+)',name)
            if m:self.r.append((exact(Fraction(m[1])),exact(Fraction(m[2])),'cos',v));continue
            m=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',name)
            if m:self.r.append((exact(Fraction(m[2])),exact(Fraction(m[3])),m[1],v));continue
            if re.fullmatch(r'X\d+(?:_\d+)?',name):continue
            m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
            if m:
                J,n=int(m[2]),int(m[3])
                if J>L:raise ValueError('Include all local Roy spins in head')
                if J==n==0:self.c0-=2*v
                continue
            m=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
            if m:
                A,B,T=map(float,m.group(1,2,3));j,n=int(m[4]),int(m[5])
                if not 0<T<=A-4:raise ValueError('Nonphysical strip angle')
                self.groups.setdefault((A,B),[]).append((j,n,v))
                if j==0 and n==0:self.c0-=v*exact(Fraction(2,3))
                if j==0 and n==2:self.c0+=v*exact(Fraction(2,15))
                continue
            raise ValueError('Unsupported row '+name)
        windows=sorted(self.groups)
        if any(a[1]>b[0] for a,b in zip(windows,windows[1:])):
            raise ValueError('Separate strip L2 sum requires disjoint source windows')
        self.ks=range(2,min(max_k,L+2)+1)
        # B_k(mu)=coef[k]*(mu-4)^k. Exact rational spin-sum coefficients.
        self.coef={k:Fraction(math.factorial(k)**2,2*(k-1)*math.prod(range(L-k+3,L+k+1)))*Fraction(3,4)**k for k in self.ks}
        self.taylor_order=6
        p=self.taylor_order
        self.poly_remainder=sum((abs(v)*exact(Fraction(2**p*math.factorial(n+p),math.factorial(2*p+1)*math.factorial(n-1-p))) for n,v in self.h if n>p),arb(0))

    def K(self,beta):
        b2=beta*beta
        # cos(4 atan beta) decreases on [0,1]; use endpoint images to avoid
        # dependency inflation near the quadratic UV endpoint.
        def argument(t):
            t2=t*t;return (1-6*t2+t2*t2)/(1+t2)**2
        x=argument(beta.lower()).union(argument(beta.upper())).intersection(arb(0,1))
        if beta.contains(0):x=x.union(arb(1))
        if beta.contains(1):x=x.union(arb(-1))
        mid=x.mid();h=arb(0,x.rad());p=self.taylor_order
        vals=[[arb(1)]+[arb(0)]*(p-1),[2*mid,arb(2)]+[arb(0)]*(p-2)]
        for j in range(2,max([n for n,_ in self.h]+[2])):
            vals.append([2*mid*vals[-1][k]+(2*vals[-1][k-1] if k else 0)-vals[-2][k] for k in range(p)])
        coeff=[sum((v*vals[n-1][k] for n,v in self.h),arb(0)) for k in range(p)]
        result=coeff[-1]
        for value in reversed(coeff[:-1]):result=result*h+value
        result+=self.c0+arb(0,(self.poly_remainder*power(x.rad(),p)).upper())
        for r,c,kind,v in self.r:
            y2=1-x*x;dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x
            if kind=='cos':nr=1-r*c*x;ni=-r*c;mminus=(1+r*c)/(1+2*r*c+r*r)
            else:nr=r*x;ni=r;mminus=-r/(1+2*r*c+r*r)
            # Positive completed-square denominator avoids cancellation.
            den=square(2*r*x-c*(1+r*r))+(1-c*c)*square(1-r*r)
            mr=(nr*dr+y2*ni*di)/den
            mi=(ni*dr-nr*di)/den
            result+=v*((2*x*x-1)*mi+2*x*(mr-mminus))
        return result

    def density_forward(self,a,b):
        beta=arb(exact((a+b)/2),exact((b-a)/2))
        b2=beta*beta
        common=exact(65536)/(27*arb.pi())*b2*(2+b2)**3/(1+b2)**10*square(self.K(beta))
        candidates=[]
        for k in self.ks:
            if b==1 and k>3:continue
            value=common*exact(self.coef[k])*4**k*power(b2,k)*power(1-b2,3-k)
            if value.is_finite():candidates.append((value.upper(),k))
        if not candidates:return (arb.pos_inf(),2),0
        return min(candidates,key=lambda x:float(x[0])),len(candidates)

    def density(self,a,b):
        # A cell overlapping a source window must include that strip's bound.
        beta=arb(exact((a+b)/2),exact((b-a)/2));b2=beta*beta
        sa=exact(4)/(1-exact(a)*exact(a))
        sb=exact(4)/(1-exact(b)*exact(b)) if b<1 else None
        active=[(A,B,terms) for (A,B),terms in self.groups.items()
                if not (sa>=exact(B) or (sb is not None and sb<=exact(A)))]
        if not self.combined_strips or not active:return self.density_forward(a,b)
        if b==1:return (arb.pos_inf(),2),0
        mu=4/(1-b2);delta=mu-4;Q=mu-exact(Fraction(4,3))
        ubound=abs(16*beta*power(1-b2,4)/(1+b2)**5*self.K(beta))
        for A,B,terms in active:
            sclip=mu.intersection(exact(A).union(exact(B)))
            x=(2*(sclip-exact(A))/exact(B-A)-1).intersection(arb(0,1))
            ps=[arb(1),x]
            for n in range(1,max(n for _,n,_ in terms)):
                ps.append(((2*n+1)*x*ps[-1]-n*ps[-2])/(n+1))
            l2=arb(0)
            for j in sorted({j for j,_,_ in terms}):
                poly=sum((v*ps[n] for jj,n,v in terms if jj==j),arb(0))
                l2+=square(poly).upper()/(2*j+1)
            taper=(1-square(x)).abs_upper()
            dd=sclip-4;T=exact(Fraction(4,3))
            # Bernstein + Cauchy-Schwarz on v∈[0,1]:
            # mean|w(v)P_l(1-2Tv/dd)| <= ||w||_2 * sqrt(2dd/[pi(l+.5)sqrt(T(dd-T))]).
            factor=2*dd/(arb.pi()*exact(Fraction(2*self.L+5,2))*(T*(dd-T)).sqrt())
            factor_upper=min(arb(1),factor.upper(),key=float)
            ubound+=arb.pi()/exact(B-A)*taper*(l2*factor_upper).sqrt()
        common=square(ubound)*Q**3/(2*arb.pi())*8/power(1-b2,2)
        candidates=[]
        for k in self.ks:
            value=common*exact(self.coef[k])*power(delta,k)
            if value.is_finite():candidates.append((value.upper(),k))
        if not candidates:return (arb.pos_inf(),2),0
        return min(candidates,key=lambda x:float(x[0])),len(candidates)

    def forward_bound(self,target,max_cells=20000):
        heap=[];serial=0;total=arb(0)
        def insert(a,b):
            nonlocal serial,total
            pair,_=self.density(a,b);height,k=pair
            if not height.is_finite():
                if b-a<Fraction(1,2**30):raise ValueError('Unresolved finite enclosure')
                c=(a+b)/2;insert(a,c);insert(c,b);return
            mass=exact(b-a)*height
            total+=mass
            heapq.heappush(heap,(-float(mass),serial,a,b,mass,k));serial+=1
        # Small initial intervals avoid nonsensical broad rational enclosures.
        for i in range(256):insert(Fraction(i,256),Fraction(i+1,256))
        target=arb(target)
        while len(heap)<max_cells and not total<target:
            _,_,a,b,mass,_=heapq.heappop(heap);total-=mass
            c=(a+b)/2;insert(a,c);insert(c,b)
        cells=[{'a':str(a),'b':str(b),'mass_upper':dyadic(mass.upper()),'k':k} for _,_,a,b,mass,k in sorted(heap,key=lambda x:x[2])]
        return total.upper(),cells

    def strip_bound(self):
        total=arb(0);rows=[]
        for (A,B),terms in self.groups.items():
            aa,bb=exact(A),exact(B)
            u=arb.pi()/(bb-aa)*sum((abs(v) for _,_,v in terms),arb(0))
            inverse,k=min(((exact(self.coef[k])*(bb-4)**k,k) for k in self.ks),key=lambda z:float(z[0]))
            value=u*u*(bb-exact(Fraction(4,3)))**3/(2*arb.pi())*inverse*(bb-aa)/(1-4/aa).sqrt()
            total+=value
            rows.append({'window':[A,B],'k':k,'A_upper':value.upper().str(25),'U_uniform_upper':u.upper().str(25)})
        return total.upper(),rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('witness',type=Path)
    ap.add_argument('--L',type=int,default=100);ap.add_argument('--target',default='5e-8')
    ap.add_argument('--max-cells',type=int,default=20000);ap.add_argument('--precision',type=int,default=160)
    ap.add_argument('--combined-strips',action='store_true',help='Integrate tighter pointwise strip L2/Bernstein bound with the forward part')
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.L%2 or args.L<2:raise ValueError('Even L>=2 required')
    ctx.prec=args.precision;tic=time.time();bank=RealTail(args.witness,args.L,combined_strips=args.combined_strips)
    forward,cells=bank.forward_bound(args.target,args.max_cells)
    strip,strip_rows=bank.strip_bound()
    total=forward if args.combined_strips else (forward.sqrt()+strip.sqrt())**2
    result={'witness':str(args.witness.resolve()),'sha256':hashlib.sha256(args.witness.read_bytes()).hexdigest(),
            'L':args.L,'precision_bits':args.precision,'arithmetic':'python-flint 0.9.0 Arb directed ball arithmetic',
            'forward_A_upper':None if args.combined_strips else forward.str(30),'strip_A_upper_separate_crude':strip.str(30),
            'combined_strip_integration':args.combined_strips,
            'combined_A_upper':total.upper().str(30),'combined_A_upper_dyadic':dyadic(total.upper()),
            'cells':len(cells),'target_forward_upper':args.target,'target_met':bool(forward<arb(args.target)),
            'strip_components':strip_rows,'seconds':time.time()-tic,
            'scope':'Rigorous real-kernel A-tail upper bound for every even spin>L and all mu>=4, with exact-binary64 witness coefficients, exact decimal rational row parameters and T=4/3; physical constraints assumed valid. NOT a full dual exclusion: V tail and finite head remain unproved.',
            'method':'Exact positive Legendre Taylor coefficient, telescoping even-spin sum, interval upper Riemann sums over beta∈[0,1], exact UV cancellation; optional pointwise strip L2/Bernstein bound or separate crude triangle/Minkowski bound.',
            'forward_partition':cells}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='forward_partition'},indent=2))


if __name__=='__main__':main()
