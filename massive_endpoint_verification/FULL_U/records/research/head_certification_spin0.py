"""Validated spin-zero disk-support integral for the fixed 150-row witness.

Source integrals are exact rational polynomial moments with explicit geometric
remainders. PV terms are exact polynomial/log expressions. Outer integration
uses convexity of the disk support and a validated linear interpolation error
from Arb second derivatives. Tiny logarithmic endpoint neighborhoods are bounded
by analytic log integrals. This certifies only the spin-zero head component.
"""
from __future__ import annotations
import argparse,hashlib,heapq,json,math,pickle,re,sys,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sy

ROOT=Path(__file__).resolve().parent.parent
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import arb,arb_poly,arb_series,ctx
from certify_real_tail import exact,square,power,RealTail

C=F(4,3)

def poly_eval(cs,x):
    """Exact polynomial enclosure via a Taylor shift at the ball midpoint."""
    if not isinstance(cs,arb_poly):cs=arb_poly([exact(q) for q in cs])
    return cs(arb_poly([x.mid(),arb(1)]))(arb(0,x.rad()))

def series_poly(poly,x):
    """Evaluate a polynomial and two derivatives without monomial cancellation."""
    mid=x[0].mid();shift=poly(arb_poly([mid,arb(1)]));delta=x-mid
    result=arb_series([0],prec=3)
    for value in reversed(shift.coeffs()):result=result*delta+value
    return result

def cp_int(i,k):
    """Integral_0^1 x^i (x-1/2)^k dx, exact and cancellation-free."""
    return sum((F(math.comb(i,j),k+j+1) for j in range(i+1)
                if (k+j)%2==0),F(0))/2**(k+i)

def source_cache(path,N=64):
    st=np.load(path);x,y=sy.symbols('x y');groups={};norms={}
    for name,val in zip(st['names'],st['multipliers']):
        name=str(name);v=F(float(val))
        m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
        if m:
            B=F(m[1]);j,n=int(m[2]),int(m[3]);key=('roy',F(4),B)
            term=2*sy.Rational(v)*sy.legendre(n,2*x-1)*sy.legendre(j,y)
            norms[key]=norms.get(key,F(0))+2*abs(v)
        else:
            m=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
            if not m:continue
            A,B=F(m[1]),F(m[2]);j,n=int(m[4]),int(m[5]);key=('strip',A,B)
            term=sy.Rational(v)*(1-(2*x-1)**2)*sy.legendre(n,2*x-1)*sy.legendre(j,2*y-1)
            norms[key]=norms.get(key,F(0))+abs(v)
        groups[key]=groups.get(key,0)+term
    result=[]
    for (kind,A,B),expr in groups.items():
        width=B-A;terms=[(i,j,F(q)) for (i,j),q in sy.Poly(sy.expand(expr),x,y).terms()]
        degree=max(i for i,j,q in terms)
        E=[F(0)]*(degree+1)
        for i,j,q in terms:E[i]+=q/F(j+1)
        while len(E)>1 and E[-1]==0:E.pop()
        R=[sum((E[i]/F(i-k) for i in range(k+1,len(E))),F(0)) for k in range(len(E)-1)]
        if not R:R=[F(0)]
        mean=sum((q/F(i+1) for i,q in enumerate(E)),F(0))
        centered_E=[width**k*sum((q*cp_int(i,k) for i,q in enumerate(E)),F(0)) for k in range(N)]
        order=N+2
        # Raw right-channel source-variable moments, needed by the UV expansion.
        raw_s=[sum((F(math.comb(k,r))*A**(k-r)*width**r
                    *sum((q/F(i+r+1) for i,q in enumerate(E)),F(0))
                    for r in range(k+1)),F(0)) for k in range(order)]
        if kind=='roy':
            h=width/2
            sr=[]
            for r in range(order):
                val=F(0)
                for i,j,q in terms:
                    minus=F(math.factorial(j)*math.factorial(r),math.factorial(j+r+1))
                    plus=sum((F(math.comb(j,k)*(-1)**(j-k)*(2**(r+k+1)-1),r+k+1)
                              for k in range(j+1)),F(0))
                    val+=q/F(i+r+1)*(minus+plus)
                sr.append(val)
            crossed=[sum((F(math.comb(k,r))*(-h)**(k-r)*h**r*sr[r]
                          for r in range(k+1)),F(0)) for k in range(N)]
            geometric=[dict(center=h,radius=h,moments=crossed,norm=2*norms[(kind,A,B)])]
            raw_cross=[(-h)**k*sr[k] for k in range(order)]
        else:
            h=C/2
            st_raw=[C**r*sum((q/F(i+1)*F(math.factorial(j)*math.factorial(r),math.factorial(j+r+1))
                               for i,j,q in terms),F(0)) for r in range(order)]
            mt=[sum((F(math.comb(k,r))*(-h)**(k-r)*st_raw[r] for r in range(k+1)),F(0))
                for k in range(N)]
            hu=(A+B-8-C)/2;ru=(width+C)/2
            cu=[sum((q*sum((F(math.comb(k,r))*width**r*C**(k-r)*cp_int(i,r)*cp_int(j,k-r)
                            for r in range(k+1)),F(0)) for i,j,q in terms),F(0))
                for k in range(order)]
            geometric=[dict(center=h,radius=h,moments=mt,norm=norms[(kind,A,B)]),
                       dict(center=hu,radius=ru,moments=cu[:N],norm=norms[(kind,A,B)])]
            u_raw=[sum((F(math.comb(k,r))*hu**(k-r)*cu[r] for r in range(k+1)),F(0))
                   for k in range(order)]
            raw_cross=[(-1)**k*(st_raw[k]+u_raw[k]) for k in range(order)]
        raw_all=[a+b for a,b in zip(raw_s,raw_cross)]
        uv=[raw_all[k+2]-2*C*raw_all[k+1]+C*C*raw_all[k] for k in range(N)]
        result.append(dict(kind=kind,A=A,B=B,width=width,E=E,R=R,mean=mean,
                           norm=norms[(kind,A,B)],pv_moments=centered_E,
                           geometric=geometric,uv_moments=uv,
                           uv_norm=3*(B+C)**2*norms[(kind,A,B)]))
    return dict(N=N,sources=result,alpha=F(float(st['direction'][0]))-sum((s['mean'] for s in result),F(0)))


def support_upper(U,V):
    u=U.abs_upper();v=V.upper()
    rad=(u*u+v*v).sqrt()
    return (u*u/(rad-v) if v<0 else v+rad).upper()

def support_lower(U,V):
    u=arb(0) if U.contains(0) else U.abs_lower();v=V.lower()
    rad=(u*u+v*v).sqrt()
    return (u*u/(rad-v) if v<0 else v+rad).lower()


class SpinZero:
    def __init__(self,path,N=64,uv_start=96):
        self.path=Path(path);self.sha=hashlib.sha256(self.path.read_bytes()).hexdigest()
        self.N=N;self.uv_start=F(uv_start)
        cache=ROOT/f'output/certification_head/spin0_source_N{N}.pkl'
        builder=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        saved=pickle.loads(cache.read_bytes()) if cache.exists() else {}
        if saved.get('witness')!=self.sha or saved.get('builder')!=builder:
            saved=dict(witness=self.sha,builder=builder,data=source_cache(path,N))
            cache.parent.mkdir(parents=True,exist_ok=True);cache.write_bytes(pickle.dumps(saved))
        self.data=saved['data'];self.sources=self.data['sources'];self.alpha=exact(self.data['alpha'])
        self.real=RealTail(path,40)
        # Form the cosine harmonic polynomial exactly, then evaluate by shifting.
        st=np.load(path);x=sy.Symbol('x');g=sy.Rational(self.data['alpha'])*x
        h=sy.Rational(self.data['alpha'])
        if not np.array_equal(st['direction'],np.array([1.,0.])):
            raise ValueError('This fixed component implementation requires direction (1,0).')
        self.rats=[]
        for name,value in zip(st['names'],st['multipliers']):
            if not re.fullmatch(r'H\d+|Rc_[^_]+_[^_]+|X4|Roy_B[^_]+_J\d+_n\d+|Taper1_Strip_A[^_]+_B[^_]+_T1\.33333_j\d+_n\d+',str(name)):
                raise ValueError(f'Unsupported witness row: {name}')
            m=re.fullmatch(r'H(\d+)',str(name))
            if m:
                j=int(m[1]);g+=sy.Rational(float(value))*sy.chebyshevt(j,x)
                h+=sy.Rational(float(value))*sy.chebyshevu(j-1,x)
            m=re.fullmatch(r'Rc_([^_]+)_([^_]+)',str(name))
            if m:self.rats.append((exact(F(m[1])),exact(F(m[2])),exact(float(value))))
        pg=sy.Poly(sy.expand(g),x)
        self.G=arb_poly([exact(F(pg.nth(i))) for i in range(pg.degree()+1)])
        ph=sy.Poly(sy.expand(h),x)
        self.H=arb_poly([exact(F(ph.nth(i))) for i in range(ph.degree()+1)])
        self.endpoint_cache={}
        for source in self.sources:
            source['Epoly']=arb_poly([exact(q) for q in source['E']])
            source['Rpoly']=arb_poly([exact(q) for q in source['R']])
            source['UVpoly']=arb_poly([exact(q) for q in source['uv_moments']])
            source['PVpoly']=arb_poly([exact(q) for q in source['pv_moments']])
            for term in source['geometric']:
                term['poly']=arb_poly([exact((-1)**k*q) for k,q in enumerate(term['moments'])])

    def forward(self,beta):
        b2=beta*beta;e=1-b2;xmu=e/4
        def arg(b):return (1-6*b*b+power(b,4))/(1+b*b)**2
        x=arg(beta.lower()).union(arg(beta.upper())).intersection(arb(0,1))
        rg=poly_eval(self.G,x)
        for r,c,v in self.rats:
            y2=1-x*x;dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x
            nr=1-r*c*x;ni=-r*c;mminus=(1+r*c)/(1+2*r*c+r*r)
            den=square(2*r*x-c*(1+r*r))+(1-c*c)*square(1-r*r)
            mr=(nr*dr+y2*ni*di)/den;mi=(ni*dr-nr*di)/den
            rg+=v*((2*x*x-1)*(mr-mminus)-2*x*y2*mi)
        U=16*beta*e*e/(1+b2)**5*self.real.K(beta)
        V=4*e/(1+b2)**3*rg-self.alpha*xmu/(2*(1-2*xmu)*(3-4*xmu))
        return U,V

    def crossed(self,source,mu):
        result=arb(0)
        for term in source['geometric']:
            den=mu+exact(term['center']);rad=exact(term['radius'])
            result+=term['poly'](1/den)/den
            err=exact(term['norm'])*power(rad/den,self.N)/(den-rad)
            result+=arb(0,err.upper())
        return result

    def pv(self,source,mu,omit_log=False):
        A,B,width=source['A'],source['B'],source['width']
        den=mu-exact((A+B)/2);rad=exact(width/2)
        if not omit_log and den.abs_lower()>rad*exact(F(3,2)):
            val=source['PVpoly'](1/den)/den
            err=exact(source['norm'])*power(rad/abs(den),self.N)/(abs(den)-rad)
            return val+arb(0,err.upper()),None
        z=(mu-exact(A))/exact(width)
        E=poly_eval(source['Epoly'],z);R=poly_eval(source['Rpoly'],z)
        if omit_log:return -R/exact(width),E/exact(width)
        logratio=abs(mu-exact(A)).log()-abs(mu-exact(B)).log()
        return (E*logratio-R)/exact(width),None

    def fields(self,a,b,omit_sources=()):
        if a==b==1:return arb(0),arb(0),[]
        beta=arb(exact((a+b)/2),exact((b-a)/2));b2=beta*beta;e=1-b2
        U,V=self.forward(beta);logs=[]
        # Use a factored source series on the UV endpoint, including beta=1.
        if 4/(1-exact(a)*exact(a))>=exact(self.uv_start):
            x=e/4
            for source in self.sources:
                val=source['UVpoly'](x)
                z=exact(source['B'])*x
                err=exact(source['uv_norm'])*power(z,self.N)/(1-z)
                val+=arb(0,err.upper())
                V-=x/(16*(1-exact(C)*x)**2)*val
            return U,V,logs
        mu=4/e;Q=mu-exact(C)
        for si,source in enumerate(self.sources):
            A,B,width=source['A'],source['B'],source['width']
            if not (mu.upper()<=exact(A) or mu.lower()>=exact(B)):
                z=(mu-exact(A))/exact(width)
                val=arb.pi()/exact(width)*poly_eval(source['Epoly'],z)/(e*e)
                if not (mu.lower()>exact(A) and mu.upper()<exact(B)):val=val.union(arb(0))
                U+=val
            pv,logcoef=self.pv(source,mu,si in omit_sources)
            V+=(-pv-self.crossed(source,mu)+3*exact(source['mean'])/Q)/(e*e)
            if logcoef is not None:logs.append((si,logcoef/(e*e)))
        return U,V,logs

    def height(self,a,b):
        U,V,_=self.fields(a,b)
        hi=4/arb.pi()*support_upper(U,V);lo=4/arb.pi()*support_lower(U,V)
        return hi.upper(),max(arb(0),lo.lower(),key=float)

    def fields_series(self,a,b):
        """Enclose derivatives of an analytic approximant, plus its V error.

        The coefficient of order two encloses half the second derivative at
        every beta in [a,b]. Geometric-series errors are separately bounded,
        never differentiated as though constant.
        """
        ball=arb(exact((a+b)/2),exact((b-a)/2))
        beta=arb_series([ball,1],prec=3);b2=beta*beta;e=1-b2;xmu=e/4
        x=(1-6*b2+b2*b2)/(1+b2)**2
        rg=series_poly(self.G,x);kg=series_poly(self.H,x)
        for r,c,v in self.rats:
            y2=1-x*x;dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x
            nr=1-r*c*x;ni=-r*c;mminus=(1+r*c)/(1+2*r*c+r*r)
            den=(2*r*x-c*(1+r*r))**2+(1-c*c)*(1-r*r)**2
            mr=(nr*dr+y2*ni*di)/den;mi=(ni*dr-nr*di)/den
            rg+=v*((2*x*x-1)*(mr-mminus)-2*x*y2*mi)
            kg+=v*((2*x*x-1)*mi+2*x*(mr-mminus))
        U=16*beta*e*e/(1+b2)**5*kg
        V=4*e/(1+b2)**3*rg-self.alpha*xmu/(2*(1-2*xmu)*(3-4*xmu))
        error=arb(0)
        if 4/(1-exact(a)*exact(a))>=exact(self.uv_start):
            fac=xmu/(16*(1-exact(C)*xmu)**2)
            for source in self.sources:
                V-=fac*series_poly(source['UVpoly'],xmu)
                z=exact(source['B'])*xmu[0]
                err=exact(source['uv_norm'])*power(z,self.N)/(1-z)
                error+=fac[0].abs_upper()*err.abs_upper()
            return U,V,error.upper()
        mu=4/e;Q=mu-exact(C)
        for source in self.sources:
            A,B,width=source['A'],source['B'],source['width']
            z=(mu-exact(A))/exact(width)
            if 4/(1-a*a)>A and 4/(1-b*b)<B:
                U+=arb.pi()/exact(width)*series_poly(source['Epoly'],z)/(e*e)
            den=mu-exact((A+B)/2);rad=exact(width/2)
            if den[0].abs_lower()>rad*exact(F(3,2)):
                pv=series_poly(source['PVpoly'],1/den)/den
                err=exact(source['norm'])*power(rad/abs(den[0]),self.N)/(abs(den[0])-rad)
                error+=err.abs_upper()/(e[0]*e[0]).abs_lower()
            else:
                da=mu-exact(A);db=mu-exact(B)
                if da[0]<0:da=-da
                if db[0]<0:db=-db
                if not (da[0]>0 and db[0]>0):
                    return U,arb_series([arb('nan')],prec=3),arb('nan')
                pv=(series_poly(source['Epoly'],z)*(da.log()-db.log())
                    -series_poly(source['Rpoly'],z))/exact(width)
            cross=arb_series([0],prec=3)
            for term in source['geometric']:
                den=mu+exact(term['center']);rad=exact(term['radius'])
                cross+=series_poly(term['poly'],1/den)/den
                err=exact(term['norm'])*power(rad/den[0],self.N)/(den[0]-rad)
                error+=err.abs_upper()/(e[0]*e[0]).abs_lower()
            V+=(-pv-cross+3*exact(source['mean'])/Q)/(e*e)
        return U,V,error.upper()

    def trapezoid(self,a,b):
        """Convex disk support and rigorous linear interpolation remainder."""
        try:U,V,err=self.fields_series(a,b)
        except ValueError:return arb('nan'),arb(0)
        if not (U[2].is_finite() and V[2].is_finite() and err.is_finite()):
            return arb('nan'),arb(0)
        for x in (a,b):
            if x not in self.endpoint_cache:
                u,v,_=self.fields(x,x)
                self.endpoint_cache[x]=support_upper(u,v)
        width=exact(b-a)
        # |F-linear(F)| <= h² sup|F''|/8. The disk support obeys
        # H(F+dF)<=H(F)+|dU|+2|dV|. Actual V differs from its
        # analytic approximant by err; interpolation adds at most 2err.
        curvature=2*U[2].abs_upper()+4*V[2].abs_upper()
        rem=width**3/8*curvature+4*width*err
        up=4/arb.pi()*(width/2*(self.endpoint_cache[a]+self.endpoint_cache[b])+rem)
        # A certified lower bound remains useful for diagnostics, while a
        # midpoint estimate gives an efficient (noncertifying) refinement key.
        mid=(a+b)/2;um,vm,_=self.fields(mid,mid)
        estimate=4/arb.pi()*width*support_lower(um,vm)
        return up.upper(),estimate.lower()

    def gaps(self,bits=26):
        unit=F(1,2**bits);out=[(F(0),unit)]
        for endpoint in sorted(set(s[k] for s in self.sources for k in ('A','B'))-{F(4)}):
            root=(1-exact(4)/exact(endpoint)).sqrt()
            q=math.floor(float(root)*2**bits)
            a,b=F(q-1,2**bits),F(q+2,2**bits)
            assert exact(a)<root and root<exact(b)
            out.append((a,b))
        return sorted(out)

    def log_mass(self,a,b,root):
        width=exact(b-a);aa,bb=exact(a),exact(b)
        if not (root<aa or root>bb):
            # Concavity in the singularity position gives the midpoint maximum.
            return width*(1-(width/2).log())
        dist=abs(arb(exact((a+b)/2),exact((b-a)/2))-root)
        return width*abs(dist.log()).upper()

    def gap_bound(self,a,b):
        omit=[]
        for i,s in enumerate(self.sources):
            roots=[(1-exact(4)/exact(s[k])).sqrt() for k in ('A','B')]
            if any(not (r<exact(a) or r>exact(b)) for r in roots):omit.append(i)
        U,V,logs=self.fields(a,b,omit)
        total=exact(b-a)*(U.abs_upper()+2*V.abs_upper())
        for si,coef in logs:
            s=self.sources[si];A,B=s['A'],s['B']
            ra=(1-exact(4)/exact(A)).sqrt();rb=(1-exact(4)/exact(B)).sqrt()
            lm=exact(b-a)*abs((exact(A)/exact(B)).log()).upper()
            for r in (ra,rb,-ra,-rb):lm+=self.log_mass(a,b,r)
            total+=2*coef.abs_upper()*lm
        return (4/arb.pi()*total).upper()


def integrate(model,target,max_cells=30000,gap_bits=26):
    gaps=model.gaps(gap_bits);gapmass=sum((model.gap_bound(a,b) for a,b in gaps),arb(0))
    heap=[];serial=0;upper=gapmass;lower=arb(0);tic=time.time()
    # Intervals use exact rational beta endpoints; a switch point above UV start
    # avoids evaluating mu at beta=1 in the finite-energy expression.
    cut=F(math.ceil(math.sqrt(1-4/float(model.uv_start))*2**20),2**20)
    points=set([F(i,256) for i in range(257)]+[cut]+[v for g in gaps for v in g])
    points=sorted(points)
    def insert(a,b):
        nonlocal serial,upper,lower
        up,low=model.trapezoid(a,b)
        if not up.is_finite():
            if b-a<F(1,2**45):raise RuntimeError(f'Nonfinite cell {a} {b}')
            c=(a+b)/2;insert(a,c);insert(c,b);return
        upper+=up;lower+=low
        heapq.heappush(heap,(-float(up-low),serial,a,b,up,low));serial+=1
    for a,b in zip(points,points[1:]):
        if any(a>=ga and b<=gb for ga,gb in gaps):continue
        insert(a,b)
    print(json.dumps({'event':'initial','cells':len(heap),'upper':str(upper.upper()),'midpoint_estimate':str(lower),'gaps':str(gapmass)}),flush=True)
    while len(heap)<max_cells and not upper<arb(target):
        _,_,a,b,up,low=heapq.heappop(heap);upper-=up;lower-=low
        c=(a+b)/2;insert(a,c);insert(c,b)
        if len(heap)%1000==0:
            print(json.dumps({'event':'refine','cells':len(heap),'upper':float(upper.upper()),'midpoint_estimate':float(lower), 'seconds':time.time()-tic}),flush=True)
    return dict(upper=str(upper.upper()),midpoint_estimate=str(lower),target=str(target),
                target_met=bool(upper<arb(target)),cells=len(heap),gap_upper=str(gapmass.upper()),
                gap_beta_intervals=[[str(a),str(b)] for a,b in gaps],seconds=time.time()-tic,
                upper_dyadic=list(map(str,upper.upper().man_exp())),
                cells_data=[{'a':str(a),'b':str(b),'upper_dyadic':list(map(str,up.upper().man_exp()))}
                            for _,_,a,b,up,low in sorted(heap,key=lambda q:q[2])])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,default=ROOT/'output/angles_resummed/angle_000p000.npz')
    ap.add_argument('--target',default='2.695');ap.add_argument('--max-cells',type=int,default=30000)
    ap.add_argument('--precision',type=int,default=192);ap.add_argument('--N',type=int,default=64)
    ap.add_argument('--gap-bits',type=int,default=26);ap.add_argument('--out',type=Path,default=ROOT/'output/certification_head/spin0_component.json')
    a=ap.parse_args();ctx.prec=a.precision;tic=time.time();model=SpinZero(a.witness,a.N)
    print(json.dumps({'event':'built','seconds':time.time()-tic,'sources':len(model.sources),'alpha':str(model.data['alpha'])}),flush=True)
    r=integrate(model,a.target,a.max_cells,a.gap_bits)
    r.update(witness_sha256=model.sha,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),N=a.N,
             real_tail_source_sha256=hashlib.sha256((ROOT/'research/certify_real_tail.py').read_bytes()).hexdigest(),
             arithmetic=f'Arb{a.precision}, exact rational source moments with analytic geometric remainders',
             method='Convex support interpolation with interval second derivatives; analytic logarithmic endpoint caps',
             claim='Spin-zero dual disk-support integral upper bound only; not full exclusion')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2))
    print(json.dumps({k:v for k,v in r.items() if k!='cells_data'},indent=2))

if __name__=='__main__':main()
