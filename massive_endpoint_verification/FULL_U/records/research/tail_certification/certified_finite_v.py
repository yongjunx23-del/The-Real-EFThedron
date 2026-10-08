"""Certified finite-spin V kernel assembled from exact polynomial moments.

No optimization. Moment geometric remainders are enclosed explicitly. Roy PV
residues vanish for the supported spins l>=8. Tapered strip PV is analytic.
"""
import sys,math,pickle,json,argparse,re,time,hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research'),str(ROOT/'tmp/certification_deps'),str(Path(__file__).parent)]
from flint import arb,acb,arb_poly,ctx
from certified_roy_moments import load
from certify_finite_real import lp,tight
from certify_real_tail import RealTail,exact,square
from certify_uv import ar

X=arb_poly([0,1]);C=F(4,3)
EXPECTED='bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c'
def coeff(l,j):return math.factorial(l+j)//(math.factorial(l-j)*math.factorial(j)**2)
def rev(p,n=None):
 n=p.degree() if n is None else n
 return arb_poly([p[n-i] for i in range(n+1)])
def psum(values):return sum(values,arb_poly([]))
def geometric_num(polys,den,N):
 result=arb_poly([])
 for k in range(N):result=result*den+polys[k]
 return result

def xlog(x):
 assert x.abs_upper()<ar(F(1,4))
 def f(a):return arb(0) if a==0 else a*abs(a).log()
 return f(x.upper()).union(f(x.lower()))
def taper_q(n,v):
 # Endpoint cancellation is evaluated as x*log|x|, including x=0.
 if v.contains(1) or v.contains(-1):
  sg=1
  if v.contains(-1):v=-v;sg=(-1)**(n+1)
  e=v-1;G=-e*(2+e)*(2+e).log()+(2+e)*xlog(e)
  pp=[arb(1),v];rr=[arb(0),arb(1)]
  for k in range(1,n):
   pp.append(((2*k+1)*v*pp[-1]-k*pp[-2])/(k+1));rr.append(((2*k+1)*v*rr[-1]-k*rr[-2])/(k+1))
  value=G*pp[n]/2-(1-v*v)*rr[n]+(v if n==0 else arb(1)/3 if n==1 else 0)
  return sg*value
 if v>1 or v<-1:
  sg=1
  if v<-1:v=-v;sg=(-1)**(n+1)
  # Exact Legendre expansion of (1-x^2)P_n.
  qq=(arb_poly([1,0,-1])*lp(n))
  # power->Legendre via exact triangular polynomial subtraction
  terms=[];work=qq
  for k in range(n+2,-1,-1):
   pk=lp(k);a=work[k]/pk[k]
   if not a.is_zero():terms.append((k,a));work=work-a*pk
  return sg*sum((a*acb(v).legendre_q(k,type=3).real for k,a in terms),arb(0))
 L=((1+v)/(1-v)).log()/2
 pp=[arb(1),v];rr=[arb(0),arb(1)]
 for k in range(1,n):
  pp.append(((2*k+1)*v*pp[-1]-k*pp[-2])/(k+1));rr.append(((2*k+1)*v*rr[-1]-k*rr[-2])/(k+1))
 return (1-v*v)*(pp[n]*L-rr[n])+(v if n==0 else arb(1)/3 if n==1 else 0)

class Source:
 def __init__(self,mom,spins,roy=False,uv=False):
  assert mom['metadata']['witness_sha256']==EXPECTED
  self.m=mom;self.roy=roy;self.uv=uv;self.N=mom['metadata']['kmax'];N=self.N
  assert N>=3 and mom['metadata']['jmax']>=max(spins)
  if roy:assert F(mom['metadata']['B'])==F(119,10) and mom['metadata']['Jmax']==6
  else:assert F(mom['metadata']['T'])==C and (F(mom['metadata']['A']),F(mom['metadata']['B'])) in [(F(12),F(24)),(F(24),F(48))]
  if roy:self.h=F(4);self.rt=F(2);self.ru=F(4);self.tm=F(2);self.norm=F('5.618')
  else:
   md=mom['metadata'];self.h=F(md['u_center_minus']);self.ru=F(md['u_radius']);self.tm=F(md['t_center_minus']);self.rt=F(md['t_radius']);self.norm=F('.2401')
  self.u_den=X+ar(4+self.h);self.t_den=X+ar(4+self.tm);self.q=X+ar(2*C)
  Du=self.u_den**N;Dt=self.t_den**N
  self.polys={}
  for l in spins:
   cs=[coeff(l,j) for j in range(l+1)]
   rs=arb_poly([ar(cs[l-i]*C**(l-i)) for i in range(l+1)])
   ir=arb_poly([cs[l-i]*mom['T'][l-i] for i in range(l+1)])
   pus=[arb_poly([cs[l-i]*mom['M_u'][l-i][k] for i in range(l+1)]) for k in range(N)]
   Nu=geometric_num(pus,self.u_den,N)
   Nt=geometric_num([arb_poly([m]) for m in mom['M_t']],self.t_den,N)
   ddcoef={}
   for j in range(1,l+1):
    for a in range(j):
     r=j-1-a
     for h in range(a+1):
      p=l-j+h;ddcoef[r,p]=ddcoef.get((r,p),F(0))-cs[j]*(-1)**a*math.comb(a,h)*C**(a-h)
   dd=arb_poly([sum((ar(v)*mom['T'][r] for (r,pp),v in ddcoef.items() if pp==p),arb(0)) for p in range(l)])
   pv=arb_poly([])
   if roy:
    for j in range(l):
     for a in range(l-j):pv+=cs[j]*mom['M_s'][j][a]*(X**(l-j-1-a))
   base=(pv+dd)*self.q*Du*Dt-(ir+2*mom['T'][0]*rs)*Du*Dt+self.q*Nu*Dt+self.q*rs*Nt*Du
   den_extra=None
   if uv and not roy:
    smid=F(mom['metadata']['s_center']);sd=X+ar(4-smid);Ds=sd**N
    pvs=[arb_poly([cs[l-i]*mom['M_s'][l-i][k] for i in range(l+1)]) for k in range(N)]
    Ns=geometric_num(pvs,sd,N)
    base=base*Ds+self.q*Du*Dt*Ns;den_extra=(sd,smid,F(mom['metadata']['s_radius']))
   if roy or uv:
    target=l+(3 if uv and not roy else 2)*N-2
    for k in range(target+1,base.degree()+1):
     assert base[k].contains(0),(l,k,base[k])
    base=arb_poly([base[k] for k in range(target+1)])
    assert base.degree()==target, 'UV evaluator assumes this exact degree'
   self.polys[l]=(base,rs,den_extra)
 def eval(self,mu,l):
  delta=mu-4;Q=mu-ar(C);num,rs,extra=self.polys[l];N=self.N
  # Homogenized variable 1/delta is stable throughout the UV.
  if delta>1:
   w=1/delta;nr=tight(rev(num),w);rr=tight(rev(rs),w)
   den=(1+ar(4+self.h)*w)**N*(1+ar(4+self.tm)*w)**N
   degree=num.degree()+2-l-2*N
   if extra:
    sd,smid,srad=extra;den*=(1+ar(4-smid)*w)**N;degree-=N
   val=-(1+2*ar(C)*w)**2*nr/rr/den
   if degree>0:val*=delta**degree
   elif degree<0:val*=w**(-degree)
   rstar=tight(rs,delta)
  else:
   rstar=tight(rs,delta);den=(mu+ar(self.h))**N*(mu+ar(self.tm))**N
   if extra:den*=(mu-ar(extra[1]))**N
   val=-Q*Q*tight(num,delta)/rstar/den
  # Exact geometric-series residual bounds, after normalizing P_l(Z).
  if self.roy and mu.lower()<ar(F(159,20)):
   tmin=F(-79,20)
   rmax=arb_poly([ar(coeff(l,l-i)*tmin**(l-i)) for i in range(l+1)])
   ratio=tight(rmax,delta).abs_upper()/rstar.lower()+delta**l/rstar.lower()
  else:ratio=arb(1)  # valid for all strips, and for physical Roy arguments
  qu=ar(self.ru)/(mu+ar(self.h));qt=ar(self.rt)/(mu+ar(self.tm))
  err=Q**3*ar(self.norm)*(ratio*qu**N/((mu+ar(self.h))*(1-qu))+qt**N/((mu+ar(self.tm))*(1-qt)))
  if extra:
   qs=ar(extra[2])/(mu-ar(extra[1]));err+=Q**3*ar(self.norm)*qs**N/((mu-ar(extra[1]))*(1-qs))
  return val+arb(0,err.upper()),rstar

 def eval_w(self,w,l):
  assert self.roy or self.uv
  num,rs,extra=self.polys[l];N=self.N;qbar=1+2*ar(C)*w
  den=(1+ar(4+self.h)*w)**N*(1+ar(4+self.tm)*w)**N
  if extra:den*=(1+ar(4-extra[1])*w)**N
  val=-qbar*qbar*tight(rev(num),w)/tight(rev(rs),w)/den
  def rem(radius,center):
   return ar(radius)**N*qbar**3*w**(N-2)/((1+ar(4+center)*w)**N*(1+ar(4+center-radius)*w))
  err=ar(self.norm)*(rem(self.ru,self.h)+rem(self.rt,self.tm))
  if extra:err+=ar(self.norm)*rem(extra[2],-extra[1])
  return val+arb(0,err.upper())

class FiniteV:
 def __init__(self,spins=range(8,42,2)):
  self.spins=list(spins);assert all(l>=8 and l%2==0 for l in self.spins);path=ROOT/'output/angles_resummed/angle_000p000.npz';assert hashlib.sha256(path.read_bytes()).hexdigest()==EXPECTED;self.f=RealTail(path,40)
  normcache=pickle.loads((ROOT/'output/certification_tail/compact_source_cache.pkl').read_bytes())
  assert normcache['witness_sha256']==EXPECTED and normcache['roy_abs']<F('5.618') and normcache['strip_abs']<F('.2401')
  self.roy=Source(load(ROOT/'output/certification_head/roy_moments.json'),self.spins,roy=True)
  self.strip=[]
  for A,B in [(12,24),(24,48)]:
   m=load(ROOT/f'output/certification_head/strip_{A}_{B}_moments.json')
   self.strip.append((A,B,Source(m,self.spins),Source(m,self.spins,uv=True)))
  st=np.load(path);self.lx=ar(F(float(st['multipliers'][list(st['names']).index('X4')])));self.hp={};self.cross={}
  for l in self.spins:
   cs=[coeff(l,j) for j in range(l+1)]
   rc=arb_poly([ar(cs[l-i]*(l-i)*C**(l-i-1)) if l-i>=1 else 0 for i in range(l+1)])
   rcc=arb_poly([ar(cs[l-i]*(l-i)*(l-i-1)*C**(l-i-2)) if l-i>=2 else 0 for i in range(l+1)])
   self.cross[l]=((X+ar(2*C))*rcc-3*rc)/arb(2) if False else ((X+ar(2*C))*rcc-3*rc)*(arb(1)/2)
   for j in range(3):
    self.hp[l,j]=arb_poly([ar(F((-1)**(k+j)*math.factorial(l+k),math.factorial(l-k)*math.factorial(k-j)*math.factorial(k+j+1))*C**k) if k>=j else 0 for k in reversed(range(l+1))])
  self.cheb={n:arb_poly([(-1)**0]) for n,v in []}
  cp=[arb_poly([1]),X]
  for n in range(1,max(n for n,v in self.f.h)):cp.append(2*X*cp[-1]-cp[-2])
  self.H=psum(v*cp[n] for n,v in self.f.h)+self.f.c0*X
 def fields(self,mu):
  delta=mu-4;Q=mu-ar(C);x=8/(mu-2)**2-1;y2=1-x*x
  qr=tight(self.H,x)
  for r,c,kind,v in self.f.r:
   assert kind=='cos'
   dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x;nr=1-r*c*x;ni=-r*c
   den=square(2*r*x-c*(1+r*r))+(1-c*c)*square(1-r*r)
   mr=(nr*dr+y2*ni*di)/den;mi=(ni*dr-nr*di)/den;mminus=(1+r*c)/(1+2*r*c+r*r)
   qr+=v*((2*x*x-1)*(mr-mminus)-2*x*y2*mi)
  vf=32/(mu-2)**3*qr+self.f.c0*8/((3*mu-8)*(3*mu-4)*(mu-2))
  vbase=-32*self.f.c0*Q**3/(mu*(3*mu-8)*(3*mu-4))
  strip_pv=[]
  if mu.upper()<101:
   for A,B,src,uvsrc in self.strip:
    zz=2*(mu-A)/(B-A)-1
    vals=[taper_q(n,zz) for n in range(6)]
    weights={j:sum((v*vals[n]*2/(B-A) for jj,n,v in self.f.groups[A,B] if jj==j),arb(0)) for j in range(3)}
    strip_pv.append(weights)
  out={}
  for l in self.spins:
   v,rs=self.roy.eval(mu,l)
   v+=vbase+Q**3*delta**l/rs*vf+self.lx*tight(self.cross[l],delta)/(Q*rs)
   for k,(A,B,src,uvsrc) in enumerate(self.strip):
    if mu.upper()<101:
     vs,_=src.eval(mu,l);pv=sum((strip_pv[k][j]*tight(self.hp[l,j],delta) for j in range(3)),arb(0))
     v+=vs-Q**3/rs*pv
    else:
     vs,_=uvsrc.eval(mu,l);v+=vs
   out[l]=v
  return out

 def uv_fields(self,w):
  x=8*w*w/(1+2*w)**2-1;y2=1-x*x;qbar=1+2*ar(C)*w
  qr=tight(self.H,x)
  for r,c,kind,v in self.f.r:
   dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x;nr=1-r*c*x;ni=-r*c
   den=square(2*r*x-c*(1+r*r))+(1-c*c)*square(1-r*r)
   mr=(nr*dr+y2*ni*di)/den;mi=(ni*dr-nr*di)/den;mminus=(1+r*c)/(1+2*r*c+r*r)
   qr+=v*((2*x*x-1)*(mr-mminus)-2*x*y2*mi)
  base=-ar(F(32,9))*self.f.c0*qbar*qbar/((1+4*w)*(1+ar(C)*w))
  vf=32*(qbar/(1+2*w))**3*qr+ar(F(8,9))*self.f.c0*qbar*qbar/((1+ar(C)*w)*(1+2*w))
  out={}
  for l in self.spins:
   rs=tight(rev(self.roy.polys[l][1]),w)
   v=self.roy.eval_w(w,l)+base+vf/rs+self.lx*w*w*tight(rev(self.cross[l],l-1),w)/(qbar*rs)
   for A,B,src,uvsrc in self.strip:v+=uvsrc.eval_w(w,l)
   out[l]=v
  return out

if __name__=='__main__':
 import numpy as np
 ctx.prec=240;tic=time.time();bank=FiniteV()
 for m in [4,4.1,5,8,12,24,48,49,100,101,1000,10000]:
  vals=bank.fields(ar(F(str(m))));print(m,{l:str(v) for l,v in vals.items() if l in[8,20,40]},flush=True)
 print('seconds',time.time()-tic)

