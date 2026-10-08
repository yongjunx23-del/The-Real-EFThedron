"""Arb certificate of a compact part of the *joint limiting* V/p function.

This does not certify finite-spin approach to the limit, nor a physical bound.
Exact binary64 multipliers; Roy endpoint119/10 and strip T=4/3 exact.
Source integrals reduced to exact moments; Taylor remainder uses |J_q(x)|
<= (x/2)^q/q!, a consequence of Poisson's integral representation.
"""
import sys,re,json,math,hashlib,argparse,time,pickle
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sy
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import arb,ctx
from crossing_kernels import crossing_catalogue

def ar(x):
 x=F(x);return arb(x.numerator)/x.denominator
def pw(x,n):
 out=arb(1)
 for _ in range(n):out*=x
 return out
def beta(b,k):return F(math.factorial(b)*math.factorial(k),math.factorial(b+k+1))
def derive(path,N=200):
 st=np.load(path);lam={str(n):F(float(v)) for n,v in zip(st['names'],st['multipliers'])}
 x,z=sy.symbols('x z');c=sy.Rational(4,3);roy={};strip={};alpha=F(float(st['direction'][0]));polys=[]
 for n,v in lam.items():
  m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',n)
  if m:
   B=F(m[1]);J,k=int(m[2]),int(m[3]);roy[B]=roy.get(B,0)+2*sy.Rational(v)*sy.legendre(k,2*x-1)*sy.legendre(J,z)
   if J==k==0:alpha-=2*v
  m=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',n)
  if m:
   A,B=F(m[1]),F(m[2]);j,k=int(m[4]),int(m[5]);key=(A,B)
   e=sy.Rational(v)*(1-(2*x-1)**2)*sy.legendre(k,2*x-1)*sy.legendre(j,2*z-1)
   strip[key]=strip.get(key,0)+e
   if j==0:alpha-=F(sy.integrate(e,(x,0,1),(z,0,1)))
 for B,W in roy.items():
  s=4+sy.Rational(B-4)*x;t=-sy.Rational((B-4)/2)*x*(1-z)
  polys.append((W,s,t,True,F(B-4)/2))
 for (A,B),W in strip.items():
  s=sy.Rational(A)+sy.Rational(B-A)*x;t=-c*(1-z)
  polys.append((W,s,t,False,F(4,3)))
 moments=[F(0)]*(N+7);absMom=[F(0)]*7;Hs=F(0);C=F(0);tmax=F(0)
 for W,s,t,isroy,q in polys:
  K=2*(s-c)*(c-(4-s-t));H=2*(t-c)**2
  terms=[(a,b,F(co)) for (a,b),co in sy.Poly(sy.expand(W*K),x,z).terms()]
  Hs+=F(sy.integrate(sy.expand(W*H),(x,0,1),(z,0,1)))
  for k in range(len(moments)):
   moments[k]+=q**k*sum((co*beta(b,k)/F(a+k+1 if isroy else a+1) for a,b,co in terms),F(0))
  for k in range(7):
   absMom[k]+=q**k*sum((abs(co)*beta(b,k)/F(a+k+1 if isroy else a+1) for a,b,co in terms),F(0))
  C+=sum((abs(co)/F((a+1)*(b+1)) for a,b,co in terms),F(0));tmax=max(tmax,q)
 return dict(lam=lam,alpha=alpha,moments=moments,absMom=absMom,C=C,tmax=tmax,Hs=Hs,N=N)

class Certificate:
 def __init__(self,data):
  self.d=data;self.lam=data['lam'];self.cross=crossing_catalogue(8)
  N=data['N'];self.co=[]
  for q in range(6):
   self.co.append([ar((-1)**(k+q)*data['moments'][k+q]/F(math.factorial(k)*math.factorial(k+q)*math.factorial(q))) for k in range(N)])
 def source(self,lo,hi):
  d=self.d;N=d['N'];mid=(lo+hi)/2;h=(hi-lo)/2;ma=ar(mid);dx=arb(0,ar(h));out=arb(0)
  for q in range(5,-1,-1):
   val=arb(0)
   for co in reversed(self.co[q]):val=val*ma+co
   a=d['tmax']*mid;ratio=a/F((N+1)*(N+q+1));assert ratio<1
   rem=d['C']*d['tmax']**q*a**N/F(math.factorial(N)*math.factorial(N+q)*math.factorial(q))/(1-ratio)
   val+=arb(0,ar(rem));out=out*dx+val
  remainder=d['absMom'][6]*h**6/F(math.factorial(6)**2)
  return out+arb(0,ar(remainder))
 def field(self,lo,hi):
  d=self.d;ta=arb(ar((lo+hi)/2),ar((hi-lo)/2));c=F(4,3);P=(ar(c)*ta).hypgeom_0f1(1)
  forward=-F(280,9)*d['alpha']+sum((32*v*(-1)**int(n[1:]) for n,v in self.lam.items() if re.fullmatch(r'H\d+',n)),F(0))
  result=ar(-F(32,9)*d['alpha']-d['Hs'])+(ar(forward)-self.source(lo,hi))/P
  for row in self.cross:
   j=row['degree']-2;b=self.lam[row['name']]*F(row['beta'][j])
   result+=ar(b/F(math.factorial(j)**2))*pw(ta,j)*(ar(c)*ta).hypgeom_0f1(j+1)/P
  for n,val in self.lam.items():
   if not n.startswith('Dt'):continue
   p=int(n[2]);j=int(n.split('_j')[1]) if '_j' in n else 1
   core=ar(F(32,math.factorial(p)))*pw(ta,p)+ar(F(16*p,math.factorial(p-1)))*pw(ta,p-1)
   if p>=2:core+=ar(F(2*p*(p-1),math.factorial(p-2)))*pw(ta,p-2)
   term=ar((-1)**j*4**(p-1))*core/P
   if n.endswith('G_regular'):
    ka=ar(F(-8,9*math.factorial(p)))*pw(ta,p)+ar(F(4*p,3*math.factorial(p-1)))*pw(ta,p-1)
    if p>=2:ka-=ar(F(p*(p-1),2*math.factorial(p-2)))*pw(ta,p-2)
    kb={1:F(-16,3),2:F(4),3:F(0)}[p]
    term-=ka/P+ar(kb)
   result+=ar(val)*term
  return result

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--prec',type=int,default=256);ap.add_argument('--max-tau',default='256');args=ap.parse_args()
 ctx.prec=args.prec;tic=time.time();w=ROOT/'output/angle_constraints/ming0_scaled_threshold_repaired.npz'
 folder=ROOT/'output/lower_g0_upper_g2_20261006';cache=folder/'tail_joint_exact_moments.pkl'
 if cache.exists():d=pickle.loads(cache.read_bytes())
 else:d=derive(w);cache.write_bytes(pickle.dumps(d));print('exact moments ready',flush=True)
 cert=Certificate(d);end=F(args.max_tau);stack=[(F(0),end,0)];leaves=[];failed=[];count=0
 while stack:
  lo,hi,depth=stack.pop();v=cert.field(lo,hi);count+=1
  if v.is_finite() and v<0:leaves.append({'lo':str(lo),'hi':str(hi),'upper':str(v.upper())});continue
  if depth>=30:failed.append({'lo':str(lo),'hi':str(hi),'value':str(v)});continue
  mid=(lo+hi)/2;stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
  if count%1000==0:print('cells',count,'pending',len(stack),'leaves',len(leaves),'failed',len(failed),flush=True)
 out={'witness_sha256':hashlib.sha256(w.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'precision_bits':args.prec,'domain':['0',str(end)],'status':'PASS' if not failed else 'FAIL','leaves':leaves,'failed':failed,
  'cell_evaluations':count,'seconds':time.time()-tic,'claim':'Compact-tau joint UV limiting function only; finite-spin convergence is NOT certified.'}
 p=folder/f'tail_joint_limit_certificate_{args.prec}.json';p.write_text(json.dumps(out,indent=2))
 print(json.dumps({k:v for k,v in out.items() if k not in ('leaves','failed')},indent=2))
