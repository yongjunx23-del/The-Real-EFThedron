"""Validated real-kernel norm for a FINITE spin block, no V claim.

Strip angular projections are exact polynomials. Each polynomial is composed
with a centered affine variable before interval evaluation, preserving
cancellations. Positive interval upper Riemann sums cover beta in[0,1].
"""
import sys,math,json,argparse,time,heapq,hashlib
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from certify_real_tail import RealTail,exact,power,square,dyadic
from flint import arb,arb_poly,ctx

def lp(n):
 ps=[arb_poly([1]),arb_poly([0,1])]
 for k in range(1,n):ps.append(((2*k+1)*ps[-1]*arb_poly([0,1])-k*ps[-2])*(arb(1)/(k+1)))
 return ps[n]
def tight(poly,x):
 return poly(arb_poly([x.mid(),1]))(arb(0,x.rad()))

class FiniteReal(RealTail):
 def __init__(self,path,first=8,last=40):
  super().__init__(path,40)
  assert first>6 and first%2==last%2==0 and first<=last
  self.spins=list(range(first,last+1,2));self.hp={};self.energy={}
  for (A,B),terms in self.groups.items():
   for j in {j for j,n,v in terms}:
    self.energy[A,B,j]=sum((v*lp(n) for jj,n,v in terms if jj==j),arb_poly([]))
    for l in self.spins:
     coeff=[exact(F((-1)**(k+j)*math.factorial(l+k),math.factorial(l-k)*math.factorial(k-j)*math.factorial(k+j+1))) if k>=j else arb(0) for k in range(l+1)]
     self.hp[l,j]=arb_poly(coeff)
 def density(self,a,b):
  beta=arb(exact((a+b)/2),exact((b-a)/2));x=beta*beta;one=1-x
  # 1/P_l(Z) increases with beta, so its value at the right endpoint is an upper bound.
  xb=exact(b)*exact(b);zb=(2+xb)/(3*xb)
  inv={l:(1/zb.legendre_p(l)).upper() for l in self.spins}
  sa=exact(4)/(1-exact(a)**2);sb=exact(4)/(1-exact(b)**2) if b<1 else None
  active=[(A,B,terms) for (A,B),terms in self.groups.items() if not(sa>=exact(B) or (sb is not None and sb<=exact(A)))]
  kval=self.K(beta)
  if not active:
   common=exact(65536)/(27*arb.pi())*x*(2+x)**3/(1+x)**10*power(one,3)*square(kval)
   return (common*sum(((2*l+1)*inv[l] for l in self.spins),arb(0))).upper()
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
    us[l]=us[l]+val
  density=Q**3*8/(2*arb.pi()*power(one,2))*sum(((2*l+1)*square(us[l])*inv[l] for l in self.spins),arb(0))
  return density.upper()
 def integrate(self,target,max_cells=12000):
  heap=[];serial=0;total=arb(0)
  def ins(a,b):
   nonlocal serial,total
   h=self.density(a,b);mass=exact(b-a)*h
   assert mass.is_finite()
   total+=mass;heapq.heappush(heap,(-float(mass),serial,a,b,mass));serial+=1
  for i in range(256):ins(F(i,256),F(i+1,256))
  while len(heap)<max_cells and not total<arb(target):
   _,_,a,b,m=heapq.heappop(heap);total-=m;c=(a+b)/2;ins(a,c);ins(c,b)
  return total.upper(),[{'a':str(a),'b':str(b),'mass_upper':dyadic(m)} for _,_,a,b,m in sorted(heap,key=lambda x:x[2])]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--first',type=int,default=8);ap.add_argument('--last',type=int,default=40);ap.add_argument('--target',default='.25');ap.add_argument('--max-cells',type=int,default=12000);ap.add_argument('--output',type=Path,default=ROOT/'output/certification_tail/finite_A8_40.json');a=ap.parse_args();ctx.prec=180
 path=ROOT/'output/angles_resummed/angle_000p000.npz';tic=time.time();bank=FiniteReal(path,a.first,a.last);value,cells=bank.integrate(a.target,a.max_cells)
 out={'witness_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'spins':bank.spins,'A_upper':str(value),'A_upper_dyadic':dyadic(value),'cells':len(cells),'target_met':bool(value<arb(a.target)),'seconds':time.time()-tic,'claim':'Finite-spin real U norm only; no V sign or complete head integral claim.','partition':cells}
 a.output.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='partition'},indent=2))
if __name__=='__main__':main()

