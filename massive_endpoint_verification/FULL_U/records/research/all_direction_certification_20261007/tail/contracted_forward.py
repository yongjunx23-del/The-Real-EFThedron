"""Exact contraction of the real harmonic and transfer-derivative numerator.

Chebyshev identities eliminate complex interval powers before coefficients are
combined. This is an evaluation refinement, not an additional physical null.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,math,re,pickle,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research/lower_g0_upper_g2_20261006')]
from flint import fmpq,fmpq_poly,arb,acb
from tail_finite_uv_region import Region,ar,pw,jet,add,scale,mul,inv,jpow

def fq(x):x=F(x);return fmpq(x.numerator,x.denominator)
W=fmpq_poly([0,1]);Q0=1+2*W;QQ=1+fq(F(8,3))*W
ZERO=fmpq_poly([])
def ra(n=0,d=0):return(fmpq_poly(n) if isinstance(n,list) else n if isinstance(n,fmpq_poly) else fmpq_poly([fq(n)]),d)
def radd(a,b):
 n=max(a[1],b[1]);return(a[0]*Q0**(n-a[1])+b[0]*Q0**(n-b[1]),n)
def rscale(a,c):return(a[0]*fq(c),a[1])
def rmul(a,b):return(a[0]*b[0],a[1]+b[1])
def jadd(a,b):return[radd(x,y)for x,y in zip(a,b)]
def jscale(a,c):return[rscale(x,c)for x in a]
def jmul(a,b):return[sumr(rmul(a[k],b[j-k])for k in range(j+1))for j in range(4)]
def sumr(items):
 out=ra()
 for v in items:out=radd(out,v)
 return out
def jconst(a):return[ra(a),ra(),ra(),ra()]
def shifted_eval(poly,lo,hi):
 mid=(lo+hi)/2;rad=(hi-lo)/2;p=poly(fmpq_poly([fq(mid),1]));x=arb(0,ar(rad));out=arb(0)
 for co in reversed(p.coeffs()):out=out*x+arb(co)
 return out
def bounds(x):
 def f(v):m,e=v.man_exp();return F(int(m))*F(2)**int(e)
 return f(x.lower()),f(x.upper())

class ContractedRegion(Region):
 def __init__(self,*args):
  super().__init__(*args)
  A2=[ra(4),ra(2),ra(F(1,4)),ra()]
  q2=[ra((-1)**n*(n+1)*(W/2)**n,n+2)for n in range(4)]
  xx=jadd(jconst(-1),jscale(jmul(A2,q2),2))
  # The correction 2 A^2 w^2/q^2 contains w^2, not a transfer power.
  xx=jadd(jconst(-1),[(a[0]*W*W,a[1])for a in jscale(jmul(A2,q2),2)])
  q3=[ra((-1)**n*math.comb(n+2,2)*(W/2)**n,n+3)for n in range(4)]
  ff=[(8*QQ**3*a[0],a[1])for a in jmul(A2,q3)]
  cheb=[jconst(1),xx]
  for j in range(1,40):cheb.append(jadd(jscale(jmul(xx,cheb[-1]),2),jscale(cheb[-2],-1)))
  fc=[jmul(ff,p)for p in cheb]
  coeff=[ra(),ra(),ra(),ra()]
  coeff[0]=rscale(fc[1][0],self.d['alpha'])
  for name,val in self.lam.items():
   mm=re.fullmatch(r'H(\d+)',name)
   if mm:coeff[0]=radd(coeff[0],rscale(fc[int(mm[1])][0],val));continue
   if not name.startswith('Dt'):continue
   p=int(name[2]);j=int(name.split('_j')[1])if'_j'in name else 1
   fj=[]
   for k in range(4):fj.append(sumr(rscale(fc[j+h][k],(-1)**h*math.comb(2*p-2,h))for h in range(2*p-1)))
   for k in range(p+1):coeff[k]=radd(coeff[k],rscale(fj[p-k],val*math.factorial(p)))
  self.contract=coeff
  # The rational rows vanish quadratically in w in their real part. Evaluating
  # omega+1 as a complex interval loses this cancellation and was the largest
  # remaining enclosure error near the joint limit.
  X=fmpq_poly([0,1]);tc=[fmpq_poly([1]),X]
  for j in range(1,6):tc.append(2*X*tc[-1]-tc[-2])
  nx=2*(2*W)**2-Q0**2;dx=Q0**2
  def homogeneous(pol):
   n=pol.degree();return sum((pol[k]*nx**k*dx**(n-k)for k in range(n+1)),ZERO),n
  self.rats=[]
  for name,val in self.lam.items():
   mm=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',name)
   if not mm:continue
   kind,r,c=mm[1],fq(F(mm[2])),fq(F(mm[3]));dm=1+2*r*c+r*r
   nn=X**2*r*(X+1)*((c+r-r*(1+r*c)*X)if kind=='cos'else(1+r*r*X))
   dd=dm*(1-2*r*c*X+r*r*X*X)
   realn=sum((nn[i]*dd[j]*tc[abs(i-j)]for i in range(nn.degree()+1)for j in range(dd.degree()+1)),ZERO)
   reald=sum((dd[i]*dd[j]*tc[abs(i-j)]for i in range(dd.degree()+1)for j in range(dd.degree()+1)),ZERO)
   np,nd=homogeneous(realn);dp,ddg=homogeneous(reald)
   # FF=32 Qbar^3/q0^3 and degrees can differ.
   powq=3+2*nd-2*ddg
   if powq>=0:dp*=Q0**powq
   else:np*=Q0**(-powq)
   self.rats.append((fq(val)*32*QQ**3*np,dp))
 def numerator(self,ta,et,w):
  lo,hi=bounds(w);q=1+2*w;Q=1+ar(F(8,3))*w
  pt=[arb(1),ta,pw(ta,2)*(1-2*et)/4,pw(ta,3)*(1-2*et)*(1-6*et)/36]
  num=sum((shifted_eval(pol,lo,hi)/pw(q,n)*pt[k]for k,(pol,n)in enumerate(self.contract)),arb(0))
  alpha=ar(self.d['alpha'])
  num+=alpha*8*pw(Q,3)/((3+4*w)*(3+8*w)*(1+2*w))
  outside=alpha*(-32*pw(Q,3)/((1+4*w)*(3+4*w)*(3+8*w)))-ar(self.d['Hs'])
  # Small rational subtraction jets remain separate from the large harmonic
  # contraction. All their exact finite-w factors are retained.
  t=[acb(0),acb(1),acb(0),acb(0)]
  poly=jpow(add(scale(t,3),jet(-4)),2);common=scale(poly,-pw(8*w+3,2)/27)
  ka=mul(common,inv(mul(add(scale(t,w),jet(4*w+2)),add(scale(t,3*w),jet(4*w+3)))))
  kb=mul(scale(common,2),inv(mul(add(scale(t,w),jet(-4*w-1)),add(scale(t,3*w),jet(4*w+3)))))
  pka=mul([acb(v)for v in pt],ka)
  for p in [1,2,3]:
   val=self.lam.get(f'Dt{p}G_regular',F(0));vv=ar(val)*math.factorial(p)
   num-=vv*pka[p].real;outside-=vv*kb[p].real
  for np,dp in self.rats:num+=shifted_eval(np,lo,hi)/shifted_eval(dp,lo,hi)
  return outside,num
