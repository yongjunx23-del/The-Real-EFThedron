"""Component certificate for V_l/p_l <= 0 at all even l>=L, mu>=M.

All source integrals are enclosed by exact rational Bernstein polynomials.
Finite-mu errors relative to the joint UV kernel are explicit rational bounds.
The compact tau interval uses Arb; tau>=8 closes by monotone Legendre bounds.
This is NOT a certificate for mu<M or for the complete dual objective.
"""
import sys,re,json,math,hashlib,argparse,pickle
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import sympy as sy
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import arb,ctx
from certify_roy_dangerous_sign import split_axis
ctx.prec=160

def ar(x):
 x=F(x);return arb(x.numerator)/x.denominator

def bernstein(expr,x,y):
 p=sy.Poly(sy.expand(expr),x,y);n,m=p.degree(x),p.degree(y)
 a=[[F(p.coeff_monomial(x**i*y**j)) for j in range(m+1)] for i in range(n+1)]
 return [[sum((a[k][l]*F(math.comb(i,k),math.comb(n,k))*F(math.comb(j,l),math.comb(m,l)) for k in range(i+1) for l in range(j+1)),F(0)) for j in range(m+1)] for i in range(n+1)]

def integral(b):return sum((sum(row,F(0)) for row in b),F(0))/len(b)/len(b[0])
def split_box(b,depth=8):
 stack=[(b,F(0),F(1),F(0),F(1),0)]
 while stack:
  bs,x0,x1,y0,y1,d=stack.pop()
  if d==depth:yield bs,x0,x1,y0,y1;continue
  l,r=split_axis(bs,d%2)
  if d%2==0:
   mid=(x0+x1)/2;stack.extend([(l,x0,mid,y0,y1,d+1),(r,mid,x1,y0,y1,d+1)])
  else:
   mid=(y0+y1)/2;stack.extend([(l,x0,x1,y0,mid,d+1),(r,x0,x1,mid,y1,d+1)])

def source_data(path,depth):
 st=np.load(path);assert float(st['direction'][1])==0; x,y=sy.symbols('x y');c=sy.Rational(4,3)
 roy={};strip={};alpha=F(st['direction'][0]);harm=[];rats=[];lx=None
 for name,val in zip(st['names'],st['multipliers']):
  v=F(float(val));m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',str(name))
  if m:
   key=F(m[1]);J,n=int(m[2]),int(m[3]);assert 4<key<=12 and J%2==0;roy[key]=roy.get(key,0)+sy.Rational(v)*sy.legendre(n,2*x-1)*sy.legendre(J,y)*2
   if J==n==0:alpha-=2*v
   continue
  m=re.fullmatch(r'Taper(\d+)_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',str(name))
  if m:
   power,A,B,j,n=int(m[1]),F(m[2]),F(m[3]),int(m[5]),int(m[6]);assert 4<A<B<=48 and abs(float(m[4])-4/3)<1e-5;key=(A,B)
   expr=(1-(2*x-1)**2)**power*sy.legendre(n,2*x-1)*sy.legendre(j,2*y-1)
   strip[key]=strip.get(key,0)+sy.Rational(v)*expr
   if j==0:alpha-=v*F(sy.integrate(expr,(x,0,1)))
   continue
  m=re.fullmatch(r'H(\d+)',str(name))
  if m:harm.append((int(m[1]),v));continue
  m=re.fullmatch(r'Rc_([^_]+)_([^_]+)',str(name))
  if m:
   assert 0<=F(m[1])<1 and abs(F(m[2]))<=1
   rats.append((F(m[1]),F(m[2]),v));continue
  if name=='X4':lx=v;continue
  raise ValueError('Unsupported '+str(name))
 specs=[]
 for B,W in roy.items():
  s=4+(sy.Rational(B)-4)*x;t=-(s-4)*(1-y)/2
  specs.append((W,s,t,lambda a,b,d,e,i,j,n,m,B=B:(B-4)/2*(a+(b-a)*F(i+1,n+2))*(1-d-(e-d)*F(j+1,m+2))))
 for (A,B),W in strip.items():
  s=sy.Rational(A)+(sy.Rational(B)-sy.Rational(A))*x;t=-c*(1-y)
  specs.append((W,s,t,lambda a,b,d,e,i,j,n,m:F(4,3)*(1-d-(e-d)*F(j+1,m+2))))
 out={'alpha':alpha,'harm':harm,'rats':rats,'lx':lx,'signed1':F(0),'signed2':F(0),'norm1':F(0),'norm2':F(0),'correction':[],'correction_integrated_t':True}
 for W,s,t,tmax in specs:
  u=4-s-t;K1=2*(s-c)*(c-u);K2=2*(t-c)**2
  for tag,K in [('1',K1),('2',K2)]:
   bs=bernstein(W*K,x,y);out['signed'+tag]+=integral(bs)
   for box,a,b,d,e in split_box(bs,depth):
    area=(b-a)*(e-d);den=len(box)*len(box[0])
    out['norm'+tag]+=area*sum((sum((abs(q) for q in row),F(0)) for row in box),F(0))/den
    if tag=='1':
     positive=area*sum((sum((max(q,F(0)) for q in row),F(0)) for row in box),F(0))/den
     if positive:
      moment=area*sum((max(q,F(0))*tmax(a,b,d,e,i,j,len(box)-1,len(box[0])-1) for i,row in enumerate(box) for j,q in enumerate(row)),F(0))/den
      out['correction'].append((positive,moment))
 return out

def errors(d,M):
 c=F(4,3);alpha=abs(d['alpha']);M=F(M)
 eps1=max(1-(1-4/M)*(1-c/M)**2/(1+44/M),1/((1-48/M)*(1-F(20,3)/M))-1)
 # Expanded numerator of R1-1 cancels the large O(s/M) terms exactly.
 eps1=min(eps1,(8/M+2153/M**2+14088/M**3)/((1-48/M)*(1-F(20,3)/M)))
 eps2=max(1-(1-c/M)**2/(1+4/M),1/(1-F(20,3)/M)-1)
 ea=F(32,9)*alpha*c*c/(M*(M-2*c))+eps2*d['norm2']
 rf=(1+F(2,3)/(M-2))**3
 signed_limit=-d['alpha']+sum((v*(-1)**j for j,v in d['harm']),F(0))
 ef=32*(rf-1)*abs(signed_limit)
 ef+=32*rf*8/(M-2)**2*(alpha+sum((abs(v)*j*j for j,v in d['harm']),F(0)))
 rg=max(1-(1-c/M)**2,1/((1-2*c/M)*(1-2/M))-1)
 ef+=F(8,9)*alpha*rg
 for r,co,v in d['rats']:
  # Real UV cancellation: Chebyshev Lipschitz bound at Re(omega)=-1.
  # sum r^k[(k+2)^2+4], using |T_k(cosine)|<=1.
  lipschitz=r*(1+r)/(1-r)**3+4*r/(1-r)**2+8/(1-r)
  ef+=abs(v)*32*rf*8/(M-2)**2*lipschitz
 return ea,ef+eps1*d['norm1']

def leg_coeffs(L):
 T=L*(L+1)
 # P_L(1+2a/T) has positive coefficient product_{j<k}(T-j(j+1))/(T^k(k!)^2).
 cs=[F(1)]
 for k in range(1,L+1):cs.append(cs[-1]*F(T-(k-1)*k,T*k*k))
 return cs

def poly(cs,a):
 out=arb(0)
 for c in reversed(cs):out=out*a+ar(c)
 return out

def certify(d,L,M,step,dynamic=False):
 assert L>=8 and L%2==0 and M>48 and d['lx']<0 and d['correction_integrated_t']
 ea,en=errors(d,M);A=-F(32,9)*d['alpha']-d['signed2']+ea
 B=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+en
 coeff=leg_coeffs(L);d2=[coeff[k]*k*(k-1) for k in range(2,len(coeff))]
 # At tau>=8 use positivity of X4 and monotonicity of P_l at fixed tau.
 tau_end=F(256) if dynamic else F(8)
 pp=poly(coeff,ar(F(4,3)*tau_end))
 Nmax=B+2*sum((v for v,t in d['correction']),F(0))
 big=ar(A)+ar(max(Nmax,F(0)))/pp
 if not big<0:raise RuntimeError('Large-tau closure failed '+str(big))
 worst=arb(-1000000);worstwhere=None;cells=0
 cuts=[F(0)];step=F(str(step))
 while cuts[-1]<8:cuts.append(min(F(8),cuts[-1]+step))
 if dynamic:cuts.extend(F(i) for i in range(9,257))
 for lo,hi in zip(cuts,cuts[1:]):
  mcell=max(F(M),F(4)+F(L*(L+1))/hi) if dynamic else F(M)
  eac,enc=errors(d,mcell) if dynamic else (ea,en)
  Ac=-F(32,9)*d['alpha']-d['signed2']+eac
  Bc=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+enc
  # Monotone factor bounds avoid dependency broadening in bivariate expressions.
  al=ar(F(4,3)*lo);ah=ar(F(4,3)*hi)
  pl=poly(coeff,al);pu=(2*ah.sqrt()).bessel_i(0)
  cor=sum((min(2*v,hi*t) for v,t in d['correction']),F(0))
  num=Bc+cor
  # If numerator is negative its worst ratio uses the largest possible denominator.
  ratio=ar(num)/(pl if num>=0 else pu)
  if lo==0 or lo>=8:xlo=arb(0)
  else:
   p2=poly(d2,al)
   # Uniform upper bound for d/da 0F1(1;a) = I1(2sqrt(a))/sqrt(a).
   pa=(2*ah.sqrt()).bessel_i(1)/ah.sqrt()
   neg=ar(F(3,2)*hi/(mcell-F(4,3)))*pa/pl
   raw=ar(lo*lo/F(2))*p2/pu-neg
   xlo=raw if raw>0 else arb(0)
  val=ar(Ac)+ratio+ar(d['lx'])*xlo
  cells+=1
  if not val<0:raise RuntimeError(f'Uncertified cell {lo}..{hi}: {val}')
  if val.upper()>worst.upper():worst=val;worstwhere=[float(lo),float(hi)]
 return dict(A_upper=str(ar(A)),B_with_errors=str(ar(B)),A_error=str(ar(ea)),N_error=str(ar(en)),
  small_tau_max_upper=str(worst.upper()),small_tau_worst_cell=worstwhere,small_tau_cells=cells,
  large_tau_upper=str(big.upper()),norm1_upper=str(ar(d['norm1'])),norm2_upper=str(ar(d['norm2'])),
  positive_K1_upper=str(ar(sum((v for v,t in d['correction']),F(0)))))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,default=ROOT/'output/angles_resummed/angle_000p000.npz');ap.add_argument('--L',type=int,default=102);ap.add_argument('--M',type=int,default=1000000);ap.add_argument('--depth',type=int,default=8);ap.add_argument('--step',default='.005');ap.add_argument('--dynamic',action='store_true');ap.add_argument('--output',type=Path,default=ROOT/'output/certification_tail/uv_component.json');a=ap.parse_args()
 cache=ROOT/f'output/certification_tail/source_cache_v2_depth{a.depth}.pkl'
 sha=hashlib.sha256(a.witness.read_bytes()).hexdigest()
 builder_sha=hashlib.sha256(Path(__file__).read_bytes()+(ROOT/'research/certify_roy_dangerous_sign.py').read_bytes()).hexdigest()
 if cache.exists():
  saved=pickle.loads(cache.read_bytes());d=saved['data'] if saved['sha']==sha and saved.get('builder_sha256')==builder_sha else source_data(a.witness,a.depth)
 else:d=source_data(a.witness,a.depth)
 cache.write_bytes(pickle.dumps({'sha':sha,'builder_sha256':builder_sha,'data':d}))
 res=certify(d,a.L,a.M,a.step,a.dynamic);res.update(witness=str(a.witness),witness_sha256=sha,L=a.L,M=a.M,builder_sha256=builder_sha,bernstein_depth=a.depth,dynamic_mu_constraint=a.dynamic,arithmetic='Exact rational source Bernstein enclosures; Arb160 scalar bounds',domain=f'all even l>={a.L}, all mu>={a.M}',claim='V_l(mu)/p_l(mu)<0 on the stated domain only',limitation='Does not control 4<mu<M or finite-spin head quadrature')
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
if __name__=='__main__':main()
