"""Scoped all-spin V/p upper enclosure in a bounded joint UV region.

All even l>=L, tau=l(l+1)/(mu-4) in an explicitly recorded interval.
This is not a complete dual coupling certificate.
"""
import sys,json,pickle,math,hashlib,argparse,time,re
from fractions import Fraction as F
from pathlib import Path
import sympy as sy
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(Path(__file__).resolve().parent),str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from flint import arb,acb,ctx
from tail_certify_joint_limit import Certificate,ar,pw
from matrix_crossing_series import coefficient_polynomial,shifted_certificate

N=3
def jet(x):return [acb(x)]+[acb(0)]*N
def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,x):return [t*x for t in a]
def mul(a,b):return [sum((a[k]*b[j-k] for k in range(j+1)),acb(0)) for j in range(N+1)]
def inv(a):
 out=[1/a[0]]
 for n in range(1,N+1):out.append(-sum((a[k]*out[n-k] for k in range(1,n+1)),acb(0))/a[0])
 return out
def jpow(a,k):
 out=jet(1)
 for _ in range(k):out=mul(out,a)
 return out
def ball(lo,hi):return (lo+hi)/2+arb(0,(hi-lo)/2)

class Region:
 def __init__(self,witness,moments,norms,L):
  self.witness=Path(witness);self.sha=hashlib.sha256(self.witness.read_bytes()).hexdigest();self.L=L;self.T=F(L*(L+1));self.K=64
  stored=pickle.loads(Path(moments).read_bytes());assert stored['witness_sha256']==self.sha;self.d=stored['data'];self.lam=self.d['lam'];self.model=Certificate(self.d)
  nr=json.loads(Path(norms).read_text());assert nr['sha256']==self.sha;self.norm1=F(nr['norm1_rational']);self.norm2=F(nr['norm2_rational'])
  p=coefficient_polynomial(self.witness);k,z=p.gens;assert shifted_certificate(p,10,1) is not None
  self.fc=[]
  for j in range(self.K):
   q=sy.Poly(p.as_expr().subs(k,j),z);self.fc.append([ar(F(q.nth(i))) for i in range(max(0,q.degree())+1)])
  self.pc=[];pr=F(1)
  for k in range(self.K):
   if k:pr*=1-F((k-1)*k,self.T)
   assert pr>0
   self.pc.append(pr*F(4,3)**k/F(math.factorial(k)**2))
 def numerator(self,ta,et,w):
  # All operations are analytic at w=0; powers of the spectral energy
  # have already canceled from the scaled kernels Q^3 U and Q^3 V.
  t=[acb(0),acb(1),acb(0),acb(0)];A=add(jet(2),scale(t,ar(F(1,2))));q=add(jet(1+2*w),scale(t,w/2))
  Q=1+ar(F(8,3))*w
  root=[]
  for k in range(4):root.append(acb((1+4*w).sqrt()*ar(F(sy.binomial(sy.Rational(1,2),k)))*pw(w,k)/pw(1+4*w,k)))
  zeta=add(scale(A,w),scale(root,acb(0,1)));om=mul(jpow(zeta,2),jpow(inv(q),2))
  FF=scale(mul(jpow(A,2),jpow(inv(q),3)),8*pw(Q,3));op=[jet(1)]
  for j in range(1,41):op.append(mul(op[-1],om))
  pt=[acb(1),acb(ta),acb(pw(ta,2)*(1-2*et)/4),acb(pw(ta,3)*(1-2*et)*(1-6*et)/36)]
  diff=add(jet(1),scale(om,-1));dp={p:jpow(diff,2*p-2) for p in [1,2,3]}
  poly=jpow(add(scale(t,3),jet(-4)),2);common=scale(poly,-pw(8*w+3,2)/27)
  ka=mul(common,inv(mul(add(scale(t,w),jet(4*w+2)),add(scale(t,3*w),jet(4*w+3)))))
  kb=mul(scale(common,2),inv(mul(add(scale(t,w),jet(-4*w-1)),add(scale(t,3*w),jet(4*w+3)))))
  pka=mul(pt,ka)
  alpha=ar(self.d['alpha'])
  outside=alpha*(-32*pw(Q,3)/((1+4*w)*(3+4*w)*(3+8*w)))-ar(self.d['Hs'])
  num=alpha*(FF[0].real*om[0].real+8*pw(Q,3)/((3+4*w)*(3+8*w)*(1+2*w)))
  for name,val in self.lam.items():
   vv=ar(val);m=re.fullmatch(r'H(\d+)',name)
   if m:num+=vv*(FF[0]*op[int(m[1])][0]).real;continue
   m=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',name)
   if m:
    kind,r,c=m[1],ar(F(m[2])),ar(F(m[3]));oo=om[0];den=1-2*r*c*oo+r*r*oo*oo;dm=1+2*r*c+r*r
    diff0=r*(oo+1)*(c+r-r*(1+r*c)*oo)/(den*dm) if kind=='cos' else r*(oo+1)*(1+r*r*oo)/(den*dm)
    num+=vv*(FF[0]*oo*oo*diff0).real;continue
   if not name.startswith('Dt'):continue
   p=int(name[2]);j=int(name.split('_j')[1]) if '_j' in name else 1
   row=mul(mul(mul(FF,pt),op[j]),dp[p])[p].real*math.factorial(p)
   if name.endswith('G_regular'):
    row-=math.factorial(p)*pka[p].real;outside-=vv*math.factorial(p)*kb[p].real
   num+=vv*row
  return outside,num
 def field(self,lo,hi):
  ta=ball(ar(lo),ar(hi));et=ball(arb(0),ar(1/self.T));w=ta*et;z=w/(1+ar(F(8,3))*w)
  outside,num=self.numerator(ta,et,w)
  # Positive Legendre-series bounds for the unknown integer l>=L.
  plo=arb(0)
  for co in reversed(self.pc):plo=plo*ar(lo)+ar(co)
  phi=(ar(F(4,3)*hi)).hypgeom_0f1(1);P=ball(plo.lower(),phi.upper())
  nc=arb(0);weight=arb(1)
  for k in range(self.K):
   if k:weight*=ar(F(4,3))*ta/(k*k)*(1-(k-1)*k*et)
   f=arb(0)
   for co in reversed(self.fc[k]):f=f*z+co
   nc+=weight*f
  # The omitted k>=64 crossing terms are nonpositive, exactly.
  source=self.model.source(lo,hi)
  # Mehler's integral error, uniform in every l>=L and source t.
  b=F(79,20)*hi;beta=b/self.T;assert beta<F(3,5)
  mehler=2*ar(b).sqrt()*((1+ar(1/(4*self.T))).sqrt()/(1-ar(beta)).sqrt()-1)+ar(4*beta/(6-10*beta))
  M=4+self.T/hi if hi else F(10)**100;e=1/M
  eps1=(8*e+2153*e*e+14088*e**3)/((1-48*e)*(1-F(20,3)*e));assert M>48
  es=ar(self.norm1)*(ar(eps1)+mehler)
  wm=hi/self.T;r=(F(79,20)+F(4,3))*wm/(1+F(8,3)*wm);assert r<1
  eh=ar(self.norm2*r*r/(1-r*r))
  return outside+eh+(num+nc-source+es)/P

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--lmin',type=int,default=1000);ap.add_argument('--lo',default='1');ap.add_argument('--hi',default='3');ap.add_argument('--prec',type=int,default=256);ap.add_argument('--test',action='store_true');args=ap.parse_args();ctx.prec=args.prec
 folder=ROOT/'output/lower_g0_upper_g2_20261006';w=folder/'blend_0.9.npz';tic=time.time()
 obj=Region(w,folder/'tail_blend90_joint_moments.pkl',folder/'matrix_blend90_source_norms.json',args.lmin)
 if args.test:
  for t in [F(0),F(1),F(2),F(57,25),F(3),F(5),F(10),F(21),F(40)]:print(str(t),str(obj.field(t,t)),flush=True)
  sys.exit()
 stack=[(F(args.lo),F(args.hi),0)];leaves=[];failed=[];evals=0
 while stack:
  lo,hi,depth=stack.pop();val=obj.field(lo,hi);evals+=1
  if val.is_finite() and val<0:leaves.append({'lo':str(lo),'hi':str(hi),'upper':str(val.upper())});continue
  if depth>=22:failed.append({'lo':str(lo),'hi':str(hi),'value':str(val)});break
  mid=(lo+hi)/2
  if depth>10 and obj.field(mid,mid)>0:failed.append({'point':str(mid),'value':str(obj.field(mid,mid))});break
  stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
  if evals%1000==0:print('evaluations',evals,'leaves',len(leaves),flush=True)
 result={'status':'PASS' if not failed and not stack else 'FAIL','witness_sha256':obj.sha,'spin_domain':f'all even l>={args.lmin}',
  'tau_definition':'l(l+1)/(mu-4)','tau_domain':[args.lo,args.hi],'precision_bits':args.prec,'leaves':leaves,'failed':failed,'evals':evals,'seconds':time.time()-tic,
  'claim':'Scoped all-spin sign component only. Does not cover complementary tau or spins, real-kernel support integrals, or certify a coupling bound.',
  'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 dest=folder/f'tail_blend90_finite_uv_L{args.lmin}_{args.lo}_{args.hi}_{args.prec}.json';dest.write_text(json.dumps(result,indent=2))
 print(json.dumps({k:v for k,v in result.items() if k!='leaves'},indent=2))
