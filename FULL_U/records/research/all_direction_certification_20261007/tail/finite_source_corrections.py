"""Exact finite-w Roy/strip source moments and rational remainder.

No Mehler approximation. R1 is expanded to finite order; the exact rational
remainder bound remains valid even outside the power series convergence disk.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,re,math,pickle,json,time,hashlib,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research/tail_certification'),str(ROOT/'research/lower_g0_upper_g2_20261006')]
import numpy as np
import sympy as sy
from flint import fmpq
from certify_uv import bernstein,split_box
WITNESS=ROOT/'output/competitive_20261007/robust/g2_robust_stage1.npz'
OUT=ROOT/'output/all_direction_certification_20261007/tail/g2_robust_stage1'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ff(q):q=F(q);return fmpq(q.numerator,q.denominator)
def frac(q):return F(int(q.numerator),int(q.denominator))
def specs():
 st=np.load(WITNESS);x,z=sy.symbols('x z');c=sy.Rational(4,3);roy=0;strips={}
 for name,v in zip(st['names'],st['multipliers']):
  v=sy.Rational(float(v));mm=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',str(name))
  if mm:
   assert sy.Rational(mm[1])==sy.Rational(119,10)
   roy+=2*v*sy.legendre(int(mm[3]),2*x-1)*sy.legendre(int(mm[2]),z);continue
  mm=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',str(name))
  if mm:
   A,B=F(mm[1]),F(mm[2]);j,n=int(mm[4]),int(mm[5]);strips[A,B]=strips.get((A,B),0)+v*(1-(2*x-1)**2)*sy.legendre(n,2*x-1)*sy.legendre(j,2*z-1)
 result=[('roy',roy,4+sy.Rational(79,10)*x,-sy.Rational(79,20)*x*(1-z),True,F(79,20))]
 for(A,B),poly in strips.items():result.append((f'strip{A}_{B}',poly,sy.Rational(A)+sy.Rational(B-A)*x,-c*(1-z),False,F(4,3)))
 return x,z,result
def build(N=12,K=80):
 tic=time.time();x,z,groups=specs();s,t=sy.symbols('s t');c=sy.Rational(4,3)
 # denominator=(1+(4-s)w)(1+(s+t)w)(1+(c+t)w)
 ds=[sy.S.One,sy.expand(4-s+s+t+c+t),sy.expand((4-s)*(s+t)+(4-s)*(c+t)+(s+t)*(c+t)),sy.expand((4-s)*(s+t)*(c+t))]
 ns=[sy.S.One,sy.Rational(22,3)+t/2,sy.Rational(160,9)+sy.Rational(8,3)*t,sy.Rational(128,9)+sy.Rational(32,9)*t]
 rs=[]
 for j in range(N+1):rs.append(sy.expand((ns[j]if j<4 else 0)-sum(ds[h]*rs[j-h]for h in range(1,min(j,3)+1))))
 rem=[rs[N],sy.expand(-ds[2]*rs[N-1]-ds[3]*rs[N-2]),sy.expand(-ds[3]*rs[N-1])]
 moments=[[fmpq(0)for k in range(K+1)]for j in range(N)];bounds=[]
 for label,W,ss,tt,isroy,scale in groups:
  kval=2*(ss-c)*(c-(4-ss-tt));base=sy.expand(W*kval)
  for j in range(N):
   pol=sy.Poly(sy.expand(base*rs[j].subs({s:ss,t:tt})),x,z);terms=[(a,b,ff(co))for(a,b),co in pol.terms()]
   bmax=max(b for a,b,co in terms);betas=[]
   for b in range(bmax+1):
    bb=[fmpq(1,b+1)]
    for k in range(1,K+1):bb.append(bb[-1]*k/(b+k+1))
    betas.append(bb)
   qq=fmpq(1);qs=ff(scale)
   for k in range(K+1):
    moments[j][k]+=qq*sum((co*betas[b][k]/(a+k+1 if isroy else a+1)for a,b,co in terms),fmpq(0));qq*=qs
   print(json.dumps(dict(group=label,row=j,terms=len(terms),seconds=time.time()-tic)),flush=True)
  # Uniform absolute rational-remainder coefficient bounds, exact Bernstein.
  rb=[]
  for expr in rem:
   bs=bernstein(sy.expand(expr.subs({s:ss,t:tt})),x,z)
   rb.append(max(abs(v)for box,*_ in split_box(bs,depth=6)for row in box for v in row))
  bounds.append(dict(source=label,coefficient_abs_upper_exact=list(map(str,rb))))
 allbounds=[max(F(g['coefficient_abs_upper_exact'][j])for g in bounds)for j in range(3)]
 return dict(witness_sha256=sha(WITNESS),N=N,K=K,moments=[[frac(q)for q in row]for row in moments],
  remainder_coefficients=allbounds,remainder_by_source=bounds,seconds=time.time()-tic,
  definition='M[j,k]=integral W*K1_infinity*r_j(s,t)*(-t)^k, with R1=sum_(j<N)r_j w^j + w^N(A+B w+C w²)/den.',
  denominator_lower='(1-44w)*(1-157w/60), valid0<=w<1/44; middle factor1+(s+t)w>=1.',
  proof='Exact rational beta moments. Remainder is a polynomial identity, so no convergent-series assumption is needed.',
  builder_sha256=sha(__file__))
def remainder(data,w):
 assert w>=0 and w<F(1,44)
 A,B,C=data['remainder_coefficients']
 return w**data['N']*(A+B*w+C*w*w)/((1-44*w)*(1-F(157,60)*w))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--order',type=int,default=12);ap.add_argument('--kmax',type=int,default=80);a=ap.parse_args()
 d=build(a.order,a.kmax);p=OUT/f'finite_source_N{a.order}_K{a.kmax}.pkl';p.write_bytes(pickle.dumps(d))
 report={k:v for k,v in d.items()if k not in ['moments','remainder_coefficients']};report['remainder_coefficients_exact']=list(map(str,d['remainder_coefficients']))
 report['sha256_manifest']={str(q.relative_to(ROOT)):sha(q)for q in[WITNESS,p,Path(__file__),ROOT/'research/tail_certification/certify_uv.py']}
 report['remainder_examples']={str(mu):float(remainder(d,F(1,mu-4)))for mu in[49,64,100,200,1000]}
 p.with_suffix('.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
