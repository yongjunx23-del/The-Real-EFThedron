"""Coarse all-spin compact-energy coercivity using fifth-degree crossing.

Unlike the old150 proof, only W_Roy>=0 on t<=-4/3 is required. The
dangerous nonphysical angular sector is discarded after the PV subtraction.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,re,math,pickle,json,hashlib,time,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/tail_certification'),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research')]
from certify_uv import ar,arb,ctx,leg_coeffs,poly
from compact_source_data import compute
from matrix_crossing_series import coefficient_polynomial,bernstein_univariate
import sympy as sy
OUT=ROOT/'output/all_direction_certification_20261007/tail/g2_robust_stage1'
WITNESS=ROOT/'output/competitive_20261007/robust/g2_robust_stage1.npz'
C=F(4,3)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def up(v):assert v.is_finite();m,e=v.upper().man_exp();return F(int(m))*F(2)**int(e)
def lo(v):assert v.is_finite();m,e=v.lower().man_exp();return F(int(m))*F(2)**int(e)
def mean(L,delta):
 cs=leg_coeffs(L);aa=ar(C*L*(L+1)/delta)
 return poly([k*v for k,v in enumerate(cs)],aa)/poly(cs,aa)
def forward_bounds(d,omega=F(4,3),complexFF=50*F(16,9)**3):
 lam=d['lam'];alpha=d['alpha'];FF=F(2048,27)
 co=[FF*abs(alpha)+F(64,27)*abs(alpha),F(0),F(0),F(0)];outside=F(128,27)*abs(alpha)
 mka=F(3136,81);mkb=2*mka
 for name,v in lam.items():
  mm=re.fullmatch(r'H(\d+)',name)
  if mm:co[0]+=abs(v)*FF;continue
  mm=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',name)
  if mm:
   kind,r,c=mm[1],F(mm[2]),F(mm[3]);dm=1+2*r*c+r*r;assert dm>0 and 0<=r<1
   top=(abs(c)+r+r*(1+r*abs(c)))if kind=='cos'else(1+r*r)
   co[0]+=abs(v)*FF*2*r*top/(dm*(1-r)**2);continue
  if not name.startswith('Dt'):continue
  p=int(name[2]);j=int(name.split('_j')[1])if'_j'in name else 1
  mm=complexFF*omega**j*(1+omega)**(2*p-2)
  if name.endswith('G_regular'):mm+=mka;outside+=abs(v)*math.factorial(p)*mkb
  for k in range(p+1):co[k]+=abs(v)*math.factorial(p)*mm/F(math.factorial(k)**2)
 return co,outside
def ratio_bound(co,L,D,K=None):
 T=L*(L+1)
 if K is None:
  term=F(1);records=[]
  for j in range(1,L+1):
   term*=C*F(T-(j-1)*j,D*j*j)
   if j>=len(co):records.append((term,j))
  K=max(records)[1]
 ak=C**K/F(math.factorial(K)**2)
 for j in range(K):ak*=T-j*(j+1)
 assert K>=len(co) and L>=K
 return sum((v*T**k*D**(K-k)/ak for k,v in enumerate(co)),F(0))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lmin',type=int,default=450);ap.add_argument('--precision',type=int,default=320);a=ap.parse_args();ctx.prec=a.precision;tic=time.time();L=a.lmin
 prep=json.loads((OUT/'prepared_components.json').read_text())
 for p,h in prep['sha256_manifest'].items():assert sha(ROOT/p)==h,p
 d=pickle.loads((OUT/'joint_moments.pkl').read_bytes())['data'];signpath=ROOT/'output/competitive_20261007/robust/stage1_audit/roy_exact.json';sg=json.loads(signpath.read_text())
 assert sg['sha256']==sha(WITNESS) and sg['passed'] and F(sg['metadata']['ymin'])==2*C and F(sg['W_lower_rational'])>0
 cache=OUT/'compact_norms.pkl'
 if cache.exists():cd=pickle.loads(cache.read_bytes());assert cd['witness_sha256']==sha(WITNESS)
 else:cd=compute(WITNESS,depth=10,N=40);cache.write_bytes(pickle.dumps(cd))
 pp=coefficient_polynomial(WITNESS);k,z=pp.gens;y=sy.Symbol('y');p=sy.Poly(sy.expand(pp.as_expr().subs(k,y+10)),y)
 coeff=[max(bernstein_univariate(p.nth(j),z,0,sy.Rational(3,8)))for j in range(p.degree()+1)]
 assert all(v<=0 for v in coeff) and coeff[5]<0;power5=-F(coeff[5])
 M=max([F(0)]+[F(max(bernstein_univariate(pp.as_expr().subs(k,j),z,0,sy.Rational(3,8))))for j in range(10)])
 co,base=forward_bounds(d);R,S=cd['roy_abs'],cd['strip_abs'];qmax=F(20,3)
 m=mean(L,F(4));rr=ar(lo(m))/L;assert rr*L>10
 slope=qmax**3*(R*F(3,2)/C+S/C)
 const=1+base+ratio_bound(co,L,F(4))+qmax**3*(R*F(13,8)+S*F(15,8))+M
 low=ar(const)+ar(slope)*L-ar(power5)*(rr*L-10)**5
 derivative=ar(slope)-5*ar(power5)*rr*(rr*L-10)**4
 assert low.is_finite() and derivative.is_finite() and low.upper()<-1 and derivative.upper()<0
 # At mu>=8 every source angular argument is physical. Separate K2 and
 # bound the remaining Roy PV and tapered-strip Hilbert terms over P_l.
 qmin=F(20,3);qmax=F(143,3);dmin=F(83,60)
 norm2=F(json.loads((OUT/'source_norms.json').read_text())['norm2_rational'])
 outside_source=norm2/(1-(F(317,60)/qmin)**2)
 extra0=qmax**3*(R*(F(1,8)+1/qmin+1/dmin)+sum((2*v/(B-A)for(A,B),v in cd['strip_D'].items()),F(0))+S*(F(1,8)+1/qmin+F(1,4)))
 # Roy PV contributes T/(2delta); denominator delta>=4.
 augmented=co.copy();augmented[0]+=extra0;augmented[1]+=qmax**3*R/2
 mm=mean(L,F(45));assert mm>10
 middle=ar(1+base+outside_source+M+ratio_bound(augmented,L,F(45)))-ar(power5)*(mm-10)**5
 assert middle.is_finite() and middle.upper()<-1
 files=[Path(__file__),cache,signpath,Path(compute.__code__.co_filename),ROOT/'research/crossing_kernels.py']
 hashes=dict(prep['sha256_manifest']);hashes.update({str(p.relative_to(ROOT)):sha(p)for p in files})
 out=dict(status='PASS_COMPACT_HIGH_SPIN',witness_sha256=sha(WITNESS),lmin=L,mu_domain=['4','49'],kappa_exact='1',precision_bits=a.precision,
  low_mu4_8_upper_exact=str(up(low)),low_spin_derivative_upper_exact=str(up(derivative)),middle_mu8_49_upper_exact=str(up(middle)),
  crossing_power5_exact=str(power5),crossing_low_k_cap_exact=str(M),low_mean_lower_exact=str(lo(m)),middle_mean_lower_exact=str(lo(mm)),
  forward_coefficients_exact=list(map(str,co)),forward_outside_exact=str(base),roy_abs_exact=str(R),strip_abs_exact=str(S),
  low_forward_ratio_exact=str(ratio_bound(co,L,F(4))),middle_residual_ratio_exact=str(ratio_bound(augmented,L,F(45))),
  claim='Every even ell>=lmin and every4<mu<=49 has V/p<-1. Requires the exact Roy-source sign and the stated dispersion/PV assumptions.',
  new_proof_review_status='Sharper omega<4/3 by analytic logarithmic derivative on complex transfer disk; independent review requested.',seconds=time.time()-tic,sha256_manifest=hashes)
 (OUT/f'compact_highspin_sharp_L{L}_{a.precision}.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['sha256_manifest','forward_coefficients_exact','roy_abs_exact','strip_abs_exact','low_forward_ratio_exact','middle_residual_ratio_exact']},indent=2))
if __name__=='__main__':main()
