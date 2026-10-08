"""Finite even-spin UV bridge using exact source correction polynomials."""
from pathlib import Path
from fractions import Fraction as F
import sys,pickle,json,math,hashlib,time,argparse
ROOT=Path(__file__).resolve().parents[3];sys.dont_write_bytecode=True
sys.path[:0]=[str(Path(__file__).parent),str(ROOT/'research/lower_g0_upper_g2_20261006')]
from contracted_forward import ContractedRegion,Region,shifted_eval,bounds,fq,fmpq_poly,arb,ar,pw
from tail_finite_uv_region import ctx
from finite_source_corrections import remainder
from matrix_crossing_series import coefficient_polynomial,bernstein_univariate
import sympy as sy
OUT=ROOT/'output/all_direction_certification_20261007/tail/g2_robust_stage1'
WITNESS=ROOT/'output/competitive_20261007/robust/g2_robust_stage1.npz'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upper(v):m,e=v.upper().man_exp();return F(int(m))*F(2)**int(e)
class FiniteUV:
 def __init__(self):
  self.f=ContractedRegion(WITNESS,OUT/'joint_moments.pkl',OUT/'source_norms.json',450)
  self.src=pickle.loads((OUT/'finite_source_N12_K80.pkl').read_bytes());assert self.src['witness_sha256']==sha(WITNESS)
  self.cd=pickle.loads((OUT/'compact_norms.pkl').read_bytes());assert self.cd['witness_sha256']==sha(WITNESS)
  assert self.cd['moments'][0]==self.f.d['Hs']
  assert self.src['moments'][0]==self.f.d['moments'][:81]
  assert all(self.src['moments'][1][k]==2*self.src['moments'][0][k]+F(3,2)*self.f.d['moments'][k+1]for k in range(81))
  pp=coefficient_polynomial(WITNESS);k,z=pp.gens;self.fc=[];self.D=0
  y=sy.Symbol('y');tail=sy.Poly(sy.expand(pp.as_expr().subs(k,y+64)),y)
  self.crossing_tail_coefficients=[max(bernstein_univariate(tail.nth(j),z,0,sy.Rational(3,8)))for j in range(tail.degree()+1)]
  assert all(v<=0 for v in self.crossing_tail_coefficients)
  for j in range(64):
   pol=sy.Poly(pp.as_expr().subs(k,j),z);co=[F(pol.nth(i))for i in range(max(0,pol.degree())+1)];self.D=max(self.D,len(co)-1);self.fc.append(co)
  self.data={}
 def prepare(self,l):
  if l in self.data:return self.data[l]
  T=l*(l+1);K=min(l,self.src['K']);aa=[F(1)]
  for k in range(1,l+1):aa.append(aa[-1]*F(T-(k-1)*k,T*k*k))
  pp=fmpq_poly([fq(a*F(4,3)**k)for k,a in enumerate(aa)])
  ss=[F(0)]*(K+self.src['N'])
  for j,row in enumerate(self.src['moments']):
   for k in range(K+1):ss[j+k]+=(-1)**k*aa[k]*row[k]/T**j
  sp=fmpq_poly(list(map(fq,ss)))
  x=fmpq_poly([0,1]);q=1+fq(F(8,3)/T)*x;cn=fmpq_poly([])
  for k in range(min(l+1,64)):
   for j,co in enumerate(self.fc[k]):cn+=fq(aa[k]*F(4,3)**k*co/T**j)*x**(k+j)*q**(self.D-j)
  self.data[l]=(T,pp,sp,cn);return self.data[l]
 def field(self,l,lo,hi):
  T,pp,sp,cn=self.prepare(l);assert lo>=0 and hi<=min(F(256),F(T,45))
  ta=ar((lo+hi)/2)+arb(0,ar((hi-lo)/2));et=ar(F(1,T));w=ta*et;q=1+ar(F(8,3))*w;zz=w/q
  outside,num=self.f.numerator(ta,et,w)
  P=arb(pp(fq(lo))).union(arb(pp(fq(hi))))
  cross=shifted_eval(cn,lo,hi)/pw(q,self.D)
  source=shifted_eval(sp,lo,hi)
  wm=hi/T;eps=remainder(self.src,wm)
  error=ar(self.f.norm1*eps)
  if l>self.src['K']:
   a=F(79,20)*hi;n=self.src['K']+1;rat=a/F((n+1)**2);assert rat<1
   leg=a**n/F(math.factorial(n)**2)/(1-rat)
   Rmax=(1+2*wm)*(1+F(8,3)*wm)**2/((1-44*wm)*(1-F(157,60)*wm))
   error+=ar(self.f.norm1*(Rmax+eps)*leg)
  # Exact signed K2 correction; do not retain the much looser O(w^2) error.
  N=12;s2=arb(0)
  for co in reversed(self.cd['moments'][:N]):s2=s2*zz*zz+ar(co)
  zmax=wm/(1+F(8,3)*wm);rr=F(317,60)*zmax;assert rr<1
  k2error=ar(self.f.norm2*rr**(2*N)/(1-rr*rr))
  return 1+outside+ar(self.f.d['Hs'])-s2+k2error+(num+cross-source+error)/P
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--spins',default='42,44,100,200,448');ap.add_argument('--precision',type=int,default=400);ap.add_argument('--screen',action='store_true');ap.add_argument('--kappa',default='1/1000');ap.add_argument('--replay');a=ap.parse_args()
 reference=json.loads(Path(a.replay).read_text())if a.replay else None
 ctx.prec=a.precision;tic=time.time();bank=FiniteUV();spins=reference['spins']if reference else (list(range(int(a.spins.split(':')[0]),int(a.spins.split(':')[1])+1,2))if ':'in a.spins else list(map(int,a.spins.split(','))));kappa=F(a.kappa);records=[];failed=[]
 assert all(l>=42 and l<=448 and l%2==0 for l in spins)
 if reference:
  assert reference['status']=='PASS_FINITE_UV'and F(reference['kappa_exact'])==kappa and reference['witness_sha256']==sha(WITNESS)
  for p,h in reference['sha256_manifest'].items():assert sha(ROOT/p)==h,p
 for l in spins:
  end=min(F(256),F(l*(l+1),45));t=time.time()
  if a.screen:
   for p in map(F,['0','1','2','9/4','5/2','3','4','8','16','32','64','128','256']):
    if p<=end:print(l,str(p),str(bank.field(l,p,p)),flush=True)
   continue
  cuts=sorted(set([F(0),end]+[F(x)for x in[1,2,3,4,8,16,32,64,128,256]if x<end]));stack=[(u,v,0)for u,v in zip(cuts,cuts[1:])];cells=[];count=0
  if reference:stack=[(F(c['a']),F(c['b']),0)for rr in reference['records']if rr['spin']==l for c in rr['cells']]
  while stack:
   lo,hi,depth=stack.pop();count+=1;v=bank.field(l,lo,hi)
   if v.is_finite()and v.upper()<-ar(kappa):cells.append(dict(a=str(lo),b=str(hi),upper_exact=str(upper(v))));continue
   if reference:failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v)));break
   mid=(lo+hi)/2
   if depth>=22:
    failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v)));break
   if depth>12:
    point=bank.field(l,mid,mid)
    if point.is_finite()and point.lower()>-ar(kappa):failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v),point=str(point)));break
   stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  cells.sort(key=lambda c:F(c['a']));passed=not stack and not failed
  if passed:
   assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==end and all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
  records.append(dict(spin=l,status='PASS'if passed else'FAILED_ENVELOPE',tau_end_exact=str(end),cells=cells,evaluations=count,seconds=time.time()-t))
  print(json.dumps(dict(spin=l,status=records[-1]['status'],cells=len(cells),evaluations=count,seconds=time.time()-t)),flush=True)
  if failed:break
 if a.screen:return
 files=[Path(__file__),Path(ContractedRegion.numerator.__code__.co_filename),Path(Region.field.__code__.co_filename),Path(coefficient_polynomial.__code__.co_filename),Path(bernstein_univariate.__code__.co_filename),OUT/'finite_source_N12_K80.pkl',OUT/'finite_source_N12_K80.json',OUT/'compact_norms.pkl',OUT/'prepared_components.json']
 if reference:files.append(Path(a.replay).resolve())
 manifest={}
 for p in[OUT/'prepared_components.json',OUT/'finite_source_N12_K80.json']:
  d=json.loads(p.read_text())
  for f,h in d['sha256_manifest'].items():assert sha(ROOT/f)==h;manifest[f]=h
 manifest.update({str(p.relative_to(ROOT)):sha(p)for p in files})
 out=dict(status='PASS_FINITE_UV'if not failed else'FAILED_ENVELOPE',witness_sha256=sha(WITNESS),spins=spins,kappa_exact=str(kappa),precision_bits=a.precision,records=records,failed=failed,seconds=time.time()-tic,
  scope='Listed spins, mu>=49, tau<=256. For ell>=108 the separately certified unbounded-tau component covers the rest.',acceptance_policy='Every listed spin and every closed cell has finite Arb value and strict directed upper below minus kappa.',crossing_omitted_kmin=64,crossing_tail_Bernstein_coefficients=list(map(str,bank.crossing_tail_coefficients)),reference_certificate_sha256=sha(a.replay)if reference else None,sha256_manifest=manifest)
 name='finite_uv_bridge_'+('_'.join(map(str,spins))if len(spins)<8 else f'{min(spins)}_{max(spins)}')+f'_{a.precision}.json'
 (OUT/name).write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
