"""Uniform finite-eta source contraction, initially for the frozen left witness.

The retained Legendre coefficients and R1 corrections are exact bivariate
polynomials in (tau,eta). Error bounds are uniform in every integer ell>=L.
"""
from finite_uv_bridge220 import *
LEFT=ROOT/'output/full_direction_20261007/left/angle180_final.npz'
FOLDER=ROOT/'output/all_direction_certification_20261007/tail/left/refined'

def add_matrix(mat,h,k,value):mat[h][k]+=value
def bivariate_value(rows,lo,hi,eta):
 out=arb(0)
 for pol in reversed(rows):out=out*eta+shifted_eval(pol,lo,hi)
 return out
class UniformSource:
 def __init__(self,L=450):
  self.L=L;self.T=L*(L+1);assert L>=108
  self.f=ContractedRegion(LEFT,FOLDER/'joint_moments.pkl',FOLDER/'source_norms.json',L)
  self.src=pickle.loads((FOLDER/'finite_source_N12_K80.pkl').read_bytes());self.cd=pickle.loads((FOLDER/'compact_norms.pkl').read_bytes())
  assert self.src['witness_sha256']==sha(LEFT)==self.cd['witness_sha256']
  assert self.cd['moments'][0]==self.f.d['Hs']and self.src['moments'][0]==self.f.d['moments'][:81]
  N=self.src['N'];K=self.src['K'];assert L>=K
  aa=[[F(1)]]
  for k in range(1,K+1):
   prev=aa[-1];row=[F(0)]*(len(prev)+1)
   for j,a in enumerate(prev):row[j]+=a/(k*k);row[j+1]-=a*F((k-1)*k,k*k)
   while row[-1]==0:row.pop()
   aa.append(row)
  mat=[[F(0)]*(K+N)for _ in range(K+N)]
  for j,row in enumerate(self.src['moments']):
   for k,m in enumerate(row):
    for h,a in enumerate(aa[k]):add_matrix(mat,h+j,k+j,(-1)**k*a*m)
  self.source=[fmpq_poly(list(map(fq,row)))for row in mat]
  pp=coefficient_polynomial(LEFT);k,z=pp.gens;y=sy.Symbol('y');tail=sy.Poly(sy.expand(pp.as_expr().subs(k,y+64)),y)
  self.tail_coeff=[max(bernstein_univariate(tail.nth(j),z,0,sy.Rational(3,8)))for j in range(tail.degree()+1)]
  assert all(v<=0 for v in self.tail_coeff)
  fc=[];self.D=0
  for kk in range(64):
   p=sy.Poly(pp.as_expr().subs(k,kk),z);co=[F(p.nth(j))for j in range(max(0,p.degree())+1)];fc.append(co);self.D=max(self.D,len(co)-1)
  mat=[[F(0)]*(64+self.D)for _ in range(64+self.D)]
  for k,co in enumerate(fc):
   for j,c in enumerate(co):
    for r in range(self.D-j+1):
     v=F(4,3)**k*c*math.comb(self.D-j,r)*F(8,3)**r
     for h,a in enumerate(aa[k]):add_matrix(mat,h+j+r,k+j+r,v*a)
  self.cross=[fmpq_poly(list(map(fq,row)))for row in mat]
 def p_bounds(self,lo,hi,elo,ehi):
  def prefix(t,e):
   val=arb(1);term=arb(1)
   for k in range(1,self.src['K']+1):term*=ar(F(4,3)*t*(1-(k-1)*k*e)/F(k*k));val+=term
   return val
  n=self.src['K']+1;a=F(4,3)*hi;ratio=a/F((n+1)**2);assert ratio<1
  tail=a**n/F(math.factorial(n)**2)/(1-ratio)
  lower=prefix(lo,ehi);upper=prefix(hi,elo)+ar(tail)
  return lower.union(upper)
 def field(self,lo,hi,elo,ehi):
  assert 0<=lo<=hi<=256 and 0<=elo<=ehi<=F(1,self.T)
  ta=ar((lo+hi)/2)+arb(0,ar((hi-lo)/2));et=ar((elo+ehi)/2)+arb(0,ar((ehi-elo)/2));w=ta*et;q=1+ar(F(8,3))*w;zz=w/q
  outside,num=self.f.numerator(ta,et,w)
  cross=bivariate_value(self.cross,lo,hi,et)/pw(q,self.D);source=bivariate_value(self.source,lo,hi,et)
  wm=hi*ehi;eps=remainder(self.src,wm);error=ar(self.f.norm1*eps)
  n=self.src['K']+1;a=F(79,20)*hi;ratio=a/F((n+1)**2);assert ratio<1
  leg=a**n/F(math.factorial(n)**2)/(1-ratio)
  Rmax=(1+2*wm)*(1+F(8,3)*wm)**2/((1-44*wm)*(1-F(157,60)*wm))
  error+=ar(self.f.norm1*(Rmax+eps)*leg)
  s2=arb(0)
  for co in reversed(self.cd['moments'][:12]):s2=s2*zz*zz+ar(co)
  zmax=wm/(1+F(8,3)*wm);rr=F(317,60)*zmax;k2error=ar(self.f.norm2*rr**24/(1-rr*rr));assert rr<1
  # The left objective b is exactly zero. Alpha already includes a=-1.
  return outside+ar(self.f.d['Hs'])-s2+k2error+(num+cross-source+error)/self.p_bounds(lo,hi,elo,ehi)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lmin',type=int,default=450);ap.add_argument('--precision',type=int,default=400);ap.add_argument('--eta-panels',type=int,default=1);a=ap.parse_args();ctx.prec=a.precision;tic=time.time();bank=UniformSource(a.lmin);records=[]
 for tau in map(F,['0','1','2','9/4','5/2','3','5','10','21','40','64','128','256']):
  vals=[bank.field(tau,tau,F(i,a.eta_panels*bank.T),F(i+1,a.eta_panels*bank.T))for i in range(a.eta_panels)]
  uu=max(upper(v)for v in vals if v.is_finite())if all(v.is_finite()for v in vals)else None
  records.append(dict(tau_exact=str(tau),upper_exact=str(uu)if uu is not None else None));print(str(tau),float(uu)if uu is not None else'nonfinite',flush=True)
 out=dict(status='SCOPED_UNIFORM_SOURCE_SCREEN',witness_sha256=sha(LEFT),lmin=a.lmin,eta_panels=a.eta_panels,precision_bits=a.precision,records=records,seconds=time.time()-tic,proof_status='New uniform bivariate extension requires independent review; listed tau points only.')
 (FOLDER/f'uniform_source_L{a.lmin}_panels{a.eta_panels}_{a.precision}_screen.json').write_text(json.dumps(out,indent=2));print('seconds',time.time()-tic)
if __name__=='__main__':main()
