"""Uniform old150 b6 envelope for all even ell>=42 and 0<=tau<=8.

Exact source/crossing bivariate contractions use K40 so every retained
Legendre coefficient is nonnegative over eta in[0,1/(42*43)]. The remaining
source terms use the same absolute Bessel-series remainder as the reviewed
left proof. Existing frozen Jensen packets handle tau>=8.
"""
from uniform_source_uv import *
from prepare_b6_sources import W,FOLDER

class B6Uniform:
 def __init__(self):
  self.T=42*43;self.K=40;self.b=F(6);self.delta=F(107,5)
  self.f=ContractedRegion(W,FOLDER/'joint_moments.pkl',FOLDER/'source_norms.json',450)
  self.src=pickle.loads((FOLDER/'finite_source_N12_K40.pkl').read_bytes());self.cd=pickle.loads((FOLDER/'compact_norms.pkl').read_bytes())
  assert self.src['witness_sha256']==sha(W)==self.cd['witness_sha256']
  assert self.src['moments'][0]==self.f.d['moments'][:41]and self.cd['moments'][0]==self.f.d['Hs']
  assert [n for n in self.f.lam if n.startswith('X')]==['X4']
  ratio=(self.f.lam['X4']+self.delta)/self.f.lam['X4'];assert 0<ratio<1
  aa=[[F(1)]]
  for k in range(1,self.K+1):
   row=[F(0)]*(len(aa[-1])+1)
   for h,v in enumerate(aa[-1]):row[h]+=v/(k*k);row[h+1]-=v*F((k-1)*k,k*k)
   while row[-1]==0:row.pop()
   aa.append(row)
  N=self.src['N'];mat=[[F(0)]*(self.K+N)for _ in range(self.K+N)]
  for j,row in enumerate(self.src['moments']):
   for k,m in enumerate(row):
    for h,v in enumerate(aa[k]):mat[h+j][k+j]+=(-1)**k*v*m
  self.source=[fmpq_poly(list(map(fq,row)))for row in mat]
  pp=coefficient_polynomial(W);k,z=pp.gens;pp=sy.Poly(pp.as_expr()*sy.Rational(ratio),k,z)
  y=sy.Symbol('y');tail=sy.Poly(sy.expand(pp.as_expr().subs(k,y+41)),y)
  self.tail_coeff=[max(bernstein_univariate(tail.nth(j),z,0,sy.Rational(3,8)))for j in range(tail.degree()+1)]
  assert all(v<=0 for v in self.tail_coeff)
  fc=[];self.D=0
  for kk in range(41):
   p=sy.Poly(pp.as_expr().subs(k,kk),z);co=[F(p.nth(j))for j in range(max(0,p.degree())+1)];fc.append(co);self.D=max(self.D,len(co)-1)
  mat=[[F(0)]*(41+self.D)for _ in range(41+self.D)]
  for k,co in enumerate(fc):
   for j,c in enumerate(co):
    for r in range(self.D-j+1):
     v=F(4,3)**k*c*math.comb(self.D-j,r)*F(8,3)**r
     for h,a in enumerate(aa[k]):mat[h+j+r][k+j+r]+=v*a
  self.cross=[fmpq_poly(list(map(fq,row)))for row in mat]
 def p_bounds(self,lo,hi,elo,ehi):
  def prefix(t,e):
   val=arb(1);term=arb(1)
   for k in range(1,41):term*=ar(F(4,3)*t*(1-(k-1)*k*e)/F(k*k));val+=term
   return val
  a=F(4,3)*hi;tail=a**41/F(math.factorial(41)**2)/(1-a/F(42**2))
  return prefix(lo,ehi).union(prefix(hi,elo)+ar(tail))
 def field(self,lo,hi,elo,ehi):
  assert 0<=lo<=hi<=8 and 0<=elo<=ehi<=F(1,self.T)
  ta=ar((lo+hi)/2)+arb(0,ar((hi-lo)/2));et=ar((elo+ehi)/2)+arb(0,ar((ehi-elo)/2));w=ta*et;q=1+ar(F(8,3))*w;zz=w/q
  outside,num=self.f.numerator(ta,et,w);cross=bivariate_value(self.cross,lo,hi,et)/pw(q,self.D);source=bivariate_value(self.source,lo,hi,et)
  wm=hi*ehi;eps=remainder(self.src,wm);a=F(79,20)*hi;leg=a**41/F(math.factorial(41)**2)/(1-a/F(42**2))
  Rmax=(1+2*wm)*(1+F(8,3)*wm)**2/((1-44*wm)*(1-F(157,60)*wm))
  error=ar(self.f.norm1*(eps+(Rmax+eps)*leg));s2=arb(0)
  for co in reversed(self.cd['moments'][:12]):s2=s2*zz*zz+ar(co)
  rr=F(317,60)*wm/(1+F(8,3)*wm);assert rr<1
  k2error=ar(self.f.norm2*rr**24/(1-rr*rr))
  return ar(self.b)+outside+ar(self.f.d['Hs'])-s2+k2error+(num+cross-source+error)/self.p_bounds(lo,hi,elo,ehi)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=320);ap.add_argument('--screen',action='store_true');ap.add_argument('--panels',type=int,default=8);ap.add_argument('--replay');a=ap.parse_args();ctx.prec=a.precision;tic=time.time();bank=B6Uniform();records=[];reference=json.loads(Path(a.replay).read_text())if a.replay else None
 if a.screen:
  for t in map(F,['0','1/2','1','2','3','4','6','8']):
   vals=[bank.field(t,t,F(i,a.panels*bank.T),F(i+1,a.panels*bank.T))for i in range(a.panels)]
   print(str(t),float(max(upper(v)for v in vals)),flush=True)
  return
 for i in range(a.panels):
  elo,ehi=F(i,a.panels*bank.T),F(i+1,a.panels*bank.T)
  stack=[(F(0),F(8),0)]if reference is None else [(F(c['a']),F(c['b']),0)for c in reference['records'][i]['cells']]
  cells=[];trials=0
  while stack:
   lo,hi,depth=stack.pop();trials+=1;v=bank.field(lo,hi,elo,ehi)
   if v.is_finite()and v.upper()<-1:cells.append(dict(a=str(lo),b=str(hi),upper_exact=str(upper(v))));continue
   assert reference is None and depth<22,(i,lo,hi,str(v))
   mid=(lo+hi)/2;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  cells.sort(key=lambda c:F(c['a']));assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==8 and all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
  records.append(dict(panel=i,eta_lo_exact=str(elo),eta_hi_exact=str(ehi),cells=cells));print('panel',i,'cells',len(cells),'seconds',time.time()-tic,flush=True)
 prep=json.loads((FOLDER/'prepared_components.json').read_text());manifest=dict(prep['sha256_manifest'])
 for p in[Path(__file__),Path(ContractedRegion.numerator.__code__.co_filename),Path(Region.numerator.__code__.co_filename),Path(coefficient_polynomial.__code__.co_filename),Path(bivariate_value.__code__.co_filename),FOLDER/'prepared_components.json']+([Path(a.replay).resolve()]if reference else[]):manifest[str(p.relative_to(ROOT))]=sha(p)
 for p,h in manifest.items():assert sha(ROOT/p)==h,p
 out=dict(status='PASS_B6_UNIFORM_UV',base_witness_sha256=sha(W),objective_exact=['1','6'],delta_X4_exact='107/5',lmin=42,tau_domain=['0','8'],mu_lower_bound_exact=str(4+F(bank.T,8)),kappa_exact='1',eta_panels=a.panels,precision_bits=a.precision,records=records,crossing_omitted_kmin=41,crossing_tail_coefficients=list(map(str,bank.tail_coeff)),reference_certificate_sha256=sha(a.replay)if reference else None,sha256_manifest=manifest,seconds=time.time()-tic)
 (FOLDER/f'uniform_uv_{a.precision}.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in['records','sha256_manifest']},indent=2))

if __name__=='__main__':main()
