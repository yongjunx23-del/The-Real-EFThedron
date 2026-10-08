"""Finite left UV bridge retaining useful local coercivity margins."""
import finite_uv_bridge220 as uv
from uniform_source_uv import LEFT,FOLDER
from finite_uv_bridge220 import *

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--spins',default='42,44,100,448');ap.add_argument('--precision',type=int,default=400);ap.add_argument('--kappa',default='1/10000');ap.add_argument('--replay');a=ap.parse_args();ctx.prec=a.precision;tic=time.time()
 uv.OUT=FOLDER;uv.WITNESS=LEFT;bank=uv.FiniteUV();kap=F(a.kappa);reference=json.loads(Path(a.replay).read_text())if a.replay else None
 if reference:
  assert reference['status']=='PASS_LEFT_FINITE_UV'and reference['witness_sha256']==sha(LEFT)and F(reference['kappa_exact'])==kap
  for p,h in reference['sha256_manifest'].items():assert sha(ROOT/p)==h,p
  spins=reference['spins']
 else:spins=list(range(int(a.spins.split(':')[0]),int(a.spins.split(':')[1])+1,2))if':'in a.spins else list(map(int,a.spins.split(',')))
 assert all(l>=42 and l<=448 and l%2==0 for l in spins);records=[];failed=[]
 for l in spins:
  t=time.time();end=min(F(256),F(l*(l+1),45));cuts=sorted(set([F(0),end]+[F(x)for x in['1','2','11/5','23/10','12/5','5/2','3','4','8','16','21','24','32','64','128','256']if F(x)<end]));stack=[(x,y,0)for x,y in zip(cuts,cuts[1:])]
  if reference:stack=[(F(c['a']),F(c['b']),0)for r in reference['records']if r['spin']==l for c in r['cells']]
  cells=[];trials=0
  while stack:
   lo,hi,depth=stack.pop();trials+=1;v=bank.field(l,lo,hi)-1;mid=(lo+hi)/2;point=None
   good=v.is_finite()and v.upper()<-ar(kap)
   if good and not reference:
    point=bank.field(l,mid,mid)-1
    good=point.is_finite()and v.upper()<=point.upper()/4
   if good:cells.append(dict(a=str(lo),b=str(hi),upper_exact=str(upper(v))));continue
   if reference or depth>=24:failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v)));break
   if depth>=12:
    if point is None:point=bank.field(l,mid,mid)-1
    if point.is_finite()and point.lower()>-ar(kap):failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v),point=str(point)));break
   stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  cells.sort(key=lambda c:F(c['a']));passed=not failed and not stack
  if passed:
   assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==end and all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
  records.append(dict(spin=l,status='PASS'if passed else'FAILED',tau_end_exact=str(end),cells=cells,evaluations=trials,seconds=time.time()-t));print(json.dumps({k:v for k,v in records[-1].items()if k!='cells'}|dict(cell_count=len(cells))),flush=True)
  if failed:break
 prep=json.loads((FOLDER/'prepared_components.json').read_text());manifest=dict(prep['sha256_manifest'])
 files=[Path(__file__),Path(uv.__file__),Path(ContractedRegion.numerator.__code__.co_filename),Path(Region.field.__code__.co_filename),Path(coefficient_polynomial.__code__.co_filename),FOLDER/'prepared_components.json']
 if reference:files.append(Path(a.replay).resolve())
 for p in files:manifest[str(p.relative_to(ROOT))]=sha(p)
 for p,h in manifest.items():assert sha(ROOT/p)==h,p
 out=dict(status='PASS_LEFT_FINITE_UV'if not failed else'FAILED_ENVELOPE',witness_sha256=sha(LEFT),objective_exact=['-1','0'],spins=spins,kappa_exact=str(kap),precision_bits=a.precision,records=records,failed=failed,reference_certificate_sha256=sha(a.replay)if reference else None,acceptance_policy='Every accepted cell finite and directed upper<-kappa. Producer also refines upper<=midpoint_upper/4 to retain useful local margins; replay tests sign and stores its fresh upper.',crossing_omitted_kmin=64,crossing_tail_Bernstein_coefficients=list(map(str,bank.crossing_tail_coefficients)),seconds=time.time()-tic,sha256_manifest=manifest)
 name=f'left_finite_uv_{min(spins)}_{max(spins)}_{a.precision}.json';(FOLDER/name).write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
