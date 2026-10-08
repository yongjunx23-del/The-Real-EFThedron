"""Finite UV local coercivity, retaining paired directed upper bounds."""
from angle265_models import *
from finite_uv_bridge220 import ctx,ar,upper
import argparse,time

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--spins',default='42:90');ap.add_argument('--precision',type=int,default=400);ap.add_argument('--replay');a=ap.parse_args();ctx.prec=a.precision;tic=time.time();bank=Finite();manifest=source_manifest();reference=json.loads(Path(a.replay).read_text())if a.replay else None
 if reference:
  assert reference['status']=='PASS_ANGLE265_FINITE_UV'and reference['witness_sha256']==SHA and reference['objective_exact']==list(map(str,OBJECTIVE));spins=reference['spins']
 else:lo,hi=map(int,a.spins.split(':'));spins=list(range(lo,hi+1,2))
 assert all(42<=l<=448 and l%2==0 for l in spins);records=[];failed=[]
 for l in spins:
  t=time.time();end=min(F(256),F(l*(l+1),45));cuts=sorted({F(0),end}|{F(x)for x in['1','2','9/4','5/2','11/4','3','4','8','16','21','24','32','64','128','256']if F(x)<end})
  stack=[(x,y,0)for x,y in zip(cuts,cuts[1:])]
  if reference:stack=[(F(c['a']),F(c['b']),0)for r in reference['records']if r['spin']==l for c in r['cells']]
  cells=[];trials=0
  while stack:
   lo,hi,depth=stack.pop();trials+=1;v=bank.field(l,lo,hi);mid=(lo+hi)/2;point=None
   good=v.is_finite()and v.upper()<-ar(F(1,10))
   if good:cells.append(dict(a=str(lo),b=str(hi),upper_exact=str(upper(v))));continue
   if reference or depth>=25:failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v)));break
   if depth>14:
    if point is None:point=bank.field(l,mid,mid)
    if point.is_finite()and point.lower()>-ar(F(1,10)):failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v),midpoint=str(point)));break
   stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  cells.sort(key=lambda c:F(c['a']));passed=not failed and not stack
  if passed:assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==end and all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
  records.append(dict(spin=l,status='PASS'if passed else'FAILED',tau_end_exact=str(end),cells=cells,trials=trials));print(json.dumps(dict(spin=l,status=records[-1]['status'],cells=len(cells),trials=trials,seconds=time.time()-t)),flush=True)
  if failed:break
 files=[Path(__file__),Path(__file__).with_name('angle265_models.py'),Path(fv.__file__),Path(fv.ContractedRegion.numerator.__code__.co_filename),Path(fv.coefficient_polynomial.__code__.co_filename)]
 if reference:files.append(Path(a.replay).resolve())
 manifest.update({str(p.relative_to(ROOT)):sha(p)for p in files})
 for p,h in manifest.items():assert sha(ROOT/p)==h,p
 out=dict(status='PASS_ANGLE265_FINITE_UV'if not failed else'FAILED_ENVELOPE',witness_sha256=SHA,objective_exact=list(map(str,OBJECTIVE)),spins=spins,kappa_exact='1/10',precision_bits=a.precision,records=records,failed=failed,crossing_tail_coefficients=list(map(str,bank.crossing_tail_coefficients)),reference_certificate_sha256=sha(a.replay)if reference else None,acceptance_policy='Every interval finite and strict directed upper<-1/10; no midpoint-relative refinement.',sha256_manifest=manifest,seconds=time.time()-tic)
 (OUT/f'finite_{min(spins)}_{max(spins)}_{a.precision}_k1_10.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
