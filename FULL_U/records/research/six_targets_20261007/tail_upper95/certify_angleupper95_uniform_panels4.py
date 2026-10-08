"""Uniform infinite-spin UV proof on rational eta panels, with exact direction."""
from angleupper95_models import *
from finite_uv_bridge220 import ctx,ar,upper
import argparse,time

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=400);ap.add_argument('--group',default='0:3');ap.add_argument('--replay');a=ap.parse_args();ctx.prec=a.precision;tic=time.time();bank=Uniform(450);manifest=source_manifest();reference=json.loads(Path(a.replay).read_text())if a.replay else None
 indices=list(range(*[int(a.group.split(':')[0]),int(a.group.split(':')[1])+1]))
 if reference:
  assert reference['status']=='PASS_ANGLEupper95_UNIFORM_UV'and reference['witness_sha256']==SHA and reference['objective_exact']==list(map(str,OBJECTIVE))
  indices=reference['panel_indices']
 records=[];failed=[]
 for i in indices:
  elo,ehi=F(i,4*bank.T),F(i+1,4*bank.T);stack=[(F(0),F(256),0)]
  if reference:stack=[(F(c['a']),F(c['b']),0)for r in reference['records']if r['panel']==i for c in r['cells']]
  cells=[];trials=0
  while stack:
   lo,hi,depth=stack.pop();trials+=1;v=bank.field(lo,hi,elo,ehi)
   if v.is_finite()and v.upper()<-ar(F(1,200)):cells.append(dict(a=str(lo),b=str(hi),upper_exact=str(upper(v))));continue
   mid=(lo+hi)/2
   if reference or depth>=25:
    failed.append(dict(panel=i,a=str(lo),b=str(hi),value=str(v)));break
   if depth>14:
    p=bank.field(mid,mid,elo,ehi)
    if p.is_finite()and p.lower()>-ar(F(1,200)):failed.append(dict(panel=i,a=str(lo),b=str(hi),value=str(v),midpoint=str(p)));break
   stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  cells.sort(key=lambda c:F(c['a']));passed=not failed and not stack
  if passed:assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==256 and all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
  records.append(dict(panel=i,status='PASS'if passed else'FAILED',eta_lo_exact=str(elo),eta_hi_exact=str(ehi),cells=cells,trials=trials));print(json.dumps(dict(panel=i,status=records[-1]['status'],cells=len(cells),trials=trials,seconds=time.time()-tic)),flush=True)
  if failed:break
 files=[Path(__file__),Path(__file__).with_name('angleupper95_models.py'),Path(us.__file__),Path(us.ContractedRegion.numerator.__code__.co_filename),Path(us.coefficient_polynomial.__code__.co_filename)]
 if reference:files.append(Path(a.replay).resolve())
 manifest.update({str(p.relative_to(ROOT)):sha(p)for p in files})
 for p,h in manifest.items():assert sha(ROOT/p)==h,p
 out=dict(status='PASS_ANGLEupper95_UNIFORM_UV'if not failed else'FAILED_ENVELOPE',witness_sha256=SHA,objective_exact=list(map(str,OBJECTIVE)),lmin=450,tau_domain=['0','256'],eta_panels=4,panel_indices=indices,kappa_exact='1/200',precision_bits=a.precision,records=records,failed=failed,crossing_tail_coefficients=list(map(str,bank.tail_coeff)),reference_certificate_sha256=sha(a.replay)if reference else None,sha256_manifest=manifest,seconds=time.time()-tic)
 (OUT/f'uniform4_p{min(indices)}_{max(indices)}_{a.precision}_k1_200.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
