"""Directed interval proof with an exact rectangular eta/tau partition."""
from uniform_source_uv import *
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lmin',type=int,default=450);ap.add_argument('--precision',type=int,default=400);ap.add_argument('--eta-panels',type=int,default=16);ap.add_argument('--panels',default='0:15');ap.add_argument('--kappa',default='1/1000');ap.add_argument('--replay');a=ap.parse_args();ctx.prec=a.precision;tic=time.time();bank=UniformSource(a.lmin);kap=F(a.kappa)
 reference=json.loads(Path(a.replay).read_text())if a.replay else None
 if reference:
  assert reference['status']=='PASS_UNIFORM_SOURCE_UV'and reference['lmin']==a.lmin and reference['eta_panels']==a.eta_panels and F(reference['kappa_exact'])==kap
  for p,h in reference['sha256_manifest'].items():assert sha(ROOT/p)==h,p
  indices=reference['panel_indices']
 else:indices=list(range(int(a.panels.split(':')[0]),int(a.panels.split(':')[1])+1))if':'in a.panels else list(map(int,a.panels.split(',')))
 assert all(0<=i<a.eta_panels for i in indices);records=[];failed=[]
 for i in indices:
  t=time.time();elo=F(i,a.eta_panels*bank.T);ehi=F(i+1,a.eta_panels*bank.T)
  cuts=list(map(F,['0','1','2','21/10','11/5','23/10','12/5','5/2','3','5','10','16','21','24','32','64','128','256']))
  stack=[(x,y,0)for x,y in zip(cuts,cuts[1:])]
  if reference:stack=[(F(c['a']),F(c['b']),0)for r in reference['records']if r['panel']==i for c in r['cells']]
  cells=[];count=0
  while stack:
   lo,hi,depth=stack.pop();count+=1;v=bank.field(lo,hi,elo,ehi)
   if v.is_finite()and v.upper()<-ar(kap):cells.append(dict(a=str(lo),b=str(hi),upper_exact=str(upper(v))));continue
   if reference or depth>=22:failed.append(dict(panel=i,a=str(lo),b=str(hi),value=str(v)));break
   mid=(lo+hi)/2
   if depth>=12:
    point=bank.field(mid,mid,elo,ehi)
    if point.is_finite()and point.lower()>-ar(kap):failed.append(dict(panel=i,a=str(lo),b=str(hi),point=str(point)));break
   stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  cells.sort(key=lambda c:F(c['a']));passed=not failed and not stack
  if passed:
   assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==256
   assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
  records.append(dict(panel=i,eta_lo_exact=str(elo),eta_hi_exact=str(ehi),status='PASS'if passed else'FAILED',cells=cells,evaluations=count,seconds=time.time()-t))
  print(json.dumps(dict(panel=i,status=records[-1]['status'],cells=len(cells),evaluations=count,seconds=time.time()-t)),flush=True)
  if failed:break
 prep=json.loads((FOLDER/'prepared_components.json').read_text());manifest=dict(prep['sha256_manifest'])
 files=[Path(__file__),Path(UniformSource.__init__.__code__.co_filename),Path(ContractedRegion.numerator.__code__.co_filename),Path(Region.field.__code__.co_filename),Path(coefficient_polynomial.__code__.co_filename),FOLDER/'prepared_components.json']
 if reference:files.append(Path(a.replay).resolve())
 for p in files:manifest[str(p.relative_to(ROOT))]=sha(p)
 for p,h in manifest.items():assert sha(ROOT/p)==h,p
 out=dict(status='PASS_UNIFORM_SOURCE_UV'if not failed else'FAILED_ENVELOPE',witness_sha256=sha(LEFT),objective_exact=['-1','0'],lmin=a.lmin,tau_domain=['0','256'],eta_panels=a.eta_panels,panel_indices=indices,kappa_exact=str(kap),precision_bits=a.precision,records=records,failed=failed,reference_certificate_sha256=sha(a.replay)if reference else None,crossing_tail_Bernstein_coefficients=list(map(str,bank.tail_coeff)),acceptance_policy='Every cell finite, directed upper strictly below negative kappa; exact rational domain partition.',seconds=time.time()-tic,sha256_manifest=manifest)
 name=f'uniform_uv_L{a.lmin}_p{min(indices)}_{max(indices)}_{a.precision}.json';(FOLDER/name).write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
