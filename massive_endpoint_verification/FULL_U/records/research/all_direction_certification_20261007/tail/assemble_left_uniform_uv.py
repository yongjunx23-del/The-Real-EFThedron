"""Bind the complete integer-spin>=450 UV theorem and exact domain seams."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/all_direction_certification_20261007/tail/left/refined'
W='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
GROUPS=[(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,14),(15,15)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 tic=time.time();manifest={};seen=set()
 def ingest(p,expected=None):
  p=Path(p).resolve();key=str(p.relative_to(ROOT));h=sha(p)
  if expected is not None:assert h==expected,key
  if key in manifest:assert manifest[key]==h
  manifest[key]=h
  if p.suffix!='.json':return None
  d=json.loads(p.read_text())
  if key not in seen:
   seen.add(key)
   for q,hh in d.get('sha256_manifest',{}).items():ingest(ROOT/q,hh)
  return d
 audit=ingest(ROOT/'output/competitive_20261007/left/uniform_source_independent_audit.json')
 assert audit['status']=='PASS_EXACT_ETA_SPECIALIZATION_AND_INTERVAL_COMPARISON'and audit['witness_sha256']==W
 allpanels=[];total=0;worst=None
 for a,b in GROUPS:
  prior=None
  for bits in[400,480]:
   path=OUT/f'uniform_uv_L450_p{a}_{b}_{bits}.json';d=ingest(path)
   assert d['status']=='PASS_UNIFORM_SOURCE_UV'and not d['failed']and d['witness_sha256']==W
   assert d['objective_exact']==['-1','0']and d['lmin']==450 and d['tau_domain']==['0','256']and d['eta_panels']==16
   assert d['panel_indices']==list(range(a,b+1))and d['precision_bits']==bits and F(d['kappa_exact'])==F(1,1000)
   assert all(F(x)<=0 for x in d['crossing_tail_Bernstein_coefficients'])
   if prior:assert d['reference_certificate_sha256']==sha(prior)
   assert len(d['records'])==b-a+1
   for i,r in zip(d['panel_indices'],d['records']):
    assert r['panel']==i and r['status']=='PASS'and F(r['eta_lo_exact'])==F(i,16*450*451)and F(r['eta_hi_exact'])==F(i+1,16*450*451)
    cells=r['cells'];assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==256
    assert all(F(c['a'])<F(c['b'])and F(c['upper_exact'])<-F(1,1000)for c in cells)
    assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
    if bits==480:
     total+=len(cells);u=max(F(c['upper_exact'])for c in cells);worst=u if worst is None else max(worst,u)
   prior=path
  allpanels.extend(d['panel_indices'])
 assert allpanels==list(range(16))and F(4)+F(450*451,256)>49
 big=ingest(OUT/'large_tau_complete.json')
 assert big['status']=='PASS_LEFT_LARGE_TAU'and big['witness_sha256']==W and big['lmin']==108 and big['tau_domain']==['256','infinity']and big['mu_domain']==['49','infinity']and F(big['kappa_exact'])>=F(1,1000)
 ingest(Path(__file__))
 out=dict(status='PASS_LEFT_INFINITE_SPIN_UV',witness_sha256=W,objective_exact=['-1','0'],lmin=450,mu_domain=['49','infinity'],kappa_exact='1/1000',tau_scheduled_margins=[dict(tau=['0','256'],kappa_exact='1/1000'),dict(tau=['256','infinity'],kappa_exact=big['kappa_exact'])],eta_panels=16,verified_cell_count=total,worst_directed_upper_exact=str(worst),precision_bits=[400,480],independent_exact_specialization_audit=True,nested_source_manifests_verified=True,seconds=time.time()-tic,sha256_manifest=manifest)
 (OUT/'infinite_uv_complete.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
