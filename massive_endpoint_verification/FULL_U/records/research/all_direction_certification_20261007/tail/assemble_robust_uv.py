"""Exact coverage/hash assembly of the robust witness's complete UV tail."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/all_direction_certification_20261007/tail/g2_robust_stage1'
WITNESS='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'
RANGES=[(42,92),(94,144),(146,196),(198,248),(250,300),(302,352),(354,400),(402,448)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 tic=time.time();manifest={};components=[]
 def read(name):
  p=OUT/name;d=json.loads(p.read_text());manifest[str(p.relative_to(ROOT))]=sha(p)
  for f,h in d['sha256_manifest'].items():
   assert sha(ROOT/f)==h,f
   if f in manifest:assert manifest[f]==h
   manifest[f]=h
  if'witness_sha256'in d:assert d['witness_sha256']==WITNESS
  return d
 rebuild=read('finite_source_exact_rebuild.json');assert rebuild['status']=='PASS_EXACT_RATIONAL_REBUILD'and rebuild['all_entries_equal']
 allspins=[];leaves=0;worst=None
 for a,b in RANGES:
  prior=None
  for bits in[400,480]:
   name=f'finite_uv_bridge_{a}_{b}_{bits}.json';d=read(name)
   assert d['status']=='PASS_FINITE_UV'and not d['failed']and d['precision_bits']==bits
   assert d['spins']==list(range(a,b+1,2))and F(d['kappa_exact'])==F(1,1000)
   assert len(d['records'])==len(d['spins'])and all(F(c)<=0 for c in d['crossing_tail_Bernstein_coefficients'])
   if prior:assert d['reference_certificate_sha256']==sha(OUT/prior)
   for l,rec in zip(d['spins'],d['records']):
    assert l==rec['spin']and rec['status']=='PASS'
    end=min(F(256),F(l*(l+1),45));assert F(rec['tau_end_exact'])==end
    cells=rec['cells'];assert cells and F(cells[0]['a'])==0 and F(cells[-1]['b'])==end
    assert all(F(c['a'])<F(c['b'])and F(c['upper_exact'])<-F(1,1000)for c in cells)
    assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
    if bits==480:
     leaves+=len(cells);uu=max(F(c['upper_exact'])for c in cells);worst=uu if worst is None else max(worst,uu)
   prior=name
  allspins.extend(d['spins']);components.append(dict(kind='finite_spin_bridge',spins=[a,b],certificates=[f'finite_uv_bridge_{a}_{b}_{bits}.json'for bits in[400,480]]))
 assert allspins==list(range(42,450,2))
 for bits in[320,400]:
  d=read(f'contracted_finite_uv_L450_proof_{bits}.json')
  assert d['status']=='PASS'and d['lmin']==450 and d['tau_domain']==['0','256']and not d['failed']
  assert F(d['kappa_exact'])==F(1,1000)and d['objective_exact']==['0','1']
  cells=d['leaves'];assert F(cells[0]['lo'])==0 and F(cells[-1]['hi'])==256
  assert all(F(x['lo'])<F(x['hi'])and F(x['upper_exact'])<-F(1,1000)for x in cells)
  assert all(F(x['hi'])==F(y['lo'])for x,y in zip(cells,cells[1:]))
  d=read(f'large_tau_L108_{bits}.json')
  assert d['status']=='PASS_UNBOUNDED_TAU'and d['lmin']==108 and d['mu_domain']==['49','infinity']and d['tau_domain']==['256','infinity']
  assert F(d['V_over_p_upper_exact'])<-1 and F(d['kappa_exact'])==1
 # Every smaller even spin has tau<=106*107/45<256 at mu>=49.
 assert F(106*107,45)<256
 components.extend([dict(kind='infinite_spin_small_tau',lmin=450,tau=['0','256']),dict(kind='all_spin_large_tau',lmin=108,tau=['256','infinity'])])
 manifest[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
 out=dict(status='PASS_COMPLETE_UV_COERCIVITY',witness_sha256=WITNESS,objective_exact=['0','1'],spin_domain='All even ell>=42',energy_domain=['49','infinity'],kappa_exact='1/1000',finite_bridge_spin_count=len(allspins),finite_bridge_leaf_count=leaves,finite_bridge_worst_upper_exact=str(worst),source_exact_rebuild=True,paired_precision_or_higher_precision_replay=True,components=components,
  proof='The finite-spin bridge covers ell42..448 and tau<=min(256,T/45). The all-spin L450 component covers tau<=256 for ell>=450. The L108 unbounded-tau component covers tau>=256; smaller physical spins have tau<256. All intervals and seams are checked with exact rational comparisons.',seconds=time.time()-tic,sha256_manifest=manifest)
 (OUT/'complete_uv_coercivity.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in['sha256_manifest','components']},indent=2))
if __name__=='__main__':main()
