"""Bind a tighter finite bridge without changing the global infinite-spin proof."""
from assemble_robust_uv import *
def main():
 tic=time.time();manifest={};seen=[];n=0;worst=None
 def read(name):
  p=OUT/name;d=json.loads(p.read_text());manifest[str(p.relative_to(ROOT))]=sha(p)
  for f,h in d['sha256_manifest'].items():assert sha(ROOT/f)==h;manifest[f]=h
  assert d.get('witness_sha256',WITNESS)==WITNESS
  return d
 base=read('complete_uv_coercivity.json');assert base['status']=='PASS_COMPLETE_UV_COERCIVITY'
 for a,b in RANGES:
  ref=None
  for bits in[400,480]:
   name=f'finite_uv_strong_{a}_{b}_{bits}.json';d=read(name)
   assert d['status']=='PASS_FINITE_UV'and not d['failed']and d['spins']==list(range(a,b+1,2))
   assert F(d['kappa_exact'])==F(1,100)and d['precision_bits']==bits
   if ref:assert d['reference_certificate_sha256']==sha(OUT/ref)and not d['refinement']
   assert len(d['records'])==len(d['spins'])
   for l,r in zip(d['spins'],d['records']):
    assert l==r['spin']and r['status']=='PASS'
    cells=r['cells'];end=min(F(256),F(l*(l+1),45))
    assert F(r['tau_end_exact'])==end and F(cells[0]['a'])==0 and F(cells[-1]['b'])==end
    assert all(F(c['a'])<F(c['b'])and F(c['upper_exact'])<-F(1,100)for c in cells)
    assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
    if bits==480:
     n+=len(cells);uu=max(F(c['upper_exact'])for c in cells);worst=uu if worst is None else max(worst,uu)
   ref=name
  seen.extend(d['spins'])
 assert seen==list(range(42,450,2))
 manifest[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
 out=dict(status='PASS_STRONG_FINITE_UV_COERCIVITY',witness_sha256=WITNESS,objective_exact=['0','1'],finite_spin_domain=[42,448],energy_domain=['49','infinity'],finite_kappa_exact='1/100',infinite_spin_domain='All even ell>=450',infinite_kappa_exact=base['kappa_exact'],finite_bridge_leaf_count=n,finite_bridge_worst_upper_exact=str(worst),
  proof='Stronger .01 bounds cover all finite bridge cells. The original L108 large-tau bound has kappa1, so extends the stronger finite bridge to all energy. Infinite spins retain the original .001 certificate. All original certificates remain immutable.',seconds=time.time()-tic,sha256_manifest=manifest)
 (OUT/'strong_uv_coercivity.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
