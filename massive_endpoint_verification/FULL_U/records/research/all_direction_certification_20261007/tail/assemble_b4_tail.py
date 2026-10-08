"""Complete old150 b4 omitted-spin tail, with deliberately conservative A reuse."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/all_direction_certification_20261007/tail/old150_extension'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 hashes={};witness=None;tail=[];norms=[]
 def read(p):
  d=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
  for path,h in d['sha256_manifest'].items():assert sha(ROOT/path)==h,path;hashes[path]=h
  return d
 for bits in [280,360]:
  p=OUT/f'partial_components_{bits}.json';d=read(p);assert d['status']=='SCOPED_COMPONENTS_PASS'
  witness=d['witness_sha256']if witness is None else witness;assert d['witness_sha256']==witness
  r=next(r for r in d['results']if F(r['b'])==4)
  assert F(r['delta_X4'])==F(72,5)and r['compact_spin_min']==46 and F(r['compact_V_upper_exact'])<-1
  assert r['UV_spin_min']==42 and r['UV_mu_domain']==['49','infinity']
  eta=F(r['UV_support_best']['eta_exact']);assert eta>0 and F(r['large_tau_upper_exact'])<-eta
  end=F(0);positive=F(0)
  for c in d['uv_cells']:
   a,b=F(c['tau_a']),F(c['tau_b']);assert a==end and b>a;end=b
   mass=F(c['p_mass_upper_exact']);assert mass>=0
   positive+=2*max(F(0),F(c['V_over_p_upper_by_slope']['4'])+eta)*mass
  assert end==256
  A=F(d['real_norm_A_all_spins42_up_exact']);value=A/(2*eta)+positive
  assert value==F(r['UV_support_best']['UV_support_upper_exact'])
  tail.append(value);norms.append(A)
 for bits in [320,400]:
  d=read(OUT/f'b4_bridge_{bits}.json')
  assert d['status']=='PASS_COMPACT_BRIDGE'and d['witness_sha256']==witness and d['spins']==[42,44]
  assert F(d['b'])==4 and F(d['delta_X4'])==F(72,5)and F(d['kappa_exact'])==1 and d['mu_domain']==['4','49']
  end=F(4)
  for c in d['cells']:
   a,b=F(c['a']),F(c['b']);assert a==end and b>a;end=b
   assert set(c['upper_by_spin'])=={'42','44'}and all(F(v)<-1 for v in c['upper_by_spin'].values())
  assert end==49
  if bits==400:assert d['replay_input_sha256']==sha(OUT/'b4_bridge_320.json')
 A=max(norms);compact=A/2;uv=max(tail);total=compact+uv
 hashes[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
 out=dict(status='PASS_COMPLETE_OMITTED_SPIN_TAIL',witness_sha256=witness,b='4',delta_X4='72/5',spins='Every even ell>=42',energy_domain='Every mu>4',
  compact_kappa_exact='1',compact_A_global_upper_exact=str(A),compact_support_upper_exact=str(compact),UV_support_upper_exact=str(uv),tail_support_upper_exact=str(total),tail_support_upper_float=float(total),
  note='The same global A40 norm is conservatively used once for the compact part and again inside the UV estimate. This double counting only weakens the upper bound. This file does not include spins0..40.',
  verification='Fresh high-precision replays and exact coverage/hash checks; new Jensen envelopes still require the stated analytic review.',sha256_manifest=hashes)
 (OUT/'b4_complete_tail.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
