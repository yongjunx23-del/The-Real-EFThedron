"""Exact complete omitted-spin tail assembly for the old150 b6 direction."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/all_direction_certification_20261007/tail/b6_extension'
W='bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 tic=time.time();manifest={};seen=set()
 def read(p):
  p=Path(p).resolve();key=str(p.relative_to(ROOT));h=sha(p)
  if key in manifest:assert manifest[key]==h,key
  manifest[key]=h
  if p.suffix!='.json':return None
  d=json.loads(p.read_text())
  if key not in seen:
   seen.add(key)
   for q,hq in d.get('sha256_manifest',{}).items():
    qp=ROOT/q;assert sha(qp)==hq,q;read(qp)
  return d
 def parameters(d):
  assert d['base_witness_sha256']==W and d['objective_exact']==['1','6']and F(d['delta_X4_exact'])==F(107,5)and F(d['kappa_exact'])==1
 fixed=read(ROOT/'output/certification_tail/fixed_tail_certificate.json')
 assert fixed['witness_sha256']==W and fixed['replayed_components_this_run']
 assert fixed['status']=='All-spin/all-energy omitted-spin tail certified under the stated amplitude assumptions.'
 A=F(fixed['A40_upper_exact']);assert 0<A<F(1,50000)
 # All spins>=42 on the two analytic compact intervals.
 for bits in[320,400]:
  d=read(OUT/f'compact_scoped_{bits}.json');parameters(d)
  assert d['status']=='PASS_B6_COMPACT_SCOPED'and d['minimum_even_spin']==42 and d['mu_domains']==[['4','5'],['8','49']]
  assert d['first_endpoint_open']and d['precision_bits']==bits and F(d['low_upper_exact'])<-1 and F(d['middle_upper_exact'])<-1
 # Finite missing compact spins42..106 on5..8.
 previous=None;bridge_cells=0
 for bits in[400,480]:
  path=OUT/f'compact_bridge_{bits}.json';d=read(path);parameters(d)
  assert d['status']=='PASS_B6_FINITE_COMPACT_BRIDGE'and d['spins']==list(range(42,108,2))and d['mu_domain']==['5','8']and d['precision_bits']==bits
  cursor=F(5)
  for c in d['cells']:
   assert cursor==F(c['a'])<F(c['b']);cursor=F(c['b'])
   assert sorted(map(int,c['upper_by_spin']))==d['spins']and all(F(v)<-1 for v in c['upper_by_spin'].values())
  assert cursor==8
  if previous:assert d['reference_certificate_sha256']==sha(previous)
  previous=path;bridge_cells=len(d['cells'])
 # Frozen Jensen packet at larger delta108/5 also bounds107/5, since X4/p>=0.
 for bits in[280,360]:
  d=read(ROOT/f'output/all_direction_certification_20261007/tail/old150_extension/partial_components_{bits}.json')
  assert d['status']=='SCOPED_COMPONENTS_PASS'and d['witness_sha256']==W and d['precision_bits']==bits
  r=next(r for r in d['results']if F(r['b'])==6)
  assert F(r['delta_X4'])==F(108,5)and r['compact_spin_min']==108 and r['compact_mu_domain']==['4','49']and F(r['compact_V_upper_exact'])<-1
  assert r['UV_spin_min']==42 and r['UV_mu_domain']==['49','infinity']and F(r['large_tau_upper_exact'])<-1
  cursor=F(8)
  for c in d['uv_cells']:
   if F(c['tau_b'])<=8:continue
   assert cursor==F(c['tau_a'])<F(c['tau_b']);cursor=F(c['tau_b'])
   assert F(c['V_over_p_upper_by_slope']['6'])<-1
  assert cursor==256
 # New uniform finite-eta source contraction covers every omitted spin at tau0..8.
 previous=None;uniform_cells=0
 for bits in[320,400]:
  path=OUT/f'uniform_uv_{bits}.json';d=read(path);parameters(d)
  assert d['status']=='PASS_B6_UNIFORM_UV'and d['lmin']==42 and d['tau_domain']==['0','8']and d['eta_panels']==8 and d['precision_bits']==bits
  assert F(d['mu_lower_bound_exact'])==F(919,4)>49 and d['crossing_omitted_kmin']==41
  assert all(F(v)<=0 for v in d['crossing_tail_coefficients'])and len(d['records'])==8
  for i,r in enumerate(d['records']):
   assert r['panel']==i and F(r['eta_lo_exact'])==F(i,8*42*43)and F(r['eta_hi_exact'])==F(i+1,8*42*43)
   cursor=F(0)
   for c in r['cells']:
    assert cursor==F(c['a'])<F(c['b'])and F(c['upper_exact'])<-1;cursor=F(c['b'])
   assert cursor==8
  if previous:assert d['reference_certificate_sha256']==sha(previous)
  previous=path;uniform_cells=sum(len(r['cells'])for r in d['records'])
 audit=read(ROOT/'output/competitive_20261007/upper_arc/b6/uniform_source_independent_audit.json')
 assert audit['status']=='PASS_EXACT_FINITE_SPIN_SPECIALIZATION_B6'and audit['base_witness_sha256']==W and audit['objective_exact']==['1','6']and F(audit['delta_X4_exact'])==F(107,5)
 assert audit['source_and_crossing_polynomial_equalities_spins']==[42,44,80,450]and all(c['P_exact_contained']for c in audit['positive_Legendre_denominator_checks'])
 read(ROOT/'research/all_direction_certification_20261007/REVIEW_B6_UNIFORM_SOURCE.md')
 read(ROOT/'research/competitive_20261007/robust/b4_tail_review.txt');read(Path(__file__))
 result=dict(status='PASS_COMPLETE_OMITTED_SPIN_TAIL',witness_sha256=W,base_witness_sha256=W,objective_exact=['1','6'],delta_X4_exact='107/5',
  minimum_even_spin=42,spins='Every even ell>=42',energy_domain='Every mu>4',kappa_exact='1',
  real_norm_upper_exact=str(A),tail_support_upper_exact=str(A/2),tail_support_upper_float=float(A/2),
  compact_bridge_cell_count=bridge_cells,uniform_UV_cell_count=uniform_cells,
  inherited_increment_lemma='For even ell>=4 and Q>0, X4/p>=ell(2ell-5)/Q^2>0, by paired Legendre roots. Hence V(delta107/5)=V(delta108/5)-(1/5)X4<=V(delta108/5).',
  coverage=dict(compact_allspins=[['4','5'],['8','49']],compact_bridge_spins=[42,106],compact_bridge_energy=['5','8'],compact_infinite_spins_min=108,UV_uniform_tau=['0','8'],UV_Jensen_tau=['8','256'],UV_unbounded_tau=['256','infinity']),
  formula='Global V<=-p implies h(U,V)<=U^2/(2p). U is unchanged from the fixed150 witness, so its global allspin>=42 real norm A40 is used exactly once.',
  scope='Complete omitted-spin tail only; finite0..40 head remains separate. Internal directed-interval computer-assisted certificate, not formal verification.',seconds=time.time()-tic,sha256_manifest=manifest)
 (OUT/'complete_tail.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items()if k not in['sha256_manifest','coverage']},indent=2))

if __name__=='__main__':main()
