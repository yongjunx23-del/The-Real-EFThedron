"""Exact final assembly of the b6 continuum support bound."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,sys,math
ROOT=Path(__file__).resolve().parents[4];sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from flint import arb,ctx,fmpq
from angle_benchmarks import support_report
OUT=ROOT/'output/competitive_20261007/upper_arc/b6'
EXPECTED='bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upward(q,n=12):
    k=-((-q.numerator*10**n)//q.denominator);return f'{k//10**n}.{k%10**n:0{n}d}'
def upper(q):
    assert q.is_finite();m,e=q.upper().man_exp();return F(int(m))*F(2)**int(e)
def main():
    ctx.prec=320;hashes={}
    def read(p):
        q=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
        for f,h in q['sha256_manifest'].items():assert sha(ROOT/f)==h,f;hashes[f]=h
        return q
    head=read(OUT/'head_bundle.json')
    tail=read(ROOT/'output/all_direction_certification_20261007/tail/b6_extension/complete_tail.json')
    assert head['status']=='PASS_FINITE_HEAD_0_40'and tail['status']=='PASS_COMPLETE_OMITTED_SPIN_TAIL'
    assert head['witness_sha256']==tail['witness_sha256']==EXPECTED
    assert head['objective_exact']==tail['objective_exact']==['1','6']
    assert head['crossing_increment_exact']==tail['delta_X4_exact']=='107/5'
    assert tail['minimum_even_spin']==42 and tail['energy_domain']=='Every mu>4'
    assert tail['coverage']==dict(compact_allspins=[['4','5'],['8','49']],compact_bridge_spins=[42,106],
      compact_bridge_energy=['5','8'],compact_infinite_spins_min=108,UV_uniform_tau=['0','8'],UV_Jensen_tau=['8','256'],UV_unbounded_tau=['256','infinity'])
    assert F(tail['kappa_exact'])==1 and 2*F(tail['tail_support_upper_exact'])==F(tail['real_norm_upper_exact'])
    audit=read(OUT/'uniform_source_independent_audit.json')
    assert audit['status']=='PASS_EXACT_FINITE_SPIN_SPECIALIZATION_B6'and audit['base_witness_sha256']==EXPECTED
    assert audit['objective_exact']==['1','6']and audit['delta_X4_exact']=='107/5'
    assert all(c['P_exact_contained']is True for c in audit['positive_Legendre_denominator_checks'])
    h=F(head['finite_head_upper_exact']);t=F(tail['tail_support_upper_exact']);assert h>0 and t>=0
    assert h==sum(map(F,head['component_upper_exact'].values()),F(0))
    total=h+t;rounded=F(upward(total,6));assert rounded>=total
    unit=upper(arb(fmpq(rounded.numerator,rounded.denominator))/arb(37).sqrt())
    angle=math.degrees(math.atan(6));dual=support_report(angle);primal=support_report(angle,'primal')
    for p in [Path(__file__),Path(__file__).with_name('HEAD_PROOF.txt'),ROOT/'research/angle_benchmarks.py',
      ROOT/'research/all_direction_certification_20261007/REVIEW_B6_UNIFORM_SOURCE.md',OUT/'uniform_source_independent_audit.json']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    out=dict(status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND',witness_sha256=EXPECTED,objective_exact=['1','6'],direction_exact=['1','6'],
      crossing_increment_exact='107/5',inequality='g0+6*g2 <= '+upward(total,6),support_upper_exact=str(total),support_upper_outward_6dp=upward(total,6),
      unit_normal_support_upper_exact=str(unit),unit_normal_support_upper_outward_12dp=upward(unit),angle_degrees_approx=angle,
      components_exact={'finite_head0_40':str(h),'omitted_spins42_up':str(t)},
      coverage='Every even spin ell>=0 and every mu>4. Finite0..40 covers allenergy; the independently assembled tail packet covers all evenell>=42 and allmu>4.',
      exact_witness_definition='Original150 multipliers interpreted as exact binary64 rationals, plus exact107/5 X4 increment; objective exactly(1,6). No rounded convenience NPZ defines this witness.',
      graphical_comparison={'dual':dual,'primal':primal,'certified_minus_dual_estimate':float(unit)-dual['support'],'certified_minus_primal_estimate':float(unit)-primal['support'],'relative_to_dual_estimate_percent':100*(float(unit)/dual['support']-1)},
      comparison_warning='The published contour is graphical and missing-cap-sensitive at this angle; its isolated maxg0 anchor may dominate the support estimate. Extraction uncertainty is not a rigorous interval. This comparison does not prove dominance of the published dual.',
      verification_scope='Arb interval bounds with analytic remainders, exact rational coefficient conventions, exact partitions and component assembly, and fresh higher-precision replays. The precision runs share kernel implementations. No formal theorem-prover or independent implementation claim is made.',
      sha256_manifest=hashes)
    (OUT/'complete_b6_certificate.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items()if k not in['sha256_manifest','graphical_comparison']},indent=2),flush=True)
    print('GRAPHICAL_RELATIVE_PERCENT',out['graphical_comparison']['relative_to_dual_estimate_percent'],flush=True)
if __name__=='__main__':main()
