"""Complete exact b=4 bound, combining independently scoped components."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,math
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from flint import arb,ctx,fmpq
from angle_benchmarks import support_report
OUT=ROOT/'output/competitive_20261007/robust/b4'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upward(x,n=12):
    k=-((-x.numerator*10**n)//x.denominator)
    return f'{k//10**n}.{k%10**n:0{n}d}'
def au(q):
    assert q.is_finite();m,e=q.upper().man_exp();return F(int(m))*F(2)**int(e)
def main():
    ctx.prec=320;hashes={}
    def read(path):
        q=json.loads(path.read_text());hashes[str(path.relative_to(ROOT))]=sha(path)
        for f,h in q['sha256_manifest'].items():assert sha(ROOT/f)==h,f;hashes[f]=h
        return q
    head=read(OUT/'head_bundle.json');finite=read(OUT/'finite_block_replayed.json')
    tail=read(ROOT/'output/all_direction_certification_20261007/tail/old150_extension/b4_complete_tail.json')
    assert head['status']=='PASS_FINITE_HEAD_0_40' and tail['status']=='PASS_COMPLETE_OMITTED_SPIN_TAIL'
    assert head['witness_sha256']==finite['witness_sha256']==tail['witness_sha256']
    assert head['direction_exact']==['1','4'] and head['crossing_increment_exact']==tail['delta_X4']=='72/5' and tail['b']=='4'
    assert finite['independent_replay_uniform_V_slack_exact']=='1/10000'
    r=finite['results'][0];assert r['status']=='PASS' and r['b']=='4' and r['delta_X4']=='72/5'
    low=sum((F(head['component_upper_exact'][str(l)])for l in [0,2,4,6]),F(0))
    finite_value=F(r['finite_support_upper_exact']);tail_value=F(tail['tail_support_upper_exact'])
    assert tail_value==F(tail['compact_support_upper_exact'])+F(tail['UV_support_upper_exact'])
    total=low+finite_value+tail_value;rounded=F(upward(total,6));assert rounded>=total
    unit=au(arb(fmpq(rounded.numerator,rounded.denominator))/arb(17).sqrt())
    theta=math.degrees(math.atan(4));dual=support_report(theta);primal=support_report(theta,'primal')
    files=[Path(__file__),ROOT/'research/competitive_20261007/robust/b4_tail_review.txt',ROOT/'research/angle_benchmarks.py']
    hashes.update({str(p.relative_to(ROOT)):sha(p)for p in files})
    out=dict(status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND',witness_sha256=head['witness_sha256'],direction_exact=['1','4'],crossing_increment_exact='72/5',
      inequality='g0+4*g2 <= '+str(rounded),support_upper_exact=str(total),support_upper_outward_6dp=upward(total,6),
      unit_normal_support_upper_exact=str(unit),unit_normal_support_upper_outward_12dp=upward(unit),angle_degrees_approx=theta,
      components_exact={'spin0_2_4_6':str(low),'spins8_40':str(finite_value),'compact_spins42_up':tail['compact_support_upper_exact'],'UV_spins42_up':tail['UV_support_upper_exact']},
      coverage='Every even spin ell>=0 and every mu>4. Finite0..40 allenergy; compact42,44 by intervalbridge; compact>=46 by allspin envelope; UV>=42 by positive-V mass estimate and large-tau closure.',
      exact_witness_definition='Original150 multipliers interpreted as exact binary64 rationals, plus exact72/5 X4 increment; objective exactly(1,4). The convenience family_b4.npz is not this exact witness definition.',
      graphical_comparison={'dual':dual,'primal':primal,'certified_minus_dual_estimate':float(unit)-dual['support'],'certified_minus_primal_estimate':float(unit)-primal['support'],'relative_to_dual_estimate_percent':100*(float(unit)/dual['support']-1)},
      comparison_warning='Published dual support is missing-cap-sensitive here; isolated maxg0 anchor dominates. Graphical uncertainty is not a rigorous published numerical interval. No claim of dominating the published bound at this angle.',
      verification_scope='Arb interval proofs, exact rational partitions and moment identities, higher-precision replays, exact component assembly and structural analytic review under the same dispersion/analyticity/unitarity assumptions. This is not formal theorem-prover verification, nor independent implementations of every kernel.',sha256_manifest=hashes)
    (OUT/'complete_b4_certificate.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['sha256_manifest','graphical_comparison']},indent=2))
if __name__=='__main__':main()
