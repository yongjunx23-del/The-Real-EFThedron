"""Read-only exact assembly of the complete compact V/p coercivity proof."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/competitive_20261007/compact'
TAIL=ROOT/'output/all_direction_certification_20261007/tail/g2_robust_stage1'
EXPECTED='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    manifest={};components=[]
    def bind(p):
        key=str(p.relative_to(ROOT));h=sha(p)
        if key in manifest:assert manifest[key]==h
        manifest[key]=h
    def load(p):
        j=json.loads(p.read_text());bind(p)
        assert j['witness_sha256']==EXPECTED
        for key,h in j.get('sha256_manifest',{}).items():
            assert sha(ROOT/key)==h,key
            bind(ROOT/key)
        components.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),status=j['status']))
        return j
    for bits in [400,512]:
        j=load(OUT/f'low_allspin42_{bits}.json')
        assert j['status']=='PASS_ALL_SPIN_LOW_ENERGY_V'
        assert j['minimum_even_spin']==42 and list(map(F,j['energy_domain_exact']))==[4,5]
        assert F(j['kappa_exact'])==1 and F(j['V_over_p_upper_exact'])<-1
        assert F(j['spin_derivative_upper_exact'])<0
    original=None
    for bits in [512,640]:
        p=OUT/f'compact_bridge_sliced_{bits}.json';j=load(p)
        assert j['status']=='PASS_FINITE_COMPACT_BRIDGE' and j['precision_bits']==bits
        assert j['spins']==list(range(42,110,2)) and list(map(F,j['mu_domain']))==[5,49]
        assert F(j['kappa_exact'])==1 and set(map(int,j['cells_by_spin']))==set(j['spins'])
        partitions={}
        for ell,rows in j['cells_by_spin'].items():
            rows=sorted(rows,key=lambda c:F(c['a']))
            intervals=[(F(c['a']),F(c['b']))for c in rows]
            assert intervals[0][0]==5 and intervals[-1][1]==49
            assert all(a<b for a,b in intervals)
            assert all(a[1]==b[0]for a,b in zip(intervals,intervals[1:]))
            assert all(F(c['V_over_p_upper_exact'])<-1 for c in rows)
            partitions[ell]=intervals
        if original is None:original=partitions
        else:
            assert partitions==original and j['every_original_spin_interval_replayed']
            assert j['replay_input_sha256']==sha(OUT/'compact_bridge_sliced_512.json')
    for bits in [320,400]:
        j=load(TAIL/f'compact_highspin_sharp_L110_{bits}.json')
        assert j['status']=='PASS_COMPACT_HIGH_SPIN' and j['lmin']==110
        assert list(map(F,j['mu_domain']))==[4,49] and F(j['kappa_exact'])==1
        assert F(j['low_mu4_8_upper_exact'])<-1 and F(j['middle_mu8_49_upper_exact'])<-1
        assert F(j['low_spin_derivative_upper_exact'])<0
    audit=load(OUT/'moments/independent_audit.json')
    assert audit['status']=='PASS_INDEPENDENT_EXACT_BETA_MOMENTS_AND_SHARED_TABLES'
    # Bind the inherited evaluated helpers as well as direct component inputs.
    for name in ['research/certified_roy_moments.py','research/head_certification_lowspin.py',
      'research/head_certification_spin0.py','research/tail_certification/certify_finite_real.py',
      'research/certify_real_tail.py','research/tail_certification/certify_uv.py',
      'research/competitive_20261007/robust/COMPACT_BRIDGE_PROOF.md']:
        bind(ROOT/name)
    bind(Path(__file__))
    out=dict(status='PASS_COMPLETE_COMPACT_COERCIVITY',witness_sha256=EXPECTED,
        objective_exact=['0','1'],minimum_even_spin=42,mu_domain_exact=['4','49'],
        lower_energy_endpoint_open=True,kappa_exact='1',components=components,
        coverage=[dict(even_spins='all >=42',mu='4<mu<=5',method='Analytic all-spin envelope'),
          dict(even_spins='42,44,...,108',mu='5<=mu<=49',method='Exact gapless interval partitions,512-bit producer and640-bit fresh replay'),
          dict(even_spins='all >=110',mu='4<mu<=49',method='Analytic all-spin envelope')],
        finite_bridge_spin_cells=sum(map(len,original.values())),
        sha256_manifest=manifest,
        scope='The frozen (0,1) functional has V_l(mu)/p_l(mu)<-1 for all even l>=42 and every4<mu<=49. This is a sign/coercivity component; it does not by itself bound the support integral or the coupling.',
        verification='Fresh paired-precision interval computations share the mathematical evaluator. Exact source-polynomial moments and remainders, analytic all-spin inequalities, and checked hashes are assumed as described in the linked proofs; not a formal proof-assistant certificate.')
    (OUT/'complete_compact_coercivity.json').write_text(json.dumps(out,indent=2))
    print(out['status'],'spin cells',out['finite_bridge_spin_cells'],'manifest',len(manifest))

if __name__=='__main__':main()
