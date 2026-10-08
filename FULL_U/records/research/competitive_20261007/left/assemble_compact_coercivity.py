"""Exact complete compact-domain assembly for the frozen left180 witness."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/competitive_20261007/left/compact'
EXPECTED='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    manifest={};components=[];global_upper=F(-10**100)
    def bind(p):
        p=p.resolve();key=str(p.relative_to(ROOT));h=sha(p)
        if key in manifest:assert manifest[key]==h
        manifest[key]=h
    def load(p):
        j=json.loads(p.read_text());bind(p);assert j['witness_sha256']==EXPECTED
        for key,h in j.get('sha256_manifest',{}).items():assert sha(ROOT/key)==h,key;bind(ROOT/key)
        components.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),status=j['status']))
        return j
    for bits in [400,512]:
        j=load(OUT/f'analytic_compact_{bits}.json')
        assert j['status']=='PASS_LEFT_ANALYTIC_COMPACT_COERCIVITY' and j['objective_exact']==['-1','0']
        assert F(j['kappa_exact'])==10 and j['precision_bits']==bits
        low=j['low_allspin42'];high=j['high_allspin100']
        assert low['minimum_even_spin']==42 and low['mu_domain_exact']==['4','5']
        assert F(low['V_over_p_upper_exact'])<-10 and F(low['spin_derivative_upper_exact'])<0
        assert high['minimum_even_spin']==100 and high['mu_domain_exact']==['4','49']
        assert F(high['low_mu4_8']['V_over_p_upper_exact'])<-10
        assert F(high['low_mu4_8']['spin_derivative_upper_exact'])<0
        assert F(high['middle_mu8_49_upper_exact'])<-10
        global_upper=max(global_upper,F(low['V_over_p_upper_exact']),F(high['low_mu4_8']['V_over_p_upper_exact']),F(high['middle_mu8_49_upper_exact']))
    original=None
    for bits in [512,640]:
        j=load(OUT/f'compact_bridge_sliced_{bits}.json')
        assert j['status']=='PASS_FINITE_COMPACT_BRIDGE' and j['precision_bits']==bits
        assert j['spins']==list(range(42,100,2)) and j['mu_domain']==['5','49']
        assert F(j['kappa_exact'])==10 and set(map(int,j['cells_by_spin']))==set(j['spins'])
        partitions={}
        for ell,rows in j['cells_by_spin'].items():
            rows=sorted(rows,key=lambda c:F(c['a']));ivs=[(F(c['a']),F(c['b']))for c in rows]
            assert ivs[0][0]==5 and ivs[-1][1]==49 and all(a<b for a,b in ivs)
            assert all(a[1]==b[0]for a,b in zip(ivs,ivs[1:]))
            assert all(F(c['V_over_p_upper_exact'])<-10 for c in rows);partitions[ell]=ivs
            global_upper=max(global_upper,max(F(c['V_over_p_upper_exact'])for c in rows))
        if original is None:original=partitions
        else:
            assert partitions==original and j['every_original_spin_interval_replayed']
            assert j['replay_input_sha256']==sha(OUT/'compact_bridge_sliced_512.json')
    audit=load(OUT/'moments/independent_audit.json')
    assert audit['status']=='PASS_INDEPENDENT_EXACT_BETA_MOMENTS_AND_SHARED_TABLES'
    for name in ['research/certified_roy_moments.py','research/head_certification_lowspin.py',
      'research/head_certification_spin0.py','research/tail_certification/certify_finite_real.py',
      'research/certify_real_tail.py','research/tail_certification/certify_uv.py',
      'research/competitive_20261007/robust/COMPACT_BRIDGE_PROOF.md',
      'research/competitive_20261007/robust/compact_finite220_combined.py',
      'research/competitive_20261007/left/COMPACT_PROOF.md']:
        bind(ROOT/name)
    bind(Path(__file__))
    out=dict(status='PASS_COMPLETE_COMPACT_COERCIVITY',witness_sha256=EXPECTED,objective_exact=['-1','0'],
      minimum_even_spin=42,mu_domain_exact=['4','49'],lower_energy_endpoint_open=True,kappa_exact='10',
      recorded_V_over_p_upper_exact=str(global_upper),
      safe_integer_kappa_from_recorded_enclosures=(-global_upper).numerator//(-global_upper).denominator,
      components=components,finite_bridge_spin_cells=sum(map(len,original.values())),sha256_manifest=manifest,
      coverage=[dict(even_spins='all >=42',mu='4<mu<=5',method='Analytic envelope,400/512 bits'),
        dict(even_spins='42,44,...,98',mu='5<=mu<=49',method='Exact partitions,512-bit producer and640-bit fresh replay'),
        dict(even_spins='all >=100',mu='4<mu<=49',method='Analytic envelope,400/512 bits')],
      scope='The exact frozen(-1,0) functional obeys V_l(mu)/p_l(mu)<-10 for every even l>=42 and4<mu<=49. Sign/coercivity only; head integration, real norms and UV remain separate components.',
      verification='Directed interval evaluations and analytic continuum inequalities with fresh higher-precision replays. Evaluators share the mathematical implementation; independent floating checks and exact moment audits supplement them. Not formal proof-assistant verification.')
    (OUT/'complete_compact_coercivity.json').write_text(json.dumps(out,indent=2))
    print(out['status'],'spin-cells',out['finite_bridge_spin_cells'],'manifest',len(manifest))
if __name__=='__main__':main()
