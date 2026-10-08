"""Exact energy/spin rectangle coverage of the angle85 compact component.

Only saved, paired interval proofs are assembled. This is a coercivity
component, not a complete support integral or an independent implementation.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/boundary_length_20261007/compact/angle85'
SHA='58fafda2eafa3a741610175f416a009c84f9a18f309d36e74d0e4007f65fd6d8'
OBJECTIVE=['24532160660097/281474976710656','4486462071114449/4503599627370496']
K=F(1,4)
manifest={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bind(p):
    p=p.resolve();rel=str(p.relative_to(ROOT));h=sha(p)
    if rel in manifest:assert manifest[rel]==h,rel
    manifest[rel]=h
def load(p):
    j=json.loads(p.read_text());bind(p)
    for q,h in j.get('sha256_manifest',{}).items():
        pp=ROOT/q;assert sha(pp)==h,q;bind(pp)
    return j
def witness(j):
    assert j['witness_sha256']==SHA
    if 'objective_exact'in j:assert list(map(F,j['objective_exact']))==list(map(F,OBJECTIVE))
    assert F(j['kappa_exact'])==K
def finite(tag,lo,hi,L):
    p=OUT/f'{tag}_512.json';q=OUT/f'{tag}_640.json';a,b=load(p),load(q)
    assert b['replay_input_sha256']==sha(p) and b['every_original_spin_interval_replayed'] is True
    assert a['precision_bits']==512 and b['precision_bits']==640
    n=0;maxup=None
    for j in[a,b]:
        witness(j);assert j['status']=='PASS_FINITE_COMPACT_BRIDGE'
        assert list(map(F,j['mu_domain']))==[lo,hi]
        assert j['spins']==list(range(42,L,2))
        assert set(map(int,j['cells_by_spin']))==set(j['spins'])
        for ell,rows in j['cells_by_spin'].items():
            cells=sorted(rows,key=lambda r:F(r['a']))
            assert F(cells[0]['a'])==lo and F(cells[-1]['b'])==hi
            end=lo
            for c in cells:
                x,y,v=F(c['a']),F(c['b']),F(c['V_over_p_upper_exact'])
                assert x==end and x<y and v<-K
                end=y;maxup=v if maxup is None else max(maxup,v)
            assert end==hi
        n+=sum(map(len,j['cells_by_spin'].values()))
    for ell in a['cells_by_spin']:
        geom=lambda j:sorted((F(c['a']),F(c['b']))for c in j['cells_by_spin'][ell])
        assert geom(a)==geom(b)
    return dict(spins=list(range(42,L,2)),producer=str(p.relative_to(ROOT)),replay=str(q.relative_to(ROOT)),paired_spin_cells=n,V_over_p_upper_exact=str(maxup))

def main():
    old=[load(OUT/f'analytic_L112_{bits}.json')for bits in[400,512]]
    new=[load(OUT/f'analytic_piecewise_{bits}.json')for bits in[400,512]]
    for j in old+new:witness(j)
    for j in old:
        assert j['status']=='PASS_ANGLE85_ANALYTIC_COMPACT_COERCIVITY'
        low=j['low_allspin42'];assert low['minimum_even_spin']==42 and list(map(F,low['mu_domain_exact']))==[F(4),F(19,4)]
        assert F(low['V_over_p_upper_exact'])<-K and F(low['spin_derivative_upper_exact'])<0
        high=j['high_allspin'];assert high['minimum_even_spin']==112 and list(map(F,high['mu_domain_exact']))==[F(4),F(49)]
        assert F(high['low_mu4_8']['V_over_p_upper_exact'])<-K and F(high['low_mu4_8']['spin_derivative_upper_exact'])<0 and F(high['middle_mu8_49_upper_exact'])<-K
    for j in new:
        assert j['status']=='PASS_PIECEWISE_ALLSPIN_COMPACT_ENVELOPES'
        for r in j['records']:
            assert F(r['V_over_p_upper_exact'])<-K
            if r['method']=='low':assert F(r['spin_derivative_upper_exact'])<0
    assert [(r['mu_domain_exact'],r['minimum_even_spin'])for r in new[0]['records']]==[(r['mu_domain_exact'],r['minimum_even_spin'])for r in new[1]['records']]
    plan=[r for r in new[1]['records']if F(r['mu_domain_exact'][1])<=8 or F(r['mu_domain_exact'][0])>=24]
    plan += [dict(mu_domain_exact=['8','12'],minimum_even_spin=112,tag='slice_8_12'),dict(mu_domain_exact=['12','24'],minimum_even_spin=112,tag='slice_12_24')]
    plan.sort(key=lambda r:F(r['mu_domain_exact'][0]))
    coverage=[dict(mu_domain_exact=['4','19/4'],finite_spins=[],analytic_minimum_even_spin=42,analytic_pair=[str((OUT/f'analytic_L112_{b}.json').relative_to(ROOT))for b in[400,512]],analytic_record='low_allspin42')]
    end=F(19,4)
    for r in plan:
        aa,bb=r['mu_domain_exact'];lo,hi=F(aa),F(bb);L=r['minimum_even_spin'];assert lo==end and lo<hi
        tag=r.get('tag','piece_'+aa.replace('/','d')+'_'+bb.replace('/','d'))
        fin=finite(tag,lo,hi,L)
        special='tag'in r
        pair=[str((OUT/f'analytic_L112_{b}.json'if special else OUT/f'analytic_piecewise_{b}.json').relative_to(ROOT))for b in[400,512]]
        coverage.append(dict(mu_domain_exact=[str(lo),str(hi)],finite_component=fin,analytic_minimum_even_spin=L,analytic_pair=pair,analytic_record='high_allspin'if special else r['mu_domain_exact']))
        end=hi
    assert end==49
    for p in[OUT/'moments/roy_manifest.json',OUT/'moments/strips_manifest.json',OUT/'moments/independent_audit.json']:
        j=load(p);assert str(j['status']).startswith('PASS')
    bind(ROOT/'research/boundary_length_20261007/REVIEW_NEW_ENVELOPES.md')
    bind(ROOT/'research/boundary_length_20261007/REVIEW_PIECEWISE_COMPACT.md');bind(Path(__file__))
    j=dict(status='PASS_COMPLETE_COMPACT_COERCIVITY',witness_sha256=SHA,objective_exact=OBJECTIVE,minimum_even_spin=42,mu_domain_exact=['4','49'],lower_energy_endpoint_open=True,kappa_exact='1/4',
      coverage_rectangles=coverage,coverage_proof='The exact energy intervals form [4,49] without gaps. In each interval, finite even spins42..L-2 have paired512/640 interval certificates and every even spin>=L has paired400/512 analytic certificates. The endpoint mu4 is excluded. Shared boundaries are covered by both adjacent closed proofs.',
      scope='Only V/p<-1/4 on all even spins>=42,4<mu<=49. No real norm, UV domain, finite-head integral or complete directional support is included.',
      replay_scope='Fresh higher-precision evaluation of every finite interval and paired evaluation of the analytic envelopes, sharing the mathematical implementations. Independent moment checks and source-structure review are bound; not proof-assistant verification.',sha256_manifest=manifest)
    p=OUT/'complete_compact_coercivity.json';p.write_text(json.dumps(j,indent=2));print(j['status'],len(coverage),'rectangles',len(manifest),'files')
if __name__=='__main__':main()
