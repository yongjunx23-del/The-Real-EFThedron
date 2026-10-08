from pathlib import Path
import sys
_ADAPT_ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(_ADAPT_ROOT/'research/weak_regions_20261007/tail269p5'),str(_ADAPT_ROOT/'research/boundary_length_20261007/compact')]
"""Exact complete compact269p5 rectangle coverage; fixed-witness component only."""
import assemble_angle85 as common
from pathlib import Path
from fractions import Fraction as F
import json
ROOT=common.ROOT
OUT=ROOT/'output/weak_regions_20261007/fresh269p5/compact'
SHA='f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae'
OBJECTIVE=['-2515252609197481/288230376151711744', '-9006856288192973/9007199254740992']
common.OUT=OUT;common.SHA=SHA;common.OBJECTIVE=OBJECTIVE;common.K=F(1)
load,bind,finite,witness=common.load,common.bind,common.finite,common.witness

def main():
    low=[load(OUT/f'low/analytic_piecewise_{bits}.json')for bits in[400,512]]
    rows=[load(OUT/f'analytic_piecewise_{bits}.json')for bits in[400,512]]
    for j in low+rows:
        witness(j);assert j['status']=='PASS_PIECEWISE_ALLSPIN_COMPACT_ENVELOPES'
        for r in j['records']:
            assert F(r['V_over_p_upper_exact'])<-F(1)
            if r['method']=='low':assert F(r['spin_derivative_upper_exact'])<0
    for j in low:
        assert len(j['records'])==1
        r=j['records'][0];assert list(map(F,r['mu_domain_exact']))==[F(4),F(9,2)] and r['minimum_even_spin']==42
    assert [(r['mu_domain_exact'],r['minimum_even_spin'])for r in rows[0]['records']]==[(r['mu_domain_exact'],r['minimum_even_spin'])for r in rows[1]['records']]
    coverage=[dict(mu_domain_exact=['4','9/2'],finite_spins=[],analytic_minimum_even_spin=42,analytic_pair=[str((OUT/f'low/analytic_piecewise_{b}.json').relative_to(ROOT))for b in[400,512]])]
    end=F(9,2)
    for r in rows[1]['records']:
        lo=max(F(r['mu_domain_exact'][0]),F(9,2));hi=F(r['mu_domain_exact'][1]);L=r['minimum_even_spin'];assert lo==end and lo<hi
        tag='piece_'+str(lo).replace('/','d')+'_'+str(hi).replace('/','d')
        fin=finite(tag,lo,hi,L)if L>42 else dict(spins=[])
        coverage.append(dict(mu_domain_exact=[str(lo),str(hi)],finite_component=fin,analytic_minimum_even_spin=L,
          analytic_pair=[str((OUT/f'analytic_piecewise_{b}.json').relative_to(ROOT))for b in[400,512]],analytic_record=r['mu_domain_exact']))
        end=hi
    assert end==49
    for p in[OUT/'moments/roy_manifest.json',OUT/'moments/strips_manifest.json',OUT/'moments/independent_audit.json']:
        j=load(p);assert str(j['status']).startswith('PASS')
    bind(ROOT/'research/boundary_length_20261007/REVIEW_NEW_ENVELOPES.md')
    bind(ROOT/'research/boundary_length_20261007/REVIEW_PIECEWISE_COMPACT.md')
    bind(ROOT/'research/weak_regions_20261007/compact269p5/REVIEW_ADAPTATION.md');bind(Path(__file__));bind(Path(common.__file__))
    out=dict(status='PASS_COMPLETE_COMPACT_COERCIVITY',witness_sha256=SHA,objective_exact=OBJECTIVE,minimum_even_spin=42,mu_domain_exact=['4','49'],lower_energy_endpoint_open=True,kappa_exact='1',
      coverage_rectangles=coverage,coverage_proof='Exact adjacent closed energy intervals cover[4,49]; mu4 excluded. Finite even42..L-2 and analytic every even>=L cover all omitted spins on each rectangle without gaps. Low allspin42 proof covers4<mu<=9/2.',
      scope='Only compact coercivity V/p<-1 for all even>=42 and4<mu<=49. No real norm, UV, head integral or complete directional support.',
      replay_scope='Every finite interval freshly recomputed at640bits after512bits. Analytic envelopes at400/512. Shared mathematical implementations, independently audited source structure and moments; not proof-assistant verification.',sha256_manifest=common.manifest)
    p=OUT/'complete_compact_coercivity.json';p.write_text(json.dumps(out,indent=2));print(out['status'],len(coverage),'rectangles',len(common.manifest),'files')
if __name__=='__main__':main()
