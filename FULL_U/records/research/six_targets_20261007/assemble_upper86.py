"""Fail-closed continuum assembly for the repaired83 peak witness."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research/boundary_length_20261007'),str(ROOT/'research')]
from assemble_complete_v3 import Inputs as BaseInputs,sha,head_cost,real_norms,ceiling_decimal
from angle_benchmarks import support_report
BASE=ROOT/'output/six_targets_20261007/upper86';HEAD=BASE/'head'
WITNESS=ROOT/'output/weak_regions_20261007/upper/angle86_refined.npz'
SHA='14e2b30df27e95ede876a8ad75ca190904547fb4ace535c14c6a11725d2df97f'

class Inputs(BaseInputs):
    def read(self,path,expected=None):
        p=self.bind(path,expected);d=json.loads(p.read_text())
        if p not in self.seen:
            self.seen.add(p)
            for sub,h in (d.get('sha256_manifest',{}).items()if isinstance(d,dict)else []):
                q=self.bind(sub,h)
                if q.suffix=='.json':self.read(q,h)
        return d

def compact(inputs,direction):
    d=inputs.read(BASE/'compact/complete_compact_coercivity.json')
    assert d['status']=='PASS_COMPLETE_COMPACT_COERCIVITY'and d['witness_sha256']==SHA and d['objective_exact']==direction
    assert d['minimum_even_spin']==42 and d['mu_domain_exact']==['4','49']and d['lower_energy_endpoint_open']
    k=F(d['kappa_exact']);assert k>0;end=F(4);nc=na=0
    for r in d['coverage_rectangles']:
        lo,hi=map(F,r['mu_domain_exact']);L=r['analytic_minimum_even_spin'];assert end==lo<hi<=49;end=hi
        assert len(r['analytic_pair'])==2
        for p,bits in zip(r['analytic_pair'],[400,512]):
            j=inputs.read(p);assert j['status']=='PASS_PIECEWISE_ALLSPIN_COMPACT_ENVELOPES'and j['precision_bits']==bits
            assert j['witness_sha256']==SHA and j['objective_exact']==direction
            matches=[a for a in j['records']if a['mu_domain_exact']==r['analytic_record']];assert len(matches)==1;a=matches[0]
            assert list(map(F,a['mu_domain_exact']))==[lo,hi]and a['minimum_even_spin']==L and F(a['V_over_p_upper_exact'])<-k
            if a['method']=='low':assert F(a['spin_derivative_upper_exact'])<0
            na+=1
        f=r['finite_component'];assert f['spins']==list(range(42,L,2))
        if not f['spins']:continue
        p=inputs.read(f['producer']);q=inputs.read(f['replay']);assert q['replay_input_sha256']==sha(ROOT/f['producer'])and q['every_original_spin_interval_replayed']
        for j,bits in [(p,512),(q,640)]:
            assert j['status']=='PASS_FINITE_COMPACT_BRIDGE'and j['witness_sha256']==SHA and j['precision_bits']==bits
            assert j['spins']==f['spins']and set(map(int,j['cells_by_spin']))==set(f['spins'])and list(map(F,j['mu_domain']))==[lo,hi]
            for ell in f['spins']:
                cursor=lo
                for c in j['cells_by_spin'][str(ell)]:
                    a,b=F(c['a']),F(c['b']);assert cursor==a<b<=hi and F(c['V_over_p_upper_exact'])<-k;cursor=b;nc+=1
                assert cursor==hi
        for ell in f['spins']:
            assert [(c['a'],c['b'])for c in p['cells_by_spin'][str(ell)]]==[(c['a'],c['b'])for c in q['cells_by_spin'][str(ell)]]
    assert end==49
    audit=inputs.read(BASE/'compact/moments/independent_audit.json');assert audit['status']=='PASS_INDEPENDENT_EXACT_BETA_MOMENTS_AND_SHARED_TABLES'and audit['witness_sha256']==SHA
    return k,dict(paired_analytic_rectangles=na,paired_finite_spin_cells=nc)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.exists()
    inputs=Inputs();inputs.bind(WITNESS,SHA);w=np.load(WITNESS);direction=[str(F(float(x)))for x in w['direction']]
    assert len(w['names'])==len(w['multipliers'])==220
    definition=inputs.read(HEAD/'exact_witness_definition.json');assert definition['witness_sha256']==SHA and definition['direction_exact']==direction
    assert len(definition['coefficients'])==220
    for row,name,value in zip(definition['coefficients'],w['names'],w['multipliers']):
        assert row['row']==str(name)and F(row['exact_binary64_rational'])==F(float(value))and row['hex']==float(value).hex()
    h=head_cost(inputs,HEAD/'head_summary.json',SHA,direction)
    prep=inputs.read(HEAD/'preparation.json');assert prep['status']=='PASS_PREPARED_FINITE_HEAD_INPUTS'and prep['witness_sha256']==SHA
    inputs.bind(ROOT/'research/six_targets_20261007/upper86/prepare_head.py',prep['script_sha256'])
    for p,digest in prep['dependencies_sha256'].items():inputs.bind(p,digest)
    assert len(prep['source_rebuild_checks'])==6 and all(x['bit_for_bit_rebuild']for x in prep['source_rebuild_checks'])
    assert len(prep['commands'])==2 and all(c['exit_code']==0 for c in prep['commands'])
    campaign=inputs.read(HEAD/'campaign.json');jobs=inputs.read(HEAD/'jobs.json')
    assert campaign['status']=='COMPLETE_REPLAYED_FINITE_HEAD_ONLY'and campaign['witness_sha256']==SHA and campaign['jobs']==jobs
    inputs.bind(ROOT/'research/six_targets_20261007/upper86/run_head.py',campaign['script_sha256'])
    assert [j['spin']for j in jobs]==list(range(0,41,2))and campaign['assembly']['exit_code']==0
    for j in jobs:
        cmds=inputs.read(HEAD/f"spin{j['spin']}/commands.json");assert cmds==j['commands']and len(cmds)==2 and all(c['exit_code']==0 for c in cmds)
    summary=inputs.read(HEAD/'head_summary.json');budgets=inputs.read(HEAD/'budgets.json');assert budgets['witness_sha256']==SHA
    assert [r['spin']for r in budgets['rows']]==list(range(0,41,2))
    for r,b in zip(summary['records'],budgets['rows']):
        assert r['spin']==b['spin']
        assert r['model']==b['model']
    ck,checks=compact(inputs,direction)
    uv=inputs.read(BASE/'tail/uv_complete.json');assert uv['status']=='PASS_COMPLETE_DIRECTIONAL_UV_COERCIVITY'and uv['witness_sha256']==SHA and uv['objective_exact']==direction
    assert (uv['finite_first_even_spin'],uv['finite_last_even_spin'],uv['infinite_first_even_spin'])==(42,448,450)and uv['energy_domain']==['49','infinity']
    fk,ik=F(uv['finite_kappa_exact']),F(uv['infinite_kappa_exact']);assert min(fk,ik)>0
    algebra=inputs.read(BASE/'tail/exact_rebuild_and_algebra_audit.json')
    assert algebra['status']=='PASS_UV_EXACT_REBUILD_AND_ALGEBRA'and algebra['witness_sha256']==SHA and algebra['objective_exact']==direction
    assert algebra['all_entries_equal']and algebra['moment_count']==972 and algebra['remainder_coefficient_count']==9
    assert algebra['source_polynomial_equality_spins']==algebra['crossing_polynomial_equality_spins']==[450,452,600]and algebra['precision_bits']==512 and len(algebra['checks'])==18
    finite=real_norms(inputs,BASE/'real/finite_split_replay.json',SHA,'finite');infinite=real_norms(inputs,BASE/'real/tail_split_replay.json',SHA,'tail')
    realcheck=inputs.read(BASE/'real/independent_field_check.json');assert realcheck['status']=='PASS_INDEPENDENT_REAL_KERNEL_COMPARISON'and realcheck['witness_sha256']==SHA
    components=dict(head_0_40=h,finite_42_448_compact=finite['compact']/(2*ck),finite_42_448_UV_and_sliver=finite['remainder']/(2*min(ck,fk)),infinite_450_up_compact=infinite['compact']/(2*ck),infinite_450_up_UV_and_sliver=infinite['remainder']/(2*min(ck,ik)))
    total=sum(components.values(),F(0));outward=ceiling_decimal(total,9);dx,dy=map(F,direction);norm2=dx*dx+dy*dy;scale=2**256
    normlo=F(math.isqrt((norm2.numerator*scale*scale)//norm2.denominator),scale);assert normlo**2<=norm2<(normlo+F(1,scale))**2
    normalized=ceiling_decimal(F(outward)/normlo,9);angle=math.degrees(math.atan2(float(dy),float(dx)))%360
    for p in[Path(__file__),ROOT/'research/boundary_length_20261007/assemble_complete_v3.py',ROOT/'research/six_targets_20261007/REVIEW_UPPER86_ASSEMBLY.md']:
        inputs.bind(p)
    inputs.bind(ROOT/'research/six_targets_20261007/REVIEW_SIGNED_SOURCE2_UPPER86.md')
    inputs.bind(ROOT/'research/six_targets_20261007/tail_upper86/certify_angleupper86_large_tau_moments.py')
    result=dict(status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND',witness=str(WITNESS.relative_to(ROOT)),witness_sha256=SHA,objective_exact=direction,nominal_angle_degrees=angle,direction_norm_squared_exact=str(norm2),direction_norm_lower_exact=str(normlo),inequality=f'({direction[0]}) g0+({direction[1]}) g2 <= {outward}',support_upper_exact=str(total),support_upper_outward_9dp=outward,support_normalized_upper_outward_9dp=normalized,support_normalized_float=float(normalized),components_exact={k:str(v)for k,v in components.items()},components_float={k:float(v)for k,v in components.items()},tail_upper_exact=str(total-h),tail_upper_float=float(total-h),coercivity_constants_exact=dict(compact=str(ck),finite_UV=str(fk),infinite_UV=str(ik)),compact_margin_rechecks=checks,coverage='Every even spin and every energy above threshold, including source caps and both infinite-domain tails.',graphical_comparison=dict(dual=support_report(angle),primal=support_report(angle,'primal')),proof_scope='Internal same-witness interval certification, exact rational coverage and higher-precision replay. Shared mathematical implementation; no formal verification, external review or optimality claim.',sha256_manifest=inputs.manifest)
    a.out.write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k]for k in['status','support_normalized_float','components_float','coercivity_constants_exact','compact_margin_rechecks']},indent=2))
if __name__=='__main__':main()
