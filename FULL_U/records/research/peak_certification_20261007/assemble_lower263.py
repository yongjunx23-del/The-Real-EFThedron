"""Fail-closed same-witness assembly of audited optimized 263.274737 degree bound."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research/boundary_length_20261007'),str(ROOT/'research')]
from assemble_complete_v3 import Inputs as BaseInputs,sha,head_cost,real_norms,ceiling_decimal
from banded_lower263 import Coercivity,dyadic_ceiling
from angle_benchmarks import support_report
BASE=ROOT/'output/peak_certification_20261007/lower263'
HEAD=BASE/'head'
WITNESS=ROOT/'output/weak_regions_20261007/lower/search/angle263p2747_from261.npz'
SHA='8931658c2a353a2bbf9b9cb074db490ba5d15f044c80837236d3fd9dd801c541'

class Inputs(BaseInputs):
    """Same hash bindings, also permitting logged command lists as JSON."""
    def read(self,path,expected=None):
        p=self.bind(path,expected);d=json.loads(p.read_text())
        if p not in self.seen:
            self.seen.add(p)
            for sub,h in (d.get('sha256_manifest',{}).items()if isinstance(d,dict)else []):
                q=self.bind(sub,h)
                if q.suffix=='.json':self.read(q,h)
        return d

def stronger_compact(inputs,direction):
    old=inputs.read(BASE/'compact/complete_compact_coercivity.json');new=inputs.read(BASE/'compact/complete_compact_coercivity_stronger.json')
    for d in [old,new]:
        assert d['status']=='PASS_COMPLETE_COMPACT_COERCIVITY' and d['witness_sha256']==SHA and d['objective_exact']==direction
        assert d['minimum_even_spin']==42 and d['mu_domain_exact']==['4','49'] and d['lower_energy_endpoint_open']
    assert new['coverage_rectangles']==old['coverage_rectangles']
    k=F(new['kappa_exact']);assert k>=F(old['kappa_exact'])>0
    cursor=F(4);analytic=finite=0
    for r in old['coverage_rectangles']:
        lo,hi=map(F,r['mu_domain_exact']);L=r['analytic_minimum_even_spin'];assert cursor==lo<hi;cursor=hi
        assert len(r['analytic_pair'])==2
        for file,bits in zip(r['analytic_pair'],[400,512]):
            j=inputs.read(file);assert j['witness_sha256']==SHA and j['objective_exact']==direction
            assert j['status']=='PASS_PIECEWISE_ALLSPIN_COMPACT_ENVELOPES' and j['precision_bits']==bits
            rr=next(x for x in j['records']if x['mu_domain_exact']==r['analytic_record'])if 'analytic_record'in r else j['records'][0]
            aa,bb=map(F,rr['mu_domain_exact']);assert aa<=lo<hi<=bb and rr['minimum_even_spin']==L
            assert F(rr['V_over_p_upper_exact'])<-k
            if rr['method']=='low':assert F(rr['spin_derivative_upper_exact'])<0
            analytic+=1
        f=r.get('finite_component',{'spins':[]});assert f['spins']==list(range(42,L,2))
        if not f['spins']:continue
        pp=inputs.read(f['producer']);rp=inputs.read(f['replay']);assert rp['replay_input_sha256']==sha(ROOT/f['producer'])
        for d,bits in [(pp,512),(rp,640)]:
            assert d['status']=='PASS_FINITE_COMPACT_BRIDGE' and d['witness_sha256']==SHA and d['precision_bits']==bits
            assert d['spins']==f['spins'] and list(map(F,d['mu_domain']))==[lo,hi]
            for ell in f['spins']:
                cells=sorted(d['cells_by_spin'][str(ell)],key=lambda x:F(x['a']));end=lo
                for c in cells:
                    a,b=F(c['a']),F(c['b']);assert end==a<b<=hi and F(c['V_over_p_upper_exact'])<-k;end=b;finite+=1
                assert end==hi
        for ell in f['spins']:
            assert [(c['a'],c['b'])for c in pp['cells_by_spin'][str(ell)]]==[(c['a'],c['b'])for c in rp['cells_by_spin'][str(ell)]]
    assert cursor==49
    return k,dict(analytic_rectangles_checked=analytic,finite_producer_and_replay_cells_checked=finite)

def local_cost(inputs,direction):
    path=BASE/'real/finite_banded_cost.json';d=inputs.read(path);real=inputs.read(BASE/'real/finite_split_replay.json')
    co=Coercivity('lower263');assert co.SHA==SHA and co.objective==direction
    for p,h in co.manifest.items():inputs.bind(p,h)
    assert d['status']=='PASS_BANDED_COERCIVITY_REAL_COST' and d['witness_sha256']==SHA and d['objective_exact']==direction
    assert d['first_even_spin']==42 and d['last_even_spin']==448 and d['energy_domain']==['4','infinity']
    assert len(d['cells'])==len(real['cells']);total=F(0)
    for c,r in zip(d['cells'],real['cells']):
        lo,hi=F(c['lo']),F(c['hi']);assert (lo,hi)==(F(r['lo']),F(r['hi'])) and c['region']==r['region']
        A=F(c['A_upper_exact']);assert A==F(r['retained_upper_exact'])
        cover=[(ell,co.coverage(ell,lo,hi))for ell in range(42,450,2)]
        limiting,entry=min(cover,key=lambda x:x[1][0]);k=entry[0]
        assert F(c['kappa_exact'])==k>0 and c['limiting_spin']==limiting
        assert c['compact_piece_included']==entry[1] and c['spins_checked']==204
        assert c['large_tau_spins']==[ell for ell,x in cover if x[2]]
        assert c['finite_tau_cell_index_ranges']==[x[3]for ell,x in cover if x[3]is not None]
        weight=F(c['reciprocal_weight_upper_exact']);assert weight==dyadic_ceiling(1/(2*k))>=1/(2*k)
        cost=A*weight;assert cost==F(c['cost_upper_exact']);total+=cost
    assert total==F(d['cost_upper_exact']);inputs.bind(ROOT/'research/peak_certification_20261007/banded_lower263.py')
    return total

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.exists()
    inputs=Inputs();inputs.bind(WITNESS,SHA);w=np.load(WITNESS);direction=[str(F(float(x)))for x in w['direction']]
    assert len(w['names'])==len(w['multipliers'])==220
    h=head_cost(inputs,HEAD/'head_summary.json',SHA,direction)
    prep=inputs.read(HEAD/'preparation.json');assert prep['status']=='PASS_PREPARED_FINITE_HEAD_INPUTS' and prep['witness_sha256']==SHA
    inputs.bind(ROOT/'research/peak_certification_20261007/lower263/certify_head.py',prep['script_sha256'])
    for p,digest in prep['dependencies_sha256'].items():inputs.bind(p,digest)
    assert len(prep['source_rebuild_checks'])==6 and all(x['bit_for_bit_rebuild']for x in prep['source_rebuild_checks'])
    budgets=inputs.read(HEAD/'budgets.json');assert budgets['witness_sha256']==SHA
    assert [r['spin']for r in budgets['rows']]==list(range(0,41,2))
    summary=inputs.read(HEAD/'head_summary.json')
    for r,b in zip(summary['records'],budgets['rows']):assert inputs.read(r['path'])['model']==b['model'] and r['spin']==b['spin']
    campaign=inputs.read(HEAD/'campaign.json');jobs=inputs.read(HEAD/'jobs.json')
    assert campaign['status']=='COMPLETE_REPLAYED_FINITE_HEAD_ONLY' and campaign['witness_sha256']==SHA
    assert campaign['jobs']==jobs and [r['spin']for r in jobs]==list(range(0,41,2))
    assert all(c['exit_code']==0 for c in campaign['commands'])
    for row in jobs:
        commands=inputs.read(HEAD/f"spin{row['spin']}/commands.json")
        assert commands==row['commands'] and len(commands)==2 and all(c['exit_code']==0 for c in commands)
    inputs.bind(ROOT/'research/peak_certification_20261007/lower263/resume_relaxed_head.py',campaign['script_sha256'])
    change=campaign['orchestration_change'];switch=inputs.read(change['path'],change['sha256'])
    assert change['old_preparation_and_proof_records_preserved'] and change['active_workers_adopted']==[6,8,10]
    assert len(switch['active_workers'])==3 and change['new_worker_limit']==6
    relaxed=inputs.read(campaign['budgets_path'],campaign['budgets_sha256']);assert relaxed['witness_sha256']==SHA
    assert [r['spin']for r in relaxed['rows']]==list(range(0,41,2))
    for old,row,rec in zip(budgets['rows'],relaxed['rows'],summary['records']):
        assert old['spin']==row['spin']==rec['spin'] and old['model']==row['model']==rec['model']
        assert old['numerical']==row['numerical'] and F(row['target'])==F(row['numerical'])+F(row['slack'])
        if row['spin']<=10:assert row['target']==old['target'] and row['max_cells']==7000
        else:assert F(row['slack'])==F(1,100000) and row['max_cells']==3000
        assert inputs.read(rec['path'])['target']==row['target']
    peer=inputs.read(BASE/'tail/peer_coverage_review.json')
    assert peer['status']=='PASS_LOWER_TAIL_PEER_SOURCE_COVERAGE_AND_ACCOUNTING' and peer['witness_sha256']==SHA
    assert F(peer['cost_upper_exact'])==F(inputs.read(BASE/'real/finite_banded_cost.json')['cost_upper_exact'])
    compact_review=inputs.read(BASE/'compact/parent_review.json');assert compact_review['status'].startswith('PASS')
    ck,compact_checks=stronger_compact(inputs,direction)
    uv=inputs.read(BASE/'tail/uv_complete.json');assert uv['status']=='PASS_COMPLETE_DIRECTIONAL_UV_COERCIVITY'
    assert uv['witness_sha256']==SHA and uv['objective_exact']==direction
    assert (uv['finite_first_even_spin'],uv['finite_last_even_spin'],uv['infinite_first_even_spin'])==(42,448,450)
    assert uv['energy_domain']==['49','infinity'];fk,ik=F(uv['finite_kappa_exact']),F(uv['infinite_kappa_exact']);assert min(fk,ik)>0
    algebra=inputs.read(BASE/'tail/exact_rebuild_and_algebra_audit.json')
    assert algebra['status']=='PASS_UV_EXACT_REBUILD_AND_ALGEBRA' and algebra['witness_sha256']==SHA and algebra['objective_exact']==direction
    assert algebra['all_entries_equal'] and algebra['moment_count']==972 and algebra['remainder_coefficient_count']==9
    assert algebra['source_polynomial_equality_spins']==algebra['crossing_polynomial_equality_spins']==[450,452,600]
    assert algebra['precision_bits']==512 and len(algebra['checks'])==18
    finite=real_norms(inputs,BASE/'real/finite_split_replay.json',SHA,'finite');infinite=real_norms(inputs,BASE/'real/tail_split_replay.json',SHA,'tail')
    local=local_cost(inputs,direction);global_finite=finite['compact']/(2*ck)+finite['remainder']/(2*min(ck,fk));assert local<=global_finite
    check=inputs.read(BASE/'real/independent_field_check.json');assert check['status']=='PASS_INDEPENDENT_REAL_KERNEL_COMPARISON' and check['witness_sha256']==SHA
    components=dict(head_0_40=h,finite_42_448_banded_all_energy=local,infinite_450_up_compact=infinite['compact']/(2*ck),infinite_450_up_UV_and_sliver=infinite['remainder']/(2*min(ck,ik)))
    total=sum(components.values(),F(0));outward=ceiling_decimal(total,9);dx,dy=map(F,direction);norm2=dx*dx+dy*dy;scale=2**256
    normlo=F(math.isqrt((norm2.numerator*scale*scale)//norm2.denominator),scale);assert normlo**2<=norm2<(normlo+F(1,scale))**2
    normalized=ceiling_decimal(F(outward)/normlo,9);angle=math.degrees(math.atan2(float(dy),float(dx)))%360
    inputs.bind(Path(__file__));inputs.bind(ROOT/'research/boundary_length_20261007/assemble_complete_v3.py')
    inputs.bind(ROOT/'research/boundary_length_20261007/REVIEW_BANDED_FINITE_COST.md')
    inputs.bind(ROOT/'research/peak_certification_20261007/REVIEW_LOWER263_ASSEMBLY.md')
    result=dict(status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND',witness=str(WITNESS.relative_to(ROOT)),witness_sha256=SHA,objective_exact=direction,nominal_angle_degrees=angle,
        direction_norm_squared_exact=str(norm2),direction_norm_lower_exact=str(normlo),inequality=f'({direction[0]}) g0+({direction[1]}) g2 <= {outward}',
        support_upper_exact=str(total),support_upper_outward_9dp=outward,support_normalized_upper_outward_9dp=normalized,support_normalized_float=float(normalized),
        components_exact={k:str(v)for k,v in components.items()},components_float={k:float(v)for k,v in components.items()},tail_upper_exact=str(total-h),tail_upper_float=float(total-h),
        coercivity_constants_exact=dict(compact=str(ck),finite_UV=str(fk),infinite_UV=str(ik)),compact_margin_rechecks=compact_checks,
        global_finite_tail_upper_exact=str(global_finite),finite_tail_method='Local exact energy-cell coercivity, using strengthened compact margin.',
        coverage='Every even spin ell>=0 and every energy mu>4, including threshold/source endpoint caps and infinite energy/spin.',
        graphical_comparison=dict(dual=support_report(angle),primal=support_report(angle,'primal')),
        proof_scope='Internal same-witness directed-interval certificate with higher-precision replays, exact rational coverage and stronger retained compact margin. Not formal verification, external review or optimality.',
        sha256_manifest=inputs.manifest)
    a.out.write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k]for k in ['status','support_normalized_float','components_float','coercivity_constants_exact','compact_margin_rechecks']},indent=2))
if __name__=='__main__':main()
