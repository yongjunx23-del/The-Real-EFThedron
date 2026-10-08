"""Fail-closed all-spin assembly for exact shifted 268-degree witnesses."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research/boundary_length_20261007'),str(ROOT/'research')]
from assemble_complete_v3 import Inputs,sha,ceiling_decimal,dyadic
from angle_benchmarks import support_report

def exact_partition(cells,gaps):
    whole=sorted(cells+gaps)
    assert whole[0][0]==0 and whole[-1][1]==1
    assert all(a<b for a,b in whole)
    assert all(x[1]==y[0]for x,y in zip(whole,whole[1:]))

def check_head(inputs,path,ell,bits,wp,base,direction,delta):
    d=inputs.read(path);assert d['status']=='PASS_EXACT_SHIFT_FINITE_HEAD_COMPONENT'
    assert d['spin']==ell and d['precision_bits']==bits
    complete=inputs.read(d['base_complete_certificate'],d['base_complete_sha256'])
    assert complete['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
    assert complete['witness_sha256']==d['base_witness_sha256']==base['witness_sha256']
    assert complete['objective_exact']==base['objective_exact']
    matches=[i for i,w in enumerate(d['witnesses'])if w['sha256']==sha(wp)]
    assert len(matches)==1;i=matches[0]
    inputs.bind(d['witnesses'][i]['path'],sha(wp))
    assert d['objective_exact'][i]==direction and F(d['g2_shift_exact'][i])==delta
    source=inputs.read(d['source_record'],d['source_record_sha256'])
    assert source['status']=='PASS_FINITE_SPIN_ENCLOSURE' and source['final']
    assert source['witness_sha256']==base['witness_sha256'] and source['spin']==ell
    assert source['direction_exact']==base['objective_exact']
    assert source['source_manifest_sha256']==d['source_manifest_sha256']
    summary_path=ROOT/'output/boundary_length_20261007/head/angle268/head_summary_final.json'
    summary=inputs.read(summary_path)
    ref=next(r for r in summary['records']if r['spin']==ell)
    assert ref['sha256']==d['source_record_sha256']
    assert summary['source_manifest_sha256']==d['source_manifest_sha256']
    inputs.read(summary['source_manifest_path'],d['source_manifest_sha256'])
    for p,h in d['dependencies'].items():inputs.bind(p,h)
    for p,h in source['dependencies'].items():inputs.bind(p,h)
    gaps=[tuple(map(F,x))for x in d['gap_intervals']]
    assert gaps==[tuple(map(F,x))for x in source['gap_beta_intervals']]
    assert F(d['base_gap_upper_exact'])==dyadic(source['gap_upper_dyadic'])>=0
    assert len(d['cells'])==len(source['cells_data'])
    cells=[];total=F(d['base_gap_upper_exact']);fallback=0
    valid={'base_monotonicity','shifted_riemann','shifted_convex_endpoints','shifted_negative_V'}
    for c,b in zip(d['cells'],source['cells_data']):
        lo,hi=F(c['a']),F(c['b']);assert (lo,hi)==(F(b['a']),F(b['b']))
        up=F(c['upper_exact'][i]);bu=dyadic(b['upper_dyadic'])
        assert 0<=up<=bu and c['methods'][i]in valid
        if c['methods'][i]=='base_monotonicity':assert up==bu;fallback+=1
        total+=up;cells.append((lo,hi))
    exact_partition(cells,gaps)
    assert d['partition_covers_beta_zero_to_one_exactly']
    assert total==F(d['upper_exact'][i])
    return total,dict(path=str(Path(path).relative_to(ROOT)),sha256=sha(path),precision_bits=bits,
        cells=len(cells),gaps=len(gaps),upper_exact=str(total),fallback_cells=fallback,elapsed_seconds=d['elapsed_seconds'])

def check_tail(inputs,path,wp,direction,delta):
    t=inputs.read(path);assert t['status']=='PASS_COMPLETE_TRANSFERRED_CONTINUUM_TAIL'
    assert t['witness_sha256']==sha(wp) and t['objective_exact']==direction and F(t['g2_shift_exact'])==delta
    assert t['first_even_spin']==42 and t['last_even_spin']=='infinity' and t['energy_domain']==['4','infinity']
    base=inputs.read(t['base_completed_certificate'],t['base_completed_certificate_sha256'])
    assert base['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND' and base['witness_sha256']==t['base_witness_sha256']
    old=np.load(inputs.bind(t['base_witness'],t['base_witness_sha256']));new=np.load(wp)
    assert np.array_equal(old['names'],new['names']) and np.array_equal(old['multipliers'],new['multipliers'])
    assert F(float(new['direction'][0]))==F(float(old['direction'][0]))
    assert F(float(new['direction'][1]))-F(float(old['direction'][1]))==delta<0
    folder=ROOT/'output/boundary_length_20261007/real/angle268'
    fin=inputs.read(folder/'finite_banded_cost.json');inf=inputs.read(folder/'tail_split_replay.json')
    assert fin['witness_sha256']==inf['witness_sha256']==base['witness_sha256']
    assert len(t['finite_tail_cells'])==len(fin['cells'])
    total=F(0);cells=[]
    for c,b in zip(t['finite_tail_cells'],fin['cells']):
        assert (c['lo'],c['hi'],c['region'])==(b['lo'],b['hi'],b['region'])
        A=F(c['A_upper_exact']);k=F(c['transferred_kappa_exact']);w=F(c['reciprocal_weight_upper_exact'])
        assert A==F(b['A_upper_exact']) and F(c['base_kappa_exact'])==F(b['kappa_exact'])
        assert k==F(b['kappa_exact'])-delta>0 and w>=1/(2*k)
        cost=A*w;assert cost==F(c['cost_upper_exact']);total+=cost;cells.append((F(c['lo']),F(c['hi'])))
    exact_partition(cells,[])
    ks=base['coercivity_constants_exact'];ck=F(ks['compact']);ik=F(ks['infinite_UV'])
    assert [x['region']for x in t['infinite_tail_regions']]==['compact','remainder']
    for c in t['infinite_tail_regions']:
        region=c['region'];k0=ck if region=='compact'else min(ck,ik)
        A=F(c['A_upper_exact']);k=F(c['transferred_kappa_exact']);w=F(c['reciprocal_weight_upper_exact'])
        assert A==F(inf['A_upper_exact'][region]) and F(c['base_kappa_exact'])==k0
        assert k==k0-delta>0 and w>=1/(2*k)
        cost=A*w;assert cost==F(c['cost_upper_exact']);total+=cost
    assert total==F(t['tail_upper_exact'])==sum(map(F,t['components_exact'].values()),F(0))
    return total,base

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--head-dir',type=Path,required=True)
    ap.add_argument('--tail',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    assert not a.out.exists();inputs=Inputs();wp=inputs.bind(a.witness);w=np.load(wp);direction=[str(F(float(x)))for x in w['direction']]
    initial=json.loads(a.tail.read_text());delta=F(initial['g2_shift_exact']);tail,base=check_tail(inputs,a.tail,wp,direction,delta)
    head=F(0);records=[]
    for ell in range(0,41,2):
        pair=[];evidence=[]
        for bits in [256,384]:
            value,rec=check_head(inputs,a.head_dir/f'head_l{ell}_{bits}.json',ell,bits,wp,base,direction,delta)
            pair.append(value);evidence.append(rec)
        retained=max(pair);head+=retained;records.append(dict(spin=ell,retained_upper_exact=str(retained),retained_upper_float=float(retained),precision_records=evidence))
    total=head+tail;outward=ceiling_decimal(total,9);dx,dy=map(F,direction);norm2=dx*dx+dy*dy
    scale=2**256;norm_lower=F(math.isqrt((norm2.numerator*scale*scale)//norm2.denominator),scale)
    assert 0<norm_lower and norm_lower*norm_lower<=norm2<(norm_lower+F(1,scale))**2
    normalized=ceiling_decimal(F(outward)/norm_lower,9)
    angle=math.degrees(math.atan2(float(dy),float(dx)))%360
    inputs.bind(Path(__file__));inputs.bind(Path(__file__).with_name('REVIEW_SHIFTED_HEAD.md'))
    result=dict(status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND',witness=str(wp.relative_to(ROOT)),witness_sha256=sha(wp),
        objective_exact=direction,nominal_angle_degrees=angle,direction_norm_squared_exact=str(norm2),direction_norm_lower_exact=str(norm_lower),
        inequality=f'({direction[0]}) g0 + ({direction[1]}) g2 <= {outward}',support_upper_exact=str(total),support_upper_outward_9dp=outward,
        support_normalized_upper_outward_9dp=normalized,support_normalized_float=float(F(normalized)),components_exact=dict(head_0_40=str(head),tail_42_infinity=str(tail)),
        components_float=dict(head_0_40=float(head),tail_42_infinity=float(tail)),tail_upper_exact=str(tail),tail_upper_float=float(tail),
        g2_shift_exact=str(delta),base_witness_sha256=base['witness_sha256'],finite_head_records=records,
        coverage='Every even spin ell>=0 and every energy mu>4, including threshold and source endpoint gaps, infinite energy and infinite spin.',
        graphical_comparison=dict(dual=support_report(angle),primal=support_report(angle,'primal')),
        proof_scope='Exact stored dyadic coefficients and objective shift; directed interval ordinary cells at256/384bits with exact monotonicity fallbacks and caps; exact all-spin tail transfer. Internal computer-assisted certificate, not formal verification, external review or optimality.',
        sha256_manifest=inputs.manifest)
    a.out.write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k]for k in ['status','nominal_angle_degrees','support_normalized_float','components_float']},indent=2))
if __name__=='__main__':main()
