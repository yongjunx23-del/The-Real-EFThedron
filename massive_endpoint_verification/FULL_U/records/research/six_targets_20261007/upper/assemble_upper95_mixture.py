"""Exact full-domain assembly of83+b6 with reintegrated low-spin supports."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math,bisect
import numpy as np
from assemble_upper55 import Inputs,headgeneric,upward
from upper95_common import *

def legacy_components(inputs,cert):
    base=ROOT/'output/competitive_20261007/upper_arc/b6'
    bundle=inputs.read(base/'head_bundle.json')
    assert bundle['status']=='PASS_FINITE_HEAD_0_40'and bundle['witness_sha256']==cert['witness_sha256']
    assert bundle['objective_exact']==cert['objective_exact']==['1','6']and F(bundle['crossing_increment_exact'])==F(107,5)
    parts={k:F(v)for k,v in bundle['component_upper_exact'].items()}
    assert sum(parts.values(),F(0))==F(bundle['finite_head_upper_exact'])==F(cert['components_exact']['finite_head0_40'])
    for ell in[0,2,4,6]:
        values=[];geometry=None
        for bits in[280,340]:
            p=base/f'head_p{bits}_spin{ell}.json';r=inputs.read(p);inputs.dependencies(r)
            assert r['status']=='PASS'and r['partition_exactly_covers_beta_0_1']and r['spin']==ell and r['precision_bits']==bits
            assert r['witness_sha256']==cert['witness_sha256']and F(r['repair_ratio_exact'])==F(107,30)
            inputs.bind(ROOT/'research/full_direction_20261007/head_family_extended.py',r['source_sha256'])
            assert len(r['directions'])==1;d=r['directions'][0]
            assert F(d['slope_exact'])==6 and F(d['crossing_increment_exact'])==F(107,5)
            cells=[(F(c['a']),F(c['b']))for c in d['cells']];gaps=[(F(c['a']),F(c['b']))for c in d['gaps']];partition(cells,gaps)
            if geometry is None:geometry=(cells,gaps)
            else:assert geometry==(cells,gaps)
            bounds=[F(c['upper_exact'])for c in d['cells']+d['gaps']];assert min(bounds)>=0
            value=sum(bounds,F(0));assert value==F(d['head_upper_exact']);values.append(value)
        assert max(values)==parts[str(ell)]
    finite=inputs.read(base/'finite_block.json')
    assert finite['status']=='PASS_FINITE_SUPPORT_8_40'and finite['spins']==list(range(8,42,2))
    assert finite['witness_sha256']==cert['witness_sha256']and finite['direction_exact']==['1','6']and F(finite['crossing_increment_exact'])==F(107,5)
    assert F(finite['finite_support_upper_exact'])==parts['8..40']
    return{ell:parts[str(ell)]for ell in[0,2]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();assert not args.out.exists();i=Inputs()
    certs=[i.read(p)for p in CERTS]
    for c in certs:
        assert c['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
        assert sum(map(F,c['components_exact'].values()),F(0))==F(c['support_upper_exact'])
    states=[np.load(i.bind(p,c['witness_sha256']))for p,c in zip(WITNESSES,certs)]
    assert list(map(lambda v:F(float(v)),states[0]['direction']))==list(map(F,certs[0]['objective_exact']))
    assert list(map(lambda v:F(float(v)),states[1]['direction']))==[F(1),F(0)]and F(certs[1]['crossing_increment_exact'])==F(107,5)
    target=np.load(i.bind(OUT/'objective.npz'));direction=[F(float(v))for v in target['direction']]
    weights=exact_weights(direction,[list(map(F,c['objective_exact']))for c in certs])
    h0,records=headgeneric(i,certs[0],'output/peak_certification_20261007/upper83/head/head_summary_selected.json')
    h1=legacy_components(i,certs[1]);costs=[h0,h1]
    prep=i.read(OUT/'preparation.json');assert prep['status']=='PASS_EXACT_LEGACY_SOURCE_CACHE_REBUILD'and prep['all_nested_exact_fields_equal']and prep['N']==64
    assert prep['witness_sha256']==certs[1]['witness_sha256']
    mp=ROOT/'output/peak_certification_20261007/upper83/head/moments128/manifest.json';manifest=i.read(mp)
    for table in manifest['tables']:
        for suffix,key in[('_moments.json','moment_sha256'),('_uv_moments.json','uv_sha256')]:i.bind(mp.parent/(table['label']+suffix),table[key])
    definition=[]
    for j,(wp,st)in enumerate(zip(WITNESSES,states)):
        definition.append(dict(base_witness=wp,base_witness_sha256=sha(ROOT/wp),names=st['names'].tolist(),base_multipliers_exact=[str(F(float(v)))for v in st['multipliers']],weight_exact=str(weights[j]),additional_X4_exact='0'if j==0 else'107/5'))
    retained=[];updates=[]
    for ell in[0,2]:
        raw=[i.read(p)for p in source_paths(ell)]
        cc,gg,gb=canonical_records(raw,ell);cells,gaps,gap_components,ledger=union_gap_data(cc,gg,gb)
        starts=[[a for a,b,u in group]for group in cc]
        pair=[i.read(OUT/f'head_l{ell}_{bits}.json')for bits in[256,384]]
        for r,bits in zip(pair,[256,384]):
            assert r['status']=='PASS_EXACT_POSITIVE_MIXTURE_HEAD_COMPONENT'and r['spin']==ell and r['precision_bits']==bits;i.dependencies(r)
            assert list(map(F,r['objective_exact']))==direction and list(map(F,r['weights_exact']))==weights
            assert r['base_certificate_paths']==CERTS and r['base_witness_sha256']==[c['witness_sha256']for c in certs]
            assert r['coefficient_definition']==definition and r['source_record_paths']==source_paths(ell)
            assert r['gap_cover_ledger']==ledger and list(map(F,r['gap_components_exact']))==gap_components
            assert[tuple(map(F,g))for g in r['gap_intervals']]==gaps
            assert[(F(c['a']),F(c['b']))for c in r['cells']]==cells
            assert F(r['gap_upper_exact'])==sum((w*g for w,g in zip(weights,gap_components)),F(0))
            assert all(F(c['upper_exact'])>=0 and c['method']in['convex_base_cells','mixture_riemann','mixture_convex_endpoints','mixture_negative_V']for c in r['cells'])
            assert all(c['method']!='mixture_negative_V'or(ell==2 and F(c['b'])**2<=F(79,119))for c in r['cells'])
            for cell in r['cells']:
                if cell['method']!='convex_base_cells':continue
                lo,hi=F(cell['a']),F(cell['b']);fallback=F(0)
                for j in[0,1]:
                    index=bisect.bisect_right(starts[j],lo)-1;a,b,u=cc[j][index]
                    assert a<=lo<hi<=b;fallback+=weights[j]*u
                assert F(cell['upper_exact'])==fallback
            assert F(r['upper_exact'])==F(r['gap_upper_exact'])+sum((F(c['upper_exact'])for c in r['cells']),F(0))
        old=sum(weights[j]*costs[j][ell]for j in[0,1]);direct=max(F(r['upper_exact'])for r in pair);keep=min(old,direct)
        retained.append(keep);updates.append(dict(spin=ell,old_convex_upper_exact=str(old),paired_direct_upper_exact=str(direct),retained_upper_exact=str(keep),cells=len(cells),gaps=len(gaps)))
    convex=sum((weights[j]*F(certs[j]['support_upper_exact'])for j in[0,1]),F(0))
    subtract=sum((weights[j]*sum((costs[j][ell]for ell in[0,2]),F(0))for j in[0,1]),F(0));other=convex-subtract;assert other>=0
    total=other+sum(retained,F(0));assert 0<=total<convex
    norm2=sum((v*v for v in direction),F(0));scale=2**256;normlo=F(math.isqrt(norm2.numerator*scale*scale//norm2.denominator),scale)
    assert 0<normlo and normlo**2<=norm2<(normlo+F(1,scale))**2
    outward=upward(total);normal=upward(F(outward)/normlo);angle=math.degrees(math.atan2(float(direction[1]),float(direction[0])))%360
    runs=i.read(OUT/'run_manifest.json');assert runs['status']=='COMPLETE_PAIRED_LOW_SPIN_MIXTURE'and len(runs['jobs'])==4 and all(j['exit_code']==0 for j in runs['jobs'])
    peer=i.read(OUT/'peer_field_and_gap_checks.json')
    assert peer['status'].startswith('PASS')
    for j in runs['jobs']:assert i.read(OUT/f"command_l{j['spin']}_{j['bits']}.json")==j
    i.bind(Path(__file__).with_name('run_upper95.py'),runs['script_sha256'])
    for p in[Path(__file__),Path(__file__).with_name('assemble_upper55.py'),Path(__file__).with_name('upper95_common.py'),ROOT/'research/six_targets_20261007/REVIEW_UPPER95_MIXTURE.md']:i.bind(p)
    full_definition=dict(witness_format='exact_rational_positive_combination_with_crossing_increment',objective_exact=list(map(str,direction)),base_certificate_paths=CERTS,weights_exact=list(map(str,weights)),components=definition)
    dp=OUT/'mixture_definition.json';assert not dp.exists();dp.write_text(json.dumps(full_definition,indent=2));i.bind(dp)
    sys_path=str(ROOT/'research')
    import sys
    if sys_path not in sys.path:sys.path.insert(0,sys_path)
    from angle_benchmarks import support_report
    result=dict(status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND',witness=str(dp.relative_to(ROOT)),witness_sha256=sha(dp),witness_format=full_definition['witness_format'],objective_exact=list(map(str,direction)),nominal_angle_degrees=angle,weights_exact=list(map(str,weights)),support_upper_exact=str(total),support_upper_outward_9dp=outward,support_normalized_upper_outward_9dp=normal,support_normalized_float=float(normal),direction_norm_squared_exact=str(norm2),direction_norm_lower_exact=str(normlo),inequality=f'({direction[0]})g0+({direction[1]})g2<={outward}',old_convex_bound_exact=str(convex),removed_components_exact=str(subtract),retained_other_spins_exact=str(other),refined_components=updates,components_exact=dict(all_other_spins_and_tails=str(other),spin0=str(retained[0]),spin2=str(retained[1])),coverage='Every even spin and every energy above threshold. Spins0,2 reintegrated as exact83+b6 mixture. All remaining head and tail domains inherited by positive homogeneity and subadditivity.',proof_scope='Internally reviewed directed intervals at256/384 bits and exact component accounting; source kernels shared, no external or formal verification.',graphical_comparison=dict(dual=support_report(angle),primal=support_report(angle,'primal')),sha256_manifest=i.manifest)
    args.out.write_text(json.dumps(result,indent=2));print(result['status'],normal,updates)
if __name__=='__main__':main()
