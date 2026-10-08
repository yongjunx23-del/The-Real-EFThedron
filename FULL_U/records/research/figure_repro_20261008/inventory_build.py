"""Export exactly the inequalities needed by the final two-panel figure.

This is an exact provenance/geometry/multiplier export, not an interval replay.
Frozen source files are read only. Run from any directory with --root.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse, csv, hashlib, json
import numpy as np

FIGURE = 'output/primal_dual_no_dagger_20261008/plot_data.json'
POLYGON = 'output/six_targets_20261007/comparison/certified_polygon.json'
BASE150 = 'output/angles_resummed/angle_000p000.npz'
LEFT220 = 'output/full_direction_20261007/left/angle180_final.npz'
G2220 = 'output/competitive_20261007/robust/g2_robust_stage1.npz'

def norm(p): return str(p).replace('\\', '/')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text(encoding='utf-8-sig'))
def decimal(f):
    with localcontext() as c:
        c.prec = 28
        return str(Decimal(f.numerator)/Decimal(f.denominator))
def exactrow(name, value, hx='', origin=''):
    f = F(value)
    return dict(name=name, exact=str(f), numerator=str(f.numerator), denominator=str(f.denominator),
                decimal_display=decimal(f), binary64_hex=hx, origin=origin)
def source_result(packet, objective):
    if 'results' not in packet: return packet, None
    matches = [(i,r) for i,r in enumerate(packet['results'])
               if list(map(F,r.get('direction_exact',['1',r.get('slope_exact','0')]))) == objective]
    assert len(matches)==1
    i,r=matches[0]
    return r,i
def result_objective(row):
    return list(map(F,row.get('objective_exact',row.get('direction_exact',['1',row.get('slope_exact','0')]))))
def constant(row): return F(row.get('support_upper_exact',row.get('certified_upper_exact')))
def components(row):
    if 'components_exact' in row: return row['components_exact']
    if 'heads_upper_exact' in row:
        return {**{'head_spin'+str(k):v for k,v in row['heads_upper_exact'].items()},
                'finite8_40':row['finite8_40_upper_exact'],'tail42_infinity':row['tail42_infinity_upper_exact']}
    return {}

def export(root, out):
    figure,poly=read(root/FIGURE),read(root/POLYGON)
    cuts=[tuple(map(F,[r['a_exact'],r['b_exact'],r['C_exact']])) for r in poly['cuts']]
    vv=[tuple(map(F,v)) for v in poly['vertices_exact']]
    assert all(a*x+b*y<=c for x,y in vv for a,b,c in cuts)
    edgecuts=[]
    for j,(x,y) in enumerate(vv):
        X,Y=vv[(j+1)%len(vv)]
        ids=[i for i,(a,b,c) in enumerate(cuts) if a*x+b*y==c and a*X+b*Y==c]
        assert ids, ('unbounded/unidentified edge',j)
        edgecuts.append(dict(edge_index=j,vertex_indices=[j,(j+1)%len(vv)],cut_indices=ids))
    top=sorted({i for e in edgecuts for i in e['cut_indices']})
    lower=sorted({r['cut_index'] for r in figure['displayed']} | {figure['minimum_g0_endpoint']['cut_index']})
    selected=sorted(set(top)|set(lower))
    seeds={FIGURE,POLYGON}
    raw={}; records=[]; multiplier_records=[]
    def npzwitness(path):
        path=norm(path);seeds.add(path)
        if path not in raw:
            z=np.load(root/path,allow_pickle=False)
            names=list(map(str,z['names'])); vals=[float(v) for v in z['multipliers']]
            assert len(names)==len(vals)==len(set(names))
            raw[path]=dict(path=path,sha256=sha(root/path),row_count=len(names),
                           objective_exact=list(map(lambda x:str(F.from_float(float(x))),z['direction'])),
                           rows=[exactrow(n,F.from_float(v),v.hex(),path) for n,v in zip(names,vals)])
        return raw[path]
    for i in selected:
        c=poly['cuts'][i];a,b,C=cuts[i]
        r=dict(cut_id=f'C{i:02}',cut_index=i,upper_panel_edge_active=i in top,lower_panel_displayed=i in lower,
               edge_indices=[e['edge_index']for e in edgecuts if i in e['cut_indices']],
               objective_exact=[str(a),str(b)],plotted_support_exact=str(C),
               lower_panel_markers=[x for x in figure['displayed']if x['cut_index']==i],
               minimum_g0_endpoint=figure['minimum_g0_endpoint'] if i==figure['minimum_g0_endpoint']['cut_index']else None)
        if str(c['source']).startswith('g2>='):
            assert (a,b,C)==(0,-1,0)
            r.update(proof_family='analytic positivity',source_certificate=None,exact_certified_support='0',
                     witness_recipe=dict(type='analytic',proof='Positive spectral weight p_l(mu) and Im f_l(mu)>=0 imply g2>=0.'),
                     bank_row_count=0,parents=[],selected_proof_ledger=dict(components_exact={},coverage='Every even spin and mu>4'))
            records.append(r);multiplier_records.append(dict(cut_id=r['cut_id'],rows=[]));continue
        path=norm(c['source']);seeds.add(path);packet=read(root/path);row,index=source_result(packet,[a,b])
        assert result_objective(row)==[a,b]
        Ce=constant(row); assert Ce<=C
        r.update(source_certificate=path,source_sha256=sha(root/path),source_result_index=index,
                 source_status=packet['status'],exact_certified_support=str(Ce),
                 plotted_rounding_slack_exact=str(C-Ce))
        parents=[]; effective=[]
        if 'results' in packet:
            w=npzwitness(BASE150);delta=F(row.get('certified_effective_crossing_increment_exact',row.get('crossing_increment_exact','0')))
            assert w['sha256']==packet['witness_sha256']
            for x in w['rows']:
                z=F(x['exact'])+(delta if x['name']=='X4'else 0)
                effective.append(exactrow(x['name'],z,x['binary64_hex'] if not(x['name']=='X4' and delta)else '',BASE150))
            mix=row.get('positive_combination')
            if mix:
                weights=[F(mix['left_weight_exact']),F(mix['right_weight_exact'])]
                assert all(x>=0 for x in weights) and sum(weights)==1
                for side,weight in zip(['left','right'],weights):
                    parent,pi=source_result(packet,[F(1),F(mix[side+'_slope_exact'])])
                    assert constant(parent)==F(mix[side+'_certified_constant_exact'])
                    parents.append(dict(certificate=path,result_index=pi,objective_exact=['1',mix[side+'_slope_exact']],
                                        weight_exact=str(weight),certified_constant_exact=str(constant(parent))))
                assert sum(F(p['weight_exact'])*F(p['objective_exact'][1])for p in parents)==b
                assert sum(F(p['weight_exact'])*F(p['certified_constant_exact'])for p in parents)==Ce
                assert sum(F(p['weight_exact'])*F(packet['results'][p['result_index']].get('certified_effective_crossing_increment_exact','0'))for p in parents)==delta
            r['proof_family']='150-row base with exact objective/X4 transformation'
            recipe=dict(type='npz_objective_and_X4_shift',base_witness=BASE150,base_sha256=w['sha256'],
                        objective_exact=[str(a),str(b)],X4_increment_exact=str(delta),positive_combination=mix)
        else:
            wp=norm(packet.get('witness',LEFT220 if '/left/'in path else G2220 if 'complete_g2_certificate'in path else ''))
            assert wp;seeds.add(wp)
            if wp.endswith('.npz'):
                w=npzwitness(wp);assert w['sha256']==packet['witness_sha256']
                assert list(map(F,w['objective_exact']))==[a,b]
                effective=w['rows'];recipe=dict(type='frozen_npz',witness=wp,witness_sha256=w['sha256'])
                r['proof_family']='220-row direct'
                if 'g2_shift_exact'in packet:
                    base='output/boundary_length_20261007/search/angle268.npz';bw=npzwitness(base)
                    assert bw['sha256']==packet['base_witness_sha256']
                    assert [(x['name'],x['exact'])for x in bw['rows']]==[(x['name'],x['exact'])for x in effective]
                    assert F(bw['objective_exact'][0])==a and F(bw['objective_exact'][1])+F(packet['g2_shift_exact'])==b
                    pp='output/boundary_length_20261007/complete_angle268_certificate.json';seeds.add(pp)
                    tp=wp.removesuffix('.npz')+'_tail_certificate.json';seeds.add(tp);tail=read(root/tp)
                    bp=norm(tail['base_completed_certificate']);seeds.add(bp)
                    parents=[dict(certificate=pp,witness=base,operation='same real field; base head used for monotonicity caps and fallbacks'),
                             dict(certificate=bp,witness=base,operation='same real norms; strengthened coercivity kappa_new=kappa_base-delta; tail costs recomputed')]
                    assert F(tail['tail_upper_exact'])==F(packet['tail_upper_exact'])
                    recipe.update(type='npz_exact_objective_shift',base_witness=base,g2_shift_exact=packet['g2_shift_exact'],
                                  transferred_tail_record=tp,transferred_tail_sha256=sha(root/tp))
                    r['proof_family']='220-row corrected head and monotone tail transfer'
            else:
                definition=read(root/wp);assert sha(root/wp)==packet['witness_sha256']
                assert list(map(F,definition['objective_exact']))==[a,b]
                effective=[exactrow(n,v,'',wp)for n,v in zip(definition['names'],definition['multipliers_exact'])]
                assert len(effective)==len(set(x['name']for x in effective))
                recipe=dict(type=definition.get('witness_format',definition.get('format')),witness=wp,
                            witness_sha256=packet['witness_sha256'],definition=definition)
                if 'base_witness_paths'in definition:
                    weights=list(map(F,definition['weights_exact']));assert all(x>=0 for x in weights)
                    bases=[npzwitness(p)for p in definition['base_witness_paths']]
                    assert all([x['name']for x in w['rows']]==[x['name']for x in effective]for w in bases)
                    assert [w['sha256']for w in bases]==definition['base_witness_sha256']
                    assert [sum(weight*F(w['objective_exact'][j])for weight,w in zip(weights,bases))for j in range(2)]==[a,b]
                    assert all(sum(weight*F(w['rows'][j]['exact'])for weight,w in zip(weights,bases))==F(x['exact'])for j,x in enumerate(effective))
                    for n,(bp,weight)in enumerate(zip(definition['base_certificate_paths'],weights)):
                        bp=norm(bp);seeds.add(bp);bpdoc=read(root/bp)
                        assert result_objective(bpdoc)==list(map(F,bases[n]['objective_exact']))
                        parents.append(dict(certificate=bp,sha256=sha(root/bp),witness=bases[n]['path'],weight_exact=str(weight),
                                            objective_exact=bases[n]['objective_exact'],certified_constant_exact=str(constant(bpdoc))))
                    assert sum(F(p['weight_exact'])*F(p['certified_constant_exact'])for p in parents)==F(packet['old_convex_bound_exact'])
                    assert F(packet['old_convex_bound_exact'])-F(packet['removed_components_exact'])==F(packet['retained_other_spins_exact'])
                    assert sum(F(x['old_convex_upper_exact'])for x in packet['refined_components'])==F(packet['removed_components_exact'])
                    assert F(packet['retained_other_spins_exact'])+sum(F(x['retained_upper_exact'])for x in packet['refined_components'])==Ce
                    r['proof_family']='220-row exact positive mixture with selected low spins reintegrated'
                else:
                    w=npzwitness(definition['base_witness']);delta=F(definition['crossing_increment_exact'])
                    assert w['sha256']==definition['base_witness_sha256']
                    assert [x['name']for x in w['rows']]==[x['name']for x in effective]
                    for x,rawx in zip(effective,w['rows']):
                        assert F(x['exact'])==F(rawx['exact'])+(delta if x['name']=='X4'else 0)
                        if F(x['exact'])==F(rawx['exact']):x['binary64_hex']=rawx['binary64_hex'];x['origin']=w['path']
                    r['proof_family']='150-row exact objective-shift family'
        comps=components(row)
        if not row.get('positive_combination'):assert sum(map(F,comps.values()),F())==Ce,(i,'component sum')
        # Keep selected row fields verbatim. Exclude only repeated bulk hash maps/plot comparisons.
        ledger={k:v for k,v in row.items()if k not in ['sha256_manifest','graphical_comparison']}
        ledger['effective_components_exact']=comps
        ledger['component_interpretation']='Positive combination of parent certified constants; direct component fields are historical provenance.' if row.get('positive_combination')else 'Exact additive sum of effective_components_exact.'
        if 'results'in packet:
            ledger['packet_scope']={k:v for k,v in packet.items()if k not in ['results','sha256_manifest']}
        r.update(witness_recipe=recipe,bank_row_count=len(effective),nonzero_multiplier_count=sum(F(x['exact'])!=0 for x in effective),
                 parents=parents,selected_proof_ledger=ledger)
        records.append(r);multiplier_records.append(dict(cut_id=r['cut_id'],objective_exact=r['objective_exact'],
                bank_row_count=len(effective),rows=effective))
    registry=dict(format='FINAL_FIGURE_EXACT_SELECTION_V1',status='PASS_EXACT_SELECTION_AND_MULTIPLIER_RECONSTRUCTION',
                  scope='All upper-panel edge-active inequalities, union all lower-panel displayed certificates. Parent proofs are dependencies; the full60-cut geometry remains preserved.',
                  source_figure=dict(path=FIGURE,sha256=sha(root/FIGURE)),source_polygon=dict(path=POLYGON,sha256=sha(root/POLYGON)),
                  upper_edge_active_cut_indices=top,lower_displayed_unique_cut_indices=lower,selected_cut_indices=selected,
                  selected_count=len(selected),numerical_certificate_count=sum(r['bank_row_count']>0 for r in records),
                  analytic_count=sum(r['bank_row_count']==0 for r in records),upper_vertex_count=len(vv),
                  lower_interior_marker_count=len(figure['displayed']),lower_distinct_certificate_count=len(lower),
                  lower_visual_marker_count=figure['visual_marker_count'],edge_cut_mapping=edgecuts,records=records,
                  verification_scope='Exact coefficient, source, mixture, component-sum and geometry checks only. No new optimization, interval replay, formal verification or independent mathematical proof.')
    out.mkdir(parents=True,exist_ok=True)
    def write(name,data):(out/name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    write('selection_registry.json',registry)
    write('multipliers_exact.json',dict(format='EXACT_MULTIPLIERS_ALL_ROWS_V1',convention='Binary64 source values interpreted exactly; mixtures and shifts performed over rationals. Decimal values are display only.',certificates=multiplier_records))
    write('raw_witnesses_exact.json',dict(witnesses=list(raw.values())))
    write('proof_seeds.json',dict(selected_cut_indices=selected,paths=sorted(seeds),
            note='These are closure seeds, not a complete package manifest. Follow explicit hashes, proof-record dependencies and source imports.'))
    write('full_polygon_geometry.json',poly)
    write('lower_panel_geometry.json',dict(displayed=figure['displayed'],minimum_g0_endpoint=figure['minimum_g0_endpoint'],
               approximate_continuation=figure['approximate_continuation'],connecting_lines_certified=False,
               normalization=figure['normalization'],branch_join=figure['branch_join']))
    fields=['cut_id','row_index','name','numerator','denominator','decimal_display','binary64_hex','origin']
    with (out/'multipliers.csv').open('w',newline='',encoding='utf-8')as f:
        wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader()
        for rec in multiplier_records:
            for j,x in enumerate(rec['rows']):wr.writerow(dict(cut_id=rec['cut_id'],row_index=j,**{k:x[k]for k in fields[2:]}))
    with (out/'raw_witness_multipliers.csv').open('w',newline='',encoding='utf-8')as f:
        fields2=['witness_path','row_index','name','numerator','denominator','decimal_display','binary64_hex']
        wr=csv.DictWriter(f,fieldnames=fields2);wr.writeheader()
        for rec in raw.values():
            for j,x in enumerate(rec['rows']):wr.writerow(dict(witness_path=rec['path'],row_index=j,**{k:x[k]for k in fields2[2:]}))
    with (out/'bounds.csv').open('w',newline='',encoding='utf-8')as f:
        keys=['cut_id','upper_panel_edge_active','lower_panel_displayed','a_exact','b_exact','plotted_support_exact',
              'certified_support_exact','certified_support_decimal_display','bank_row_count','source_certificate','source_result_index']
        wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader()
        for r in records:
            wr.writerow(dict(cut_id=r['cut_id'],upper_panel_edge_active=r['upper_panel_edge_active'],lower_panel_displayed=r['lower_panel_displayed'],
                 a_exact=r['objective_exact'][0],b_exact=r['objective_exact'][1],plotted_support_exact=r['plotted_support_exact'],
                 certified_support_exact=r['exact_certified_support'],certified_support_decimal_display=decimal(F(r['exact_certified_support'])),
                 bank_row_count=r['bank_row_count'],source_certificate=r['source_certificate'],source_result_index=r.get('source_result_index')))
    with (out/'support_components.csv').open('w',newline='',encoding='utf-8')as f:
        keys=['cut_id','interpretation','component','numerator','denominator','decimal_display']
        wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader()
        for r in records:
            ledger=r['selected_proof_ledger'];mix=ledger.get('positive_combination')
            if mix:
                cc={f"parent_{j}_weight_times_bound":str(F(p['weight_exact'])*F(p['certified_constant_exact']))for j,p in enumerate(r['parents'])}
            else:cc=ledger.get('effective_components_exact',{})
            for name,value in cc.items():
                v=F(value);wr.writerow(dict(cut_id=r['cut_id'],interpretation='selected exact additive contribution',component=name,
                      numerator=v.numerator,denominator=v.denominator,decimal_display=decimal(v)))
    (out/'README.md').write_text('''# Exact data for the final two-panel figure

The proof selection is the union of all edge-active upper-panel cuts and all lower-panel certificate markers. It contains **39 inequalities:38 numerical certificates plus the analytic inequality g2>=0**. All22 distinct lower-panel certificates already belong to the39 upper-panel cuts. The21 interior markers plus two copies of the minimum-g0 endpoint are23 visible circles. The full original60-cut polygon is retained only as a geometry/provenance input; redundant cuts are not additional selected bounds.

## Files

- `selection_registry.json`: authoritative exact selection, objective, plotted versus unrounded certified support, panel membership, original source/result identity, witness recipe, necessary parent certificates, selected proof ledger, and exact edge correspondence.
- `bounds.csv`: compact scalar summary; every row is the halfplane a*g0+b*g2<=C. The exact certified C may be slightly below the outward-rounded C drawn in the figure.
- `multipliers_exact.json` and `multipliers.csv`: all7030 effective multiplier rows for the38 numerical certificates, including every zero. Do not drop the row names or change their order.
- `raw_witnesses_exact.json` and `raw_witness_multipliers.csv`:15 frozen NPZ witnesses, their exact binary64 fractions and hexadecimal representations. Effective mixture coefficients are exact rational linear combinations, not rounded binary64 values.
- `support_components.csv`: exact additive contributions to each selected bound. C06 uses weighted parent constants; its historical direct-component fields do not sum to the selected bound.
- `full_polygon_geometry.json`: unchanged full60-cut geometry with exact rational vertices.
- `lower_panel_geometry.json`: exact representative positions, normal intersections and display normalization. Interpolating connections are not certificates, and the minimum-g0 endpoint uses extremal-g0 normalization.
- `proof_seeds.json`: dependency-collection seeds, **not** a self-contained manifest. The complete package also needs all recursively bound source moments, interval records, source code and pinned arithmetic dependencies.

## Exact arithmetic and interpretation

Numerator/denominator fields are integers encoded as text; use arbitrary-precision rationals. Decimal columns are for reading only. An available `binary64_hex` means the effective value itself equals that exact binary64 value; when blank, reconstruct it from the rational fields and recorded parents. Parent NPZ bit patterns remain available in the raw-witness files. The stored direction is an exact rational pair, not the transcendental cosine/sine of a rounded degree label.

All selected150-row certificates have zero X4 increment and share the same base multipliers; their objectives differ. C06 is the exact positive combination of the slope−1/2 and slope0 parent bounds. The modern220-row mixture certificates C51,C55,C57,C58,C59 retain necessary parent proof components and replace specified low-spin contributions with direct interval enclosures. C39 uses the corrected shifted head and a strengthened transferred tail: kappa_new=kappa_base−delta, not an unchanged tail constant.

The points are locations on certified outer halfplanes/polygon edges. They are not witnesses of physically feasible amplitudes, and certification does not establish optimizer optimality. Figure-reference coordinates and approximate interpolations do not enter the physical certificate proof.

## Read-only recheck

From the root of the extracted figure package, run:

```text
python inventory_verify.py --root records --data data
```

In the original workspace the equivalent command is `python research/figure_repro_20261008/inventory_verify.py --root .`. The verifier checks source hashes, exact objectives and every effective/parent coefficient, component sums and positive mixtures, C39's transferred tail arithmetic, both CSV/JSON representations, panel membership and exact reconstruction of all39 polygon vertices. It writes no files. It requires NumPy to read the original frozen NPZs. The exported rational JSON/CSV can also be read with other arbitrary-precision implementations.

This is a saved-data and rational-assembly verification, not a new interval-kernel replay. The appendix and package replay instructions describe the additional directed-Arb recomputation and mathematical assumptions. Optimizer search records are provenance; changing a witness invalidates its existing interval records until recertified.
''',encoding='utf-8')
    print(json.dumps(dict(status=registry['status'],selected=len(selected),numerical=registry['numerical_certificate_count'],
          analytic=registry['analytic_count'],witness_npz_count=len(raw),multiplier_rows=sum(len(r['rows'])for r in multiplier_records),seeds=len(seeds))))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);ap.add_argument('--out',type=Path)
    a=ap.parse_args();root=a.root.resolve();export(root,a.out.resolve()if a.out else root/'output/FIGURE_REPRO_20261008/data')
