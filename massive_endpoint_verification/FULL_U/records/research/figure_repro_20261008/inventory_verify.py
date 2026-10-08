"""Read-only check of the final-figure exact selection and coefficient exports.

Usage: python inventory_verify.py --root EXTRACTED_WORKSPACE [--data DATA_DIR]
This does not call an optimizer or evaluate interval kernels.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse, csv, hashlib, json
import numpy as np

def n(p):return str(p).replace('\\','/')
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def objective(r):return list(map(F,r.get('objective_exact',r.get('direction_exact',['1',r.get('slope_exact','0')]))))
def support(r):return F(r.get('support_upper_exact',r.get('certified_upper_exact')))
def select(p,a):
    if 'results'not in p:return p,None
    matches=[(i,x)for i,x in enumerate(p['results'])if objective(x)==a]
    assert len(matches)==1
    i,x=matches[0];return x,i
def fields(p):
    if 'components_exact'in p:return p['components_exact']
    if 'heads_upper_exact'in p:
        return {**{'head_spin'+str(k):v for k,v in p['heads_upper_exact'].items()},
                'finite8_40':p['finite8_40_upper_exact'],'tail42_infinity':p['tail42_infinity_upper_exact']}
    return {}
def decimal(f):
    with localcontext()as c:c.prec=28;return str(Decimal(f.numerator)/Decimal(f.denominator))

def verify(root,data):
    reg=read(data/'selection_registry.json')
    coeff=read(data/'multipliers_exact.json'); exported={x['cut_id']:x for x in coeff['certificates']}
    raw=read(data/'raw_witnesses_exact.json');raw_export={n(x['path']):x for x in raw['witnesses']}
    cache={}
    def npz(p):
        p=n(p)
        if p not in cache:
            z=np.load(root/p,allow_pickle=False)
            names=list(map(str,z['names']));values=[F.from_float(float(v))for v in z['multipliers']]
            hx=[float(v).hex()for v in z['multipliers']];obj=[F.from_float(float(v))for v in z['direction']]
            assert len(names)==len(set(names))==len(values)
            saved=raw_export[p];assert digest(root/p)==saved['sha256']
            assert saved['row_count']==len(names)and list(map(F,saved['objective_exact']))==obj
            assert names==[r['name']for r in saved['rows']]
            assert values==[F(r['exact'])for r in saved['rows']]
            assert hx==[r['binary64_hex']for r in saved['rows']]
            cache[p]=(names,values,hx,obj)
        return cache[p]
    figref,polyref=reg['source_figure'],reg['source_polygon']
    assert digest(root/figref['path'])==figref['sha256']
    assert digest(root/polyref['path'])==polyref['sha256']
    fig,poly=read(root/figref['path']),read(root/polyref['path'])
    assert read(data/'full_polygon_geometry.json')==poly
    cuts=[tuple(map(F,[x['a_exact'],x['b_exact'],x['C_exact']]))for x in poly['cuts']]
    vertices=[tuple(map(F,v))for v in poly['vertices_exact']]
    assert len(set(vertices))==len(vertices)
    assert all(a*x+b*y<=c for x,y in vertices for a,b,c in cuts)
    edges=[]
    for j,(x,y)in enumerate(vertices):
        X,Y=vertices[(j+1)%len(vertices)]
        ids=[i for i,(a,b,c)in enumerate(cuts)if a*x+b*y==c and a*X+b*Y==c]
        assert ids
        edges.append(dict(edge_index=j,vertex_indices=[j,(j+1)%len(vertices)],cut_indices=ids))
    top=sorted({i for e in edges for i in e['cut_indices']})
    bottom=sorted({r['cut_index']for r in fig['displayed']}|{fig['minimum_g0_endpoint']['cut_index']})
    union=sorted(set(top)|set(bottom))
    assert reg['upper_edge_active_cut_indices']==top
    assert reg['lower_displayed_unique_cut_indices']==bottom
    assert reg['selected_cut_indices']==union
    assert reg['edge_cut_mapping']==edges
    assert [r['cut_index']for r in reg['records']]==union
    assert len(exported)==len(coeff['certificates'])==len(union)
    assert set(exported)=={f'C{i:02}'for i in union}
    assert len(raw_export)==len(raw['witnesses'])
    # All finite pair intersections of the selected halfplanes give precisely
    # the same rational vertex set as the full60-cut polygon.
    picked=[cuts[i]for i in union];rebuilt=set()
    for i,(a,b,c)in enumerate(picked):
        for A,B,C in picked[i+1:]:
            determinant=a*B-b*A
            if not determinant:continue
            p=((c*B-b*C)/determinant,(a*C-c*A)/determinant)
            if all(u*p[0]+v*p[1]<=w for u,v,w in picked):rebuilt.add(p)
    assert rebuilt==set(vertices)
    all_csv=[];seed_needed={figref['path'],polyref['path']}
    for r in reg['records']:
        i=r['cut_index'];a,b,C=cuts[i];rec=exported[r['cut_id']]
        assert r['cut_id']==f'C{i:02}'and r['objective_exact']==[str(a),str(b)]and F(r['plotted_support_exact'])==C
        assert r['upper_panel_edge_active']==(i in top)and r['lower_panel_displayed']==(i in bottom)
        assert r['edge_indices']==[e['edge_index']for e in edges if i in e['cut_indices']]
        assert r['lower_panel_markers']==[x for x in fig['displayed']if x['cut_index']==i]
        expected_endpoint=fig['minimum_g0_endpoint']if i==fig['minimum_g0_endpoint']['cut_index']else None
        assert r['minimum_g0_endpoint']==expected_endpoint
        if r['source_certificate']is None:
            assert (a,b,C)==(0,-1,0)and r['bank_row_count']==0 and rec['rows']==[]
            assert r['exact_certified_support']=='0'
            continue
        path=n(r['source_certificate']);seed_needed.add(path)
        assert path==n(poly['cuts'][i]['source'])and digest(root/path)==r['source_sha256']
        packet=read(root/path);row,index=select(packet,[a,b])
        assert index==r['source_result_index']and objective(row)==[a,b]
        Ce=support(row)
        assert F(r['exact_certified_support'])==Ce<=C
        assert F(r['plotted_rounding_slack_exact'])==C-Ce
        assert r['source_status']==packet['status']
        savedledger=r['selected_proof_ledger']
        for key,value in row.items():
            if key not in ['sha256_manifest','graphical_comparison']:assert savedledger[key]==value
        assert savedledger['effective_components_exact']==fields(row)
        recipe=r['witness_recipe'];parents=r['parents']
        if 'results'in packet:
            base='output/angles_resummed/angle_000p000.npz';seed_needed.add(base)
            names,vals,hx,baseobj=npz(base)
            assert digest(root/base)==packet['witness_sha256']
            delta=F(row.get('certified_effective_crossing_increment_exact',row.get('crossing_increment_exact','0')))
            assert F(recipe['X4_increment_exact'])==delta and recipe['base_witness']==base
            vals=[v+(delta if name=='X4'else 0)for name,v in zip(names,vals)]
            if row.get('positive_combination'):
                mix=row['positive_combination'];assert recipe['positive_combination']==mix
                total=F();direction=F();repair=F();weights=[]
                assert len(parents)==2
                for parent,side in zip(parents,['left','right']):
                    weight=F(mix[side+'_weight_exact']);weights.append(weight)
                    prow,pi=select(packet,[F(1),F(mix[side+'_slope_exact'])])
                    assert parent['certificate']==path and parent['result_index']==pi
                    assert F(parent['weight_exact'])==weight and F(parent['certified_constant_exact'])==support(prow)
                    assert support(prow)==F(mix[side+'_certified_constant_exact'])
                    total+=weight*support(prow);direction+=weight*objective(prow)[1]
                    repair+=weight*F(prow.get('certified_effective_crossing_increment_exact',prow.get('crossing_increment_exact','0')))
                assert min(weights)>=0 and sum(weights)==1 and total==Ce and direction==b and repair==delta
        else:
            wp=n(packet.get('witness','output/full_direction_20261007/left/angle180_final.npz'if '/left/'in path else 'output/competitive_20261007/robust/g2_robust_stage1.npz'))
            seed_needed.add(wp)
            assert digest(root/wp)==packet['witness_sha256']
            if wp.endswith('.npz'):
                names,vals,hx,obj=npz(wp);assert obj==[a,b]
                assert n(recipe['witness'])==wp
                if 'g2_shift_exact'in packet:
                    base=n(recipe['base_witness']);seed_needed.add(base)
                    bn,bv,bh,bo=npz(base)
                    assert bn==names and bv==vals and bo[0]==a and bo[1]+F(packet['g2_shift_exact'])==b
                    assert digest(root/base)==packet['base_witness_sha256']
                    assert len(parents)==2
                    for parent in parents:seed_needed.add(n(parent['certificate']))
                    tp=n(recipe['transferred_tail_record']);seed_needed.add(tp);tail=read(root/tp)
                    assert digest(root/tp)==recipe['transferred_tail_sha256']
                    assert tail['status']=='PASS_COMPLETE_TRANSFERRED_CONTINUUM_TAIL'
                    assert tail['witness_sha256']==packet['witness_sha256']and objective(tail)==[a,b]
                    assert n(tail['base_completed_certificate'])==n(parents[1]['certificate'])
                    baseproof=read(root/n(tail['base_completed_certificate']))
                    assert digest(root/n(tail['base_completed_certificate']))==tail['base_completed_certificate_sha256']
                    delta=F(packet['g2_shift_exact']);assert delta<0 and F(tail['g2_shift_exact'])==delta
                    finite=read(root/'output/boundary_length_20261007/real/angle268/finite_banded_cost.json')
                    infinite=read(root/'output/boundary_length_20261007/real/angle268/tail_split_replay.json')
                    assert len(tail['finite_tail_cells'])==len(finite['cells'])
                    intervals=[];tail_sum=F()
                    for tc,bc in zip(tail['finite_tail_cells'],finite['cells']):
                        assert (tc['lo'],tc['hi'],tc['region'])==(bc['lo'],bc['hi'],bc['region'])
                        assert F(tc['A_upper_exact'])==F(bc['A_upper_exact'])
                        assert F(tc['base_kappa_exact'])==F(bc['kappa_exact'])
                        k=F(tc['transferred_kappa_exact']);w=F(tc['reciprocal_weight_upper_exact'])
                        assert k==F(bc['kappa_exact'])-delta>0 and w>=1/(2*k)
                        cost=F(tc['A_upper_exact'])*w;assert cost==F(tc['cost_upper_exact']);tail_sum+=cost
                        intervals.append((F(tc['lo']),F(tc['hi'])))
                    intervals.sort();assert intervals[0][0]==0 and intervals[-1][1]==1
                    assert all(lo<hi for lo,hi in intervals)and all(x[1]==y[0]for x,y in zip(intervals,intervals[1:]))
                    kc=F(baseproof['coercivity_constants_exact']['compact']);ki=F(baseproof['coercivity_constants_exact']['infinite_UV'])
                    assert [x['region']for x in tail['infinite_tail_regions']]==['compact','remainder']
                    for tc in tail['infinite_tail_regions']:
                        region=tc['region'];k0=kc if region=='compact'else min(kc,ki)
                        assert F(tc['A_upper_exact'])==F(infinite['A_upper_exact'][region])and F(tc['base_kappa_exact'])==k0
                        k=F(tc['transferred_kappa_exact']);w=F(tc['reciprocal_weight_upper_exact'])
                        assert k==k0-delta>0 and w>=1/(2*k)
                        cost=F(tc['A_upper_exact'])*w;assert cost==F(tc['cost_upper_exact']);tail_sum+=cost
                    assert tail_sum==F(tail['tail_upper_exact'])==F(row['tail_upper_exact'])
                    assert tail_sum==sum(map(F,tail['components_exact'].values()),F())
            else:
                definition=read(root/wp);assert recipe['definition']==definition
                assert objective(definition)==[a,b]
                names=definition['names'];vals=list(map(F,definition['multipliers_exact']))
                if 'base_witness_paths'in definition:
                    weights=list(map(F,definition['weights_exact']));assert min(weights)>=0
                    bases=[]
                    for p in definition['base_witness_paths']:
                        p=n(p);seed_needed.add(p);bases.append(npz(p))
                    assert all(x[0]==names for x in bases)
                    assert vals==[sum(w*base[1][j]for w,base in zip(weights,bases))for j in range(len(names))]
                    assert [sum(w*base[3][j]for w,base in zip(weights,bases))for j in range(2)]==[a,b]
                    assert len(parents)==len(weights)
                    convex=F()
                    for j,(p,w)in enumerate(zip(definition['base_certificate_paths'],weights)):
                        p=n(p);seed_needed.add(p);pp=read(root/p)
                        assert parents[j]['certificate']==p and F(parents[j]['weight_exact'])==w
                        assert parents[j]['sha256']==digest(root/p)
                        assert objective(pp)==bases[j][3]
                        assert F(parents[j]['certified_constant_exact'])==support(pp)
                        convex+=w*support(pp)
                    assert convex==F(row['old_convex_bound_exact'])
                    removed=sum(F(x['old_convex_upper_exact'])for x in row['refined_components'])
                    assert removed==F(row['removed_components_exact'])
                    assert convex-removed==F(row['retained_other_spins_exact'])
                    assert convex-removed+sum(F(x['retained_upper_exact'])for x in row['refined_components'])==Ce
                else:
                    base=n(definition['base_witness']);seed_needed.add(base)
                    bn,bv,bh,bo=npz(base);assert bn==names
                    delta=F(definition['crossing_increment_exact'])
                    assert vals==[v+(delta if name=='X4'else 0)for name,v in zip(bn,bv)]
        if not row.get('positive_combination'):assert sum(map(F,fields(row).values()),F())==Ce
        assert rec['objective_exact']==r['objective_exact']
        assert [x['name']for x in rec['rows']]==names
        assert list(map(lambda x:F(x['exact']),rec['rows']))==vals
        assert len(vals)==r['bank_row_count']==rec['bank_row_count']
        assert sum(v!=0 for v in vals)==r['nonzero_multiplier_count']
        for j,(x,v)in enumerate(zip(rec['rows'],vals)):
            assert F(int(x['numerator']),int(x['denominator']))==v
            assert x['decimal_display']==decimal(v)
            if x['binary64_hex']:assert F.from_float(float.fromhex(x['binary64_hex']))==v
            all_csv.append(dict(cut_id=r['cut_id'],row_index=str(j),**{k:x[k]for k in ['name','numerator','denominator','decimal_display','binary64_hex','origin']}))
    with (data/'multipliers.csv').open(newline='',encoding='utf-8')as f:assert list(csv.DictReader(f))==all_csv
    raw_csv=[]
    for saved in raw['witnesses']:
        npz(saved['path'])
        for j,x in enumerate(saved['rows']):
            v=F(x['exact']);assert F(int(x['numerator']),int(x['denominator']))==v and x['decimal_display']==decimal(v)
            raw_csv.append(dict(witness_path=saved['path'],row_index=str(j),**{k:x[k]for k in ['name','numerator','denominator','decimal_display','binary64_hex']}))
    with (data/'raw_witness_multipliers.csv').open(newline='',encoding='utf-8')as f:assert list(csv.DictReader(f))==raw_csv
    expected_bounds=[];expected_components=[]
    for r in reg['records']:
        expected_bounds.append(dict(cut_id=r['cut_id'],upper_panel_edge_active=str(r['upper_panel_edge_active']),
             lower_panel_displayed=str(r['lower_panel_displayed']),a_exact=r['objective_exact'][0],b_exact=r['objective_exact'][1],
             plotted_support_exact=r['plotted_support_exact'],certified_support_exact=r['exact_certified_support'],
             certified_support_decimal_display=decimal(F(r['exact_certified_support'])),bank_row_count=str(r['bank_row_count']),
             source_certificate=r['source_certificate']or'',source_result_index=''if r.get('source_result_index')is None else str(r['source_result_index'])))
        ledger=r['selected_proof_ledger']
        if ledger.get('positive_combination'):
            cc={f'parent_{j}_weight_times_bound':str(F(p['weight_exact'])*F(p['certified_constant_exact']))for j,p in enumerate(r['parents'])}
        else:cc=ledger.get('effective_components_exact',{})
        if cc:assert sum(map(F,cc.values()),F())==F(r['exact_certified_support'])
        for key,value in cc.items():
            v=F(value);expected_components.append(dict(cut_id=r['cut_id'],interpretation='selected exact additive contribution',
                       component=key,numerator=str(v.numerator),denominator=str(v.denominator),decimal_display=decimal(v)))
    with (data/'bounds.csv').open(newline='',encoding='utf-8')as f:assert list(csv.DictReader(f))==expected_bounds
    with (data/'support_components.csv').open(newline='',encoding='utf-8')as f:assert list(csv.DictReader(f))==expected_components
    assert seed_needed<=set(read(data/'proof_seeds.json')['paths'])
    assert reg['selected_count']==len(union)==39 and reg['numerical_certificate_count']==38 and reg['analytic_count']==1
    assert reg['upper_vertex_count']==len(vertices)==39
    assert reg['lower_distinct_certificate_count']==len(bottom)==22
    assert reg['lower_interior_marker_count']==len(fig['displayed'])==21
    assert reg['lower_visual_marker_count']==fig['visual_marker_count']==23
    return dict(status='PASS_READ_ONLY_EXACT_FIGURE_EXPORT_RECHECK',selected_cuts=len(union),numerical_cuts=38,analytic_cuts=1,
                exact_multiplier_rows=len(all_csv),source_npz_witnesses=len(cache),exact_polygon_vertices=len(vertices),
                finite_pair_intersections_reproduce_full_polygon=True,
                scope='Source hashes, exact selection/geometry, all effective and parent coefficients, mixture/objective recipes, exact selected constants and CSV/JSON agreement. No interval-kernel replay or formal proof.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--data',type=Path)
    a=ap.parse_args();root=a.root.resolve();data=a.data.resolve()if a.data else root/'output/FIGURE_REPRO_20261008/data'
    print(json.dumps(verify(root,data),indent=2))
