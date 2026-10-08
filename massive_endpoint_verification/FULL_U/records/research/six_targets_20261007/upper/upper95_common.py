"""Exact83+b6 mixture definitions and conservative unions of unequal source gaps."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/six_targets_20261007/upper95'
CERTS=['output/peak_certification_20261007/complete_upper83_certificate.json','output/competitive_20261007/upper_arc/b6/complete_b6_certificate.json']
WITNESSES=['output/weak_regions_20261007/upper/angle83_repaired.npz','output/angles_resummed/angle_000p000.npz']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dy(p):return F(int(p[0]))*F(2)**int(p[1])
def partition(cells,gaps):
    z=sorted(cells+gaps);assert z[0][0]==0 and z[-1][1]==1
    assert all(a<b for a,b in z)and all(x[1]==y[0]for x,y in zip(z,z[1:]))
def source_paths(ell):
    return [f'output/peak_certification_20261007/upper83/head/{"spin2_refined"if ell==2 else "spin0"}/head_l{ell}.json',f'output/competitive_20261007/upper_arc/b6/head_p280_spin{ell}.json']
def canonical_records(raw,ell):
    a,b=raw;assert a['spin']==b['spin']==ell
    assert a['status']=='PASS_FINITE_SPIN_ENCLOSURE'and a['final']
    assert b['status']=='PASS'and b['partition_exactly_covers_beta_0_1']and F(b['repair_ratio_exact'])==F(107,30)
    assert len(b['directions'])==1;old=b['directions'][0]
    assert F(old['slope_exact'])==6 and F(old['crossing_increment_exact'])==F(107,5)
    cells=[[(F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in a['cells_data']],[(F(c['a']),F(c['b']),F(c['upper_exact']))for c in old['cells']]]
    gaps=[[tuple(map(F,g))for g in a['gap_beta_intervals']],[(F(g['a']),F(g['b']))for g in old['gaps']]]
    bounds=[dy(a['gap_upper_dyadic']),sum((F(g['upper_exact'])for g in old['gaps']),F(0))]
    for cc,gg in zip(cells,gaps):partition([(x,y)for x,y,v in cc],gg);assert all(v>=0 for x,y,v in cc)
    assert bounds[0]+sum((v for x,y,v in cells[0]),F(0))<=dy(a['upper_dyadic'])
    assert bounds[1]+sum((v for x,y,v in cells[1]),F(0))==F(old['head_upper_exact'])
    return cells,gaps,bounds
def union_gap_data(cells,gaps,bounds):
    merged=[]
    for a,b in sorted(gaps[0]+gaps[1]):
        if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else:merged.append((a,b))
    enlarged=[];ledgers=[]
    for cc,gg,bound in zip(cells,gaps,bounds):
        used=[j for j,(a,b,u)in enumerate(cc)if any(a<hi and lo<b for lo,hi in merged)]
        extra=sum((cc[j][2]for j in used),F(0));enlarged.append(bound+extra)
        ledgers.append(dict(original_gap_exact=str(bound),intersecting_cell_indices=used,added_whole_cell_bounds_exact=str(extra),enlarged_gap_exact=str(bound+extra)))
    knots=sorted({v for cc in cells for a,b,u in cc for v in[a,b]}|{v for a,b in merged for v in[a,b]})
    ordinary=[(a,b)for a,b in zip(knots,knots[1:])if not any(a>=lo and b<=hi for lo,hi in merged)]
    partition(ordinary,merged)
    return ordinary,merged,enlarged,ledgers
def exact_weights(direction,base):
    a,b=base;det=a[0]*b[1]-a[1]*b[0]
    x=(direction[0]*b[1]-direction[1]*b[0])/det;y=(a[0]*direction[1]-a[1]*direction[0])/det
    assert min(x,y)>0 and[x*a[j]+y*b[j]for j in[0,1]]==direction
    return[x,y]
