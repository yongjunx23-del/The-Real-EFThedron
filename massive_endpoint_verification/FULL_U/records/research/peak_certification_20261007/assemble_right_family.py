"""Eight full continuum bounds with the fixed original 150-row dual witness."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research'))
from angle_benchmarks import support_report
OUT=ROOT/'output/peak_certification_20261007/right_family'
SLOPES=list(map(F,['-20','-24','-28','-40','-48','-96','-192','-256']))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dy(p):return F(int(p[0]))*F(2)**int(p[1])
def outward(q,n=12):
    k=-((-q.numerator*10**n)//q.denominator);return f'{k//10**n}.{k%10**n:0{n}d}'
def main():
    assert not (OUT/'complete_right_family.json').exists()
    manifest={};seen=set()
    def bind(p,h=None):
        p=Path(p);p=p if p.is_absolute()else ROOT/p
        key=str(p.relative_to(ROOT));hh=sha(p);assert h is None or hh==h,key
        assert key not in manifest or manifest[key]==hh,key
        manifest[key]=hh;return p
    def read(p):
        p=bind(p);d=json.loads(p.read_text())
        if p not in seen:
            seen.add(p)
            for name in ['sha256_manifest','dependency_sha256']:
                for q,h in d.get(name,{}).items():bind(q,h)
        return d
    base=read('output/full_direction_20261007/certification/complete_fifteen_direction_certificate.json')
    tail=read('output/certification_tail/fixed_tail_certificate.json')
    assert base['status']=='PASS'
    assert tail['status']=='All-spin/all-energy omitted-spin tail certified under the stated amplitude assumptions.'
    assert tail['replayed_components_this_run']and F(tail['V_over_p_upper_rational'])==-F(31,50)
    witness=bind('output/angles_resummed/angle_000p000.npz',base['witness_sha256'])
    assert tail['witness_sha256']==sha(witness)
    st=np.load(witness);assert np.array_equal(st['direction'],np.array([1.,0.]))and len(st['names'])==150
    finite=[read(OUT/f)for f in ['finite_family.json','finite_family_replay.json']]
    for d,bits in zip(finite,[280,340]):
        assert d['status']=='PASS'and d['witness_sha256']==sha(witness)and d['precision_bits']==bits
        assert d['spins']==list(range(8,42,2))and [F(r['b'])for r in d['results']]==SLOPES
    head={}
    for ell in [0,2,4,6]:
        pair=[read(OUT/f'head_p{bits}_spin{ell}.json')for bits in [280,340]]
        for d,bits in zip(pair,[280,340]):
            assert d['status']=='PASS'and d['spin']==ell and d['precision_bits']==bits and d['witness_sha256']==sha(witness)
            assert d['source_sha256']==sha(ROOT/'research/full_direction_20261007/head_family_extended.py')
            assert [F(r['slope_exact'])for r in d['directions']]==SLOPES
        head[ell]=pair
    records=[]
    for b in SLOPES:
        components={};checks=[]
        for ell in [0,2,4,6]:
            values=[];partitions=[]
            for d in head[ell]:
                row=next(r for r in d['directions']if F(r['slope_exact'])==b)
                assert F(row['crossing_increment_exact'])==0
                cells=row['cells']+row['gaps'];part=sorted((F(c['a']),F(c['b']))for c in cells)
                assert part[0][0]==0 and part[-1][1]==1 and all(a<z for a,z in part)and all(z==c for(a,z),(c,e)in zip(part,part[1:]))
                assert all(F(c['upper_exact'])>=0 for c in cells)
                value=sum((F(c['upper_exact'])for c in cells),F(0));assert value==F(row['head_upper_exact'])
                if ell:assert dy(row['threshold']['new_V_over_p_upper_dyadic'])<0
                values.append(value);partitions.append(part)
            assert partitions[0]==partitions[1]
            components[f'head_spin{ell}']=str(max(values));checks.append({'spin':ell,'producer_exact':str(values[0]),'replay_exact':str(values[1]),'cells_and_gaps':len(partitions[0])})
        vals=[]
        for d in finite:
            row=next(r for r in d['results']if F(r['b'])==b)
            assert row['status']=='PASS'and F(row['delta_X4'])==0
            assert sorted(map(int,row['support_upper_by_spin_dyadic']))==list(range(8,42,2))
            assert F(row['finite_support_upper_exact'])>=sum((dy(p)for p in row['support_upper_by_spin_dyadic'].values()),F(0))
            vals.append(F(row['finite_support_upper_exact']))
        components['finite_spins8_40']=str(max(vals))
        k=F(31,50)-b;assert k>0
        components['tail_spins42_infinity']=str(F(tail['A40_upper_exact'])/(2*k))
        total=sum(map(F,components.values()),F(0));norm2=1+b*b;scale=2**256
        normlo=F(math.isqrt(norm2.numerator*scale*scale//norm2.denominator),scale)
        assert normlo**2<=norm2<(normlo+F(1,scale))**2
        bound=outward(total);normalized=outward(F(bound)/normlo)
        angle=math.degrees(math.atan2(float(b),1))%360
        definition={'format':'exact_objective_shift_of_frozen_witness','base_witness':str(witness.relative_to(ROOT)),'base_witness_sha256':sha(witness),'objective_exact':['1',str(b)],'names':st['names'].tolist(),'multipliers_exact':[str(F(float(v)))for v in st['multipliers']],'crossing_increment_exact':'0'}
        dp=OUT/f'witness_b{abs(b.numerator)}.json';assert not dp.exists();dp.write_text(json.dumps(definition,indent=2));bind(dp)
        rec={'status':'PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND','witness':str(dp.relative_to(ROOT)),'witness_format':definition['format'],'witness_sha256':sha(dp),'objective_exact':['1',str(b)],'nominal_angle_degrees':angle,'support_upper_exact':str(total),'support_upper_outward_12dp':bound,'support_normalized_upper_outward_12dp':normalized,'support_normalized_float':float(normalized),'components_exact':components,'head_checks':checks,'finite_replay_upper_exact':list(map(str,vals)),'tail_kappa_exact':str(k),'coverage':'All even spins and all mu>4. Spins0,2,4,6 have fresh full-energy interval integration; spins8..40 use exact overlaps of authenticated real mass and spectral sign partitions; every even spin>=42 uses the frozen full-domain real norm and the strictly strengthened sign bound.','proof_scope':'Directed interval arithmetic with higher-precision same-implementation replays and exact rational accounting. Not formal verification.','graphical_comparison':{'dual':support_report(angle),'primal':support_report(angle,'primal')},'sha256_manifest':{}}
        records.append(rec)
    bind(Path(__file__));bind(ROOT/'research/full_direction_20261007/head_family_extended.py')
    for rec in records:
        rec['sha256_manifest']=manifest.copy();tag=abs(int(F(rec['objective_exact'][1])))
        p=OUT/f'complete_b{tag}.json';assert not p.exists();p.write_text(json.dumps(rec,indent=2))
        print(tag,rec['support_normalized_float'],rec['graphical_comparison']['dual']['support'],flush=True)
    (OUT/'complete_right_family.json').write_text(json.dumps({'status':'PASS_EIGHT_FULL_CONTINUUM_DIRECTIONAL_BOUNDS','certificates':[str((OUT/f'complete_b{abs(b.numerator)}.json').relative_to(ROOT))for b in SLOPES],'results':records,'sha256_manifest':manifest},indent=2))
if __name__=='__main__':main()
