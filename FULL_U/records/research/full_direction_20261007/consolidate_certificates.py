"""Consolidate exact complete bounds after further finite-block refinement."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,time
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/full_direction_20261007/certification'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def outward(q,n=6):
    k=-((-q.numerator*10**n)//q.denominator);return f'{k//10**n}.{k%10**n:0{n}d}'
def main():
    tic=time.time();hashes={}
    def read(path):
        d=json.loads(path.read_text());assert d['status']=='PASS';hashes[str(path.relative_to(ROOT))]=sha(path)
        for p,h in d.get('sha256_manifest',{}).items():assert sha(ROOT/p)==h,p;hashes[p]=h
        return d
    inputs=[read(ROOT/'output/directional_certification/complete_directional_certificate.json'),read(OUT/'complete_extended_certificate.json'),read(OUT/'complete_positive_wide_certificate.json')]
    witness=inputs[0]['witness_sha256'];assert all(d['witness_sha256']==witness for d in inputs)
    finite=read(OUT/'dense_finite_family.json');assert finite['witness_sha256']==witness
    records=[]
    for data in inputs:
        for r in data['results']:
            b,d=F(r['slope_exact']),F(r['crossing_increment_exact'])
            f=next(q for q in finite['results'] if F(q['b'])==b and F(q['delta_X4'])==d);assert f['status']=='PASS'
            original_finite=F(r['finite8_40_upper_exact']);new_finite=min(original_finite,F(f['finite_support_upper_exact']))
            hs=sum(map(F,r['heads_upper_exact'].values()),F(0));tail=F(r['tail42_infinity_upper_exact'])
            new_direct=hs+new_finite+tail;previous=F(r['certified_upper_exact']);best=min(previous,new_direct)
            records.append(dict(slope_exact=str(b),crossing_increment_exact=str(d),previous_upper_exact=str(previous),direct_component_upper_exact=str(new_direct),
                original_finite8_40_upper_exact=str(original_finite),finite8_40_upper_exact=str(new_finite),heads_upper_exact=r['heads_upper_exact'],tail42_infinity_upper_exact=str(tail),
                certified_upper_exact=str(best),improved_by_refined_block=new_direct<previous,positive_combination=None))
    assert len(records)==15 and len(set(r['slope_exact'] for r in records))==15
    before=[(F(r['slope_exact']),F(r['certified_upper_exact'])) for r in records]
    increments={F(r['slope_exact']):F(r['crossing_increment_exact']) for r in records}
    for r in records:
        b=F(r['slope_exact']);best=F(r['certified_upper_exact']);mix=None
        for a,ca in before:
            for z,cz in before:
                if a<b<z:
                    wa=(z-b)/(z-a);candidate=wa*ca+(1-wa)*cz
                    if candidate<best:
                        best=candidate;mix=dict(left_slope_exact=str(a),right_slope_exact=str(z),left_weight_exact=str(wa),right_weight_exact=str(1-wa),
                            left_certified_constant_exact=str(ca),right_certified_constant_exact=str(cz))
        effective=increments[b] if mix is None else F(mix['left_weight_exact'])*increments[F(mix['left_slope_exact'])]+F(mix['right_weight_exact'])*increments[F(mix['right_slope_exact'])]
        r.update(certified_upper_exact=str(best),positive_combination=mix,direct_crossing_increment_exact=str(increments[b]),certified_effective_crossing_increment_exact=str(effective),upper_outward_6dp=outward(best),inequality=f'g0 + ({b}) g2 <= {outward(best)}',
            angle_label_degrees=math.degrees(math.atan(float(b))),normalized_support_plot_float=float(best)/math.sqrt(1+float(b)**2))
    records.sort(key=lambda r:F(r['slope_exact']))
    hashes[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    result=dict(status='PASS',claim='Fifteen complete continuum inequalities, with sharpened finite blocks and exact positive-combination improvements.',witness_sha256=witness,
        results=records,seconds=time.time()-tic,scope='All mu>4 and all even spins, under manuscript amplitude assumptions; strong coverage only in these right-facing directions and their positive cone.',
        interpolation='Every nonnegative combination of the exact stored inequalities is valid. Normals and constants are rational; angle and normalized-support fields are plot labels only. crossing_increment_exact and the component fields describe the direct family; certified_effective_crossing_increment_exact describes the possibly interpolated bound. All other150-row coefficients remain those of the fixed base witness.',
        left_and_upward_warning='Use whole_plane_baseline.json only as the explicitly weak derived Gram baseline; no competitive all-circle certificate is asserted.',sha256_manifest=hashes)
    (OUT/'complete_fifteen_direction_certificate.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(status='PASS',bounds=[r['inequality'] for r in records],seconds=result['seconds']),indent=2))
if __name__=='__main__':main()
