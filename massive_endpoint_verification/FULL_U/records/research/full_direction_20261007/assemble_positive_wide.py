"""Hash-bound complete bounds using the integrable positive-UV extension."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,sys,time
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/full_direction_20261007/certification'
sys.dont_write_bytecode=True
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(p):return F(int(p[0]))*F(2)**int(p[1])
def outward(q,n=6):
    k=-((-q.numerator*10**n)//q.denominator);return f'{k//10**n}.{k%10**n:0{n}d}'
def main():
    tic=time.time();hashes={}
    def read(path):
        d=json.loads(path.read_text());hashes[str(path.relative_to(ROOT))]=sha(path)
        for key in ('sha256_manifest','dependency_sha256'):
            for p,h in d.get(key,{}).items():assert sha(ROOT/p)==h,p;hashes[p]=h
        return d
    old=read(ROOT/'output/directional_certification/complete_directional_certificate.json');assert old['status']=='PASS'
    witness=old['witness_sha256'];finite=read(OUT/'finite_positive_refined_family.json');assert finite['status']=='PASS'
    tails=[read(OUT/p) for p in ['positive_uv_tail.json','positive_uv_tail_replay.json']]
    fresh_vectors=read(OUT/'finite_v_positive_refined.json');assert fresh_vectors['status']=='PASS'
    assert fresh_vectors['witness_sha256']==finite['witness_sha256']==witness
    for c in fresh_vectors['partition']:
        assert set(map(int,c['upper_by_spin']))==set(range(8,42,2))
        for b in map(F,['1/2','1','2']):
            assert all(dyad(q['V_upper'])+b+F(18,5)*b*dyad(q['X4_upper'])<0 for q in c['upper_by_spin'].values())
    # The exact finite merge already checks gapless domain coverage and every
    # overlap. Bind that producer and both underlying vector tables as inputs.
    results=[];head_checks=[];tail_checks=[]
    for b in map(F,['1/2','1','2']):
        delta=F(18,5)*b;heads={}
        for l in [0,2,4,6]:
            first=f'half_head_spin{l}.json' if b==F(1,2) else f'wide_head_spin{l}.json'
            values=[];parts=[]
            for path in [OUT/first,OUT/f'wide_replay_head_spin{l}.json']:
                data=read(path);assert data['status']=='PASS' and data['witness_sha256']==witness and data['spin']==l
                assert data['source_sha256']==sha(ROOT/'research/full_direction_20261007/head_family_extended.py')
                r=next(r for r in data['directions'] if F(r['slope_exact'])==b);assert F(r['crossing_increment_exact'])==delta
                cells=r['cells']+r['gaps'];part=sorted((F(c['a']),F(c['b'])) for c in cells)
                assert part[0][0]==0 and part[-1][1]==1 and all(a<bb for a,bb in part)
                assert all(bb==c for (a,bb),(c,d) in zip(part,part[1:]))
                assert all(F(c['upper_exact'])>=0 for c in cells)
                value=sum((F(c['upper_exact']) for c in cells),F(0));assert value==F(r['head_upper_exact'])
                if l:assert dyad(r['threshold']['new_V_over_p_upper_dyadic'])<0
                values.append(value);parts.append(part)
            assert parts[0]==parts[1];heads[str(l)]=max(values)
            head_checks.append(dict(b=str(b),spin=l,producer_upper_exact=str(values[0]),replay_upper_exact=str(values[1]),difference_exact=str(values[1]-values[0])))
        tailvals=[]
        for data in tails:
            assert data['witness_sha256']==witness and data['domain']=='All even ell>=42 and all mu>4'
            r=next(r for r in data['results'] if F(r['b'])==b);assert r['status']=='PASS' and F(r['delta_X4'])==delta
            ranges=[(F(c['tau_a']),F(c['tau_b'])) for c in data['uv_cells']]
            assert ranges[0][0]==0 and ranges[-1][1]==256 and all(a<bb for a,bb in ranges)
            assert all(bb==c for (a,bb),(c,d) in zip(ranges,ranges[1:]))
            eta=F(r['best']['eta_exact']);assert eta>0
            assert F(r['compact_V_over_p_upper_exact'])+eta<=0 and F(r['large_tau_V_over_p_upper_exact'])+eta<=0
            assert all(F(c['p_mass_upper_exact'])>=0 for c in data['uv_cells'])
            extra=sum((2*max(F(c['V_over_p_upper_by_slope'][str(b)])+eta,0)*F(c['p_mass_upper_exact']) for c in data['uv_cells']),F(0))
            value=F(data['all_spin_real_norm_upper_exact'])/(2*eta)+extra
            assert extra==F(r['best']['positive_V_term_exact']) and value==F(r['best']['tail_support_upper_exact'])
            tailvals.append(value)
        f=next(r for r in finite['results'] if F(r['b'])==b);assert f['status']=='PASS' and F(f['delta_X4'])==delta
        fi=F(f['finite_support_upper_exact']);ta=max(tailvals);total=sum(heads.values(),F(0))+fi+ta
        tail_checks.append(dict(b=str(b),producer_tail_exact=str(tailvals[0]),replay_tail_exact=str(tailvals[1])))
        results.append(dict(slope_exact=str(b),crossing_increment_exact=str(delta),certified_upper_exact=str(total),upper_outward_6dp=outward(total),
            inequality=f'g0 + ({b}) g2 <= {outward(total)}',heads_upper_exact={k:str(v) for k,v in heads.items()},finite8_40_upper_exact=str(fi),tail42_infinity_upper_exact=str(ta),
            normalized_support_plot_float=float(total)/math.sqrt(1+float(b)**2),angle_label_degrees=math.degrees(math.atan(float(b))),certification_kind='Direct fixed-witness continuum support upper bound with integrable positive-V UV remainder'))
    for p in [Path(__file__),ROOT/'research/full_direction_20261007/head_family_extended.py']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(status='PASS',claim='Three additional continuum inequalities at slopes1/2,1,2; all energies and all even spins.',witness_sha256=witness,
        results=results,head_replay=head_checks,tail_replay=tail_checks,seconds=time.time()-tic,
        review_scope='Fresh same-implementation interval replays: heads320/360bits and UV280/360bits. Elementary new tail inequality and spin-sum normalization reviewed by parent agent; no independent implementation or formal proof claim.',
        sha256_manifest=hashes)
    (OUT/'complete_positive_wide_certificate.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(status='PASS',bounds=[r['inequality'] for r in results],seconds=result['seconds']),indent=2))
if __name__=='__main__':main()
