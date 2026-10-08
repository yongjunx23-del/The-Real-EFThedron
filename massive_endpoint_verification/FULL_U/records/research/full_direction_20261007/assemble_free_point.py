"""Three complete negative-slope extensions, preserving the frozen15 ledger."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,time
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/full_direction_20261007/certification'
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
    frozen=read(OUT/'complete_fifteen_direction_certificate.json');assert frozen['status']=='PASS'
    tail=read(ROOT/'output/certification_tail/fixed_tail_certificate.json')
    assert F(tail['V_over_p_upper_rational'])==-F(31,50) and tail['replayed_components_this_run']
    finite=read(OUT/'free_point_finite_family.json');assert finite['status']=='PASS'
    witness=frozen['witness_sha256'];assert finite['witness_sha256']==tail['witness_sha256']==witness
    previous=F(next(r['certified_upper_exact'] for r in frozen['results'] if F(r['slope_exact'])==-16))
    records=[];head_checks=[]
    for b in map(F,['-32','-64','-128']):
        heads={}
        for l in [0,2,4,6]:
            vals=[];parts=[]
            for precision in [280,340]:
                data=read(OUT/f'free_point_head_p{precision}_spin{l}.json')
                assert data['status']=='PASS' and data['witness_sha256']==witness and data['spin']==l and data['precision_bits']==precision
                assert data['source_sha256']==sha(ROOT/'research/full_direction_20261007/head_family_extended.py')
                row=next(r for r in data['directions'] if F(r['slope_exact'])==b)
                assert F(row['crossing_increment_exact'])==0
                cells=row['cells']+row['gaps'];part=sorted((F(c['a']),F(c['b'])) for c in cells)
                assert part[0][0]==0 and part[-1][1]==1 and all(a<bb for a,bb in part)
                assert all(bb==c for (a,bb),(c,d) in zip(part,part[1:])) and all(F(c['upper_exact'])>=0 for c in cells)
                value=sum((F(c['upper_exact']) for c in cells),F(0));assert value==F(row['head_upper_exact'])
                if l:assert dyad(row['threshold']['new_V_over_p_upper_dyadic'])<0
                vals.append(value);parts.append(part)
            assert parts[0]==parts[1];heads[str(l)]=max(vals)
            head_checks.append(dict(b=str(b),spin=l,producer_upper_exact=str(vals[0]),replay_upper_exact=str(vals[1]),difference_exact=str(vals[1]-vals[0])))
        f=next(r for r in finite['results'] if F(r['b'])==b);assert f['status']=='PASS' and F(f['delta_X4'])==0
        fi=F(f['finite_support_upper_exact']);k=F(31,50)-b;ta=F(tail['A40_upper_exact'])/(2*k)
        direct=sum(heads.values(),F(0))+fi+ta;best=min(direct,previous)
        records.append(dict(slope_exact=str(b),crossing_increment_exact='0',direct_crossing_increment_exact='0',certified_effective_crossing_increment_exact='0',certified_upper_exact=str(best),direct_witness_upper_exact=str(direct),upper_outward_6dp=outward(best),
            inequality=f'g0 + ({b}) g2 <= {outward(best)}',heads_upper_exact={l:str(q) for l,q in heads.items()},finite8_40_upper_exact=str(fi),tail42_infinity_upper_exact=str(ta),tail_kappa_exact=str(k),
            previous_minus16_consequence_exact=str(previous),strictly_improves_previous_consequence=direct<previous,
            normalized_support_plot_float=float(best)/math.sqrt(1+float(b)**2),angle_label_degrees=360+math.degrees(math.atan(float(b)))))
    hashes[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    result=dict(status='PASS',claim='Three further complete continuum bounds near the free-point direction; all energies and all even spins.',witness_sha256=witness,results=records,
        tail_proof='U unchanged; V_new/p=V_old/p+b<=-(31/50-b). Use the old all-domain A40 and h<=U^2/(2 kappa p).',
        head_replay=dict(producer_bits=280,replay_bits=340,records=head_checks,scope='Fresh same-implementation replay, not independent code.'),
        seconds=time.time()-tic,sha256_manifest=hashes)
    (OUT/'complete_free_point_certificate.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(status='PASS',bounds=[r['inequality'] for r in records],seconds=result['seconds']),indent=2))
if __name__=='__main__':main()
