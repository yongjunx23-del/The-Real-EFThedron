"""Exact coverage/provenance assembly for new fixed150 directional bounds."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,pickle,sys,time
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/full_direction_20261007/certification'
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'research/directional_certification'))
import certify_directional_tail as tail_kernel

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(p):return F(int(p[0]))*F(2)**int(p[1])
def outward(q,n=6):
    k=-((-q.numerator*10**n)//q.denominator);return f'{k//10**n}.{k%10**n:0{n}d}'
def main():
    tic=time.time();hashes={}
    def read(path):
        data=json.loads(path.read_text());hashes[str(path.relative_to(ROOT))]=sha(path)
        for key in ('sha256_manifest','dependency_sha256'):
            for p,h in data.get(key,{}).items():assert sha(ROOT/p)==h,p;hashes[p]=h
        return data
    old=read(ROOT/'output/directional_certification/complete_directional_certificate.json');assert old['status']=='PASS'
    witness=old['witness_sha256'];tail=read(OUT/'tail_extended_18_5.json');finite=read(OUT/'finite_extended_18_5.json')
    source=pickle.loads((ROOT/'output/certification_tail/source_cache_v2_depth10.pkl').read_bytes())['data']
    cd=pickle.loads((ROOT/'output/certification_tail/compact_source_cache.pkl').read_bytes())
    tail_kernel.ctx.prec=340;tail_kernel.CASES=[(F(3,10),F(27,25))]
    uv,_=tail_kernel.uv_family(source);modified=dict(source);modified['lx']+=F(27,25)
    compact=tail_kernel.compact_run(42,modified,cd,'.02')
    fresh=max(F(uv[F(3,10)]['upper_exact']),tail_kernel.dyad(tail_kernel.arb(compact['low_energy_upper'])+tail_kernel.ar(F(3,10))),tail_kernel.dyad(tail_kernel.arb(compact['middle_energy_upper'])+tail_kernel.ar(F(3,10))))
    assert fresh<0
    records=[];head_records=[]
    for b in map(F,['-2','-4','-8','-16','3/10']):
        delta=0 if b<0 else F(18,5)*b;prefix='negative' if b<0 else 'positive';hs={}
        for l in [0,2,4,6]:
            values=[];partitions=[]
            for path in [OUT/f'{prefix}_head_spin{l}.json',OUT/f'replay_extended_head_spin{l}.json']:
                data=read(path);assert data['status']=='PASS' and data['witness_sha256']==witness and data['spin']==l
                assert data['source_sha256']==sha(ROOT/'research/full_direction_20261007/head_family_extended.py')
                row=next(r for r in data['directions'] if F(r['slope_exact'])==b)
                assert F(row['crossing_increment_exact'])==delta
                cells=row['cells']+row['gaps'];ranges=sorted((F(c['a']),F(c['b'])) for c in cells)
                assert ranges[0][0]==0 and ranges[-1][1]==1 and all(a<bb for a,bb in ranges)
                assert all(bb==c for (a,bb),(c,d) in zip(ranges,ranges[1:]))
                assert all(F(c['upper_exact'])>=0 for c in cells)
                total=sum((F(c['upper_exact']) for c in cells),F(0));assert total==F(row['head_upper_exact'])
                if l:assert dyad(row['threshold']['new_V_over_p_upper_dyadic'])<0
                values.append(total);partitions.append(ranges)
            assert partitions[0]==partitions[1];hs[str(l)]=max(values)
            head_records.append(dict(b=str(b),spin=l,producer_upper_exact=str(values[0]),replay_upper_exact=str(values[1]),difference_exact=str(values[1]-values[0]),cells_with_gaps=len(partitions[0])))
        t=next(r for r in tail['results'] if F(r['b'])==b);f=next(r for r in finite['results'] if F(r['b'])==b)
        assert t['status']==f['status']=='PASS' and F(t['delta_X4'])==F(f['delta_X4'])==delta
        assert tail['witness_sha256']==finite['witness_sha256']==witness
        k=F(t['tail_kappa_exact'])
        if b>0:k=min(k,-fresh)
        assert k>0;tail_bound=F(t['A_tail_upper_exact'])/(2*k);finite_bound=F(f['finite_support_upper_exact'])
        total=sum(hs.values(),F(0))+finite_bound+tail_bound
        records.append(dict(slope_exact=str(b),crossing_increment_exact=str(delta),certified_upper_exact=str(total),upper_outward_6dp=outward(total),
            inequality=f'g0 + ({b}) g2 <= {outward(total)}',heads_upper_exact={l:str(q) for l,q in hs.items()},finite8_40_upper_exact=str(finite_bound),tail42_infinity_upper_exact=str(tail_bound),
            normalized_support_plot_float=float(total)/math.sqrt(1+float(b)**2),angle_label_degrees=math.degrees(math.atan(float(b))),certification_kind='Direct fixed-witness continuum support upper bound'))
    for p in [Path(__file__),ROOT/'research/full_direction_20261007/head_family_extended.py']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(status='PASS',claim='Five additional continuum inequalities under the manuscript amplitude assumptions; all mu>4, all even spins.',
        witness_sha256=witness,results=records,head_replay=dict(producer_bits=280,replay_bits=340,scope='Fresh executions of the same reviewed interval implementation, not an independent implementation.',records=head_records),
        positive_tail_replay=dict(precision_bits=340,b='3/10',delta_X4='27/25',fresh_V_over_p_upper_exact=str(fresh)),
        method='Exact sum of replayed heads0/2/4/6, per-spin finite8..40 coercivity merge, and all-spin42+ analytic envelope.',
        seconds=time.time()-tic,sha256_manifest=hashes)
    (OUT/'complete_extended_certificate.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(status='PASS',bounds=[r['inequality'] for r in records],seconds=result['seconds']),indent=2))
if __name__=='__main__':main()
