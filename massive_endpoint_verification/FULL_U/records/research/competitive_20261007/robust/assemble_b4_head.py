"""Read-only exact assembly of the b4 finite0..40 head components.

Checks both precision runs and their exact partitions/provenance. This does
not freshly evaluate kernels, and it does not include the l>=42 tail.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/competitive_20261007/robust/b4'
EXPECTED='bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ceil_decimal(x,n=12):
    v=-((-x.numerator*10**n)//x.denominator)
    return f'{v//10**n}.{v%10**n:0{n}d}'
def main():
    hashes={};chosen={};comparisons={}
    for ell in [0,2,4,6]:
        results=[]
        for prec in [280,340]:
            p=OUT/f'head_p{prec}_spin{ell}.json';q=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
            assert q['status']=='PASS' and q['spin']==ell and q['precision_bits']==prec and q['witness_sha256']==EXPECTED
            assert q['repair_ratio_exact']=='18/5' and len(q['directions'])==1
            for f,h in q['dependency_sha256'].items():assert sha(ROOT/f)==h,f;hashes[f]=h
            source=ROOT/'research/full_direction_20261007/head_family_extended.py';assert sha(source)==q['source_sha256'];hashes[str(source.relative_to(ROOT))]=sha(source)
            d=q['directions'][0];assert d['slope_exact']=='4' and d['crossing_increment_exact']=='72/5'
            parts=d['cells']+d['gaps'];ranges=sorted((F(c['a']),F(c['b']))for c in parts)
            assert ranges[0][0]==0 and ranges[-1][1]==1 and all(a<b for a,b in ranges)
            assert all(b==c for(a,b),(c,e)in zip(ranges,ranges[1:]))
            values=[F(c['upper_exact'])for c in parts];assert all(v>=0 for v in values)
            total=sum(values,F(0));assert total==F(d['head_upper_exact'])
            if ell:
                dy=d['threshold']['new_V_over_p_upper_dyadic'];assert F(int(dy[0]))*F(2)**int(dy[1])<0
            results.append((total,ranges))
        assert results[0][1]==results[1][1]
        diff=abs(results[0][0]-results[1][0]);assert diff<F(1,10**12)
        chosen[str(ell)]=str(max(r[0]for r in results));comparisons[str(ell)]=str(diff)
    p=OUT/'finite_block.json';q=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
    assert q['status']=='PASS' and q['witness_sha256']==EXPECTED and q['spins']==list(range(8,42,2))
    for f,h in q['sha256_manifest'].items():assert sha(ROOT/f)==h,f;hashes[f]=h
    d=q['results'][0];assert len(q['results'])==1 and d['status']=='PASS' and d['b']=='4' and d['delta_X4']=='72/5'
    chosen['8..40']=d['finite_support_upper_exact']
    total=sum(map(F,chosen.values()),F(0));hashes[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    result=dict(status='PASS_FINITE_HEAD_0_40',witness_sha256=EXPECTED,direction_exact=['1','4'],crossing_increment_exact='72/5',
      finite_head_upper_exact=str(total),finite_head_upper_outward_12dp=ceil_decimal(total),component_upper_exact=chosen,
      paired_precision_absolute_difference_exact=comparisons,sha256_manifest=hashes,
      scope='Exact record, partition, parameter and dependency audit plus sum; kernels were evaluated by the two head producer runs. Infinite tail l>=42 is not included.')
    (OUT/'head_bundle.json').write_text(json.dumps(result,indent=2));print(result['status'],result['finite_head_upper_outward_12dp'])
if __name__=='__main__':main()
