"""Read-only exact assembly of paired b6 head0..40 components."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'output/competitive_20261007/upper_arc/b6'
EXPECTED='bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upward(q,n=12):
    k=-((-q.numerator*10**n)//q.denominator);return f'{k//10**n}.{k%10**n:0{n}d}'
def main():
    hashes={};components={};differences={}
    for ell in [0,2,4,6]:
        results=[]
        for prec in [280,340]:
            p=OUT/f'head_p{prec}_spin{ell}.json';q=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
            assert q['status']=='PASS'and q['spin']==ell and q['precision_bits']==prec and q['witness_sha256']==EXPECTED
            assert q['repair_ratio_exact']=='107/30'and len(q['directions'])==1
            for f,h in q['dependency_sha256'].items():assert sha(ROOT/f)==h,f;hashes[f]=h
            source=ROOT/'research/full_direction_20261007/head_family_extended.py';assert sha(source)==q['source_sha256'];hashes[str(source.relative_to(ROOT))]=sha(source)
            d=q['directions'][0];assert d['slope_exact']=='6'and d['crossing_increment_exact']=='107/5'
            parts=d['cells']+d['gaps'];ranges=sorted((F(c['a']),F(c['b']))for c in parts)
            assert ranges[0][0]==0 and ranges[-1][1]==1 and all(a<b for a,b in ranges)
            assert all(b==c for(a,b),(c,e)in zip(ranges,ranges[1:]))
            vals=[F(c['upper_exact'])for c in parts];assert all(v>=0 for v in vals)
            total=sum(vals,F(0));assert total==F(d['head_upper_exact'])
            if ell:
                m,e=d['threshold']['new_V_over_p_upper_dyadic'];assert F(int(m))*F(2)**int(e)<0
            results.append((total,ranges))
        assert results[0][1]==results[1][1]
        diff=abs(results[0][0]-results[1][0]);assert diff<F(1,10**12)
        components[str(ell)]=str(max(r[0]for r in results));differences[str(ell)]=str(diff)
    p=OUT/'finite_block.json';q=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
    assert q['status']=='PASS_FINITE_SUPPORT_8_40'and q['witness_sha256']==EXPECTED
    assert q['direction_exact']==['1','6']and q['crossing_increment_exact']=='107/5'and q['spins']==list(range(8,42,2))
    for f,h in q['sha256_manifest'].items():assert sha(ROOT/f)==h,f;hashes[f]=h
    components['8..40']=q['finite_support_upper_exact'];total=sum(map(F,components.values()),F(0))
    for p in [Path(__file__),Path(__file__).with_name('run_heads.py')]:hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(status='PASS_FINITE_HEAD_0_40',witness_sha256=EXPECTED,objective_exact=['1','6'],crossing_increment_exact='107/5',
      finite_head_upper_exact=str(total),finite_head_upper_outward_12dp=upward(total),component_upper_exact=components,
      paired_precision_absolute_difference_exact=differences,sha256_manifest=hashes,
      scope='All even0..40 and every mu>4. The record audit verifies exact partitions, sums, parameters, and hashes. Kernel integrals were freshly computed by paired280/340 low-head runs and paired320/380 finite-sign runs. No omitted spin tail is included.')
    (OUT/'head_bundle.json').write_text(json.dumps(result,indent=2));print(result['status'],result['finite_head_upper_outward_12dp'],flush=True)
if __name__=='__main__':main()
