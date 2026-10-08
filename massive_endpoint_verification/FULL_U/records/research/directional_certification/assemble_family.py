"""Assemble continuum inequalities from compatible, hash-bound components."""
import hashlib,json,math,sys
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent))
from head_family import outward,dyad
OUT=ROOT/'output/directional_certification'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    hashes={}
    def read(name):
        path=OUT/name;d=json.loads(path.read_text());hashes[str(path.relative_to(ROOT))]=sha(path)
        for key in ('sha256_manifest','dependency_sha256'):
            for p,h in d.get(key,{}).items():assert sha(ROOT/p)==h,(name,p);hashes[p]=h
        return d
    oldpath=ROOT/'output/certification_head/complete_upper_certificate.json'
    old=json.loads(oldpath.read_text());witness=old['witness_sha256']
    hashes[str(oldpath.relative_to(ROOT))]=sha(oldpath)
    heads={l:(read(f'head_spin{l}.json'),read(f'replay_head_spin{l}.json')) for l in (0,2,4,6)}
    refined={l:(read(f'refined_spin{l}.json'),read(f'refined_replay_spin{l}.json')) for l in (0,2,4,6)}
    for d in [x for pair in list(heads.values())+list(refined.values()) for x in pair]:
        assert d['status']=='PASS' and d['witness_sha256']==witness and d['partition_exactly_covers_beta_0_1']
        assert d['source_sha256']==sha(Path(__file__).parent/'head_family.py')
        for r in d['directions']:
            cells=r['cells']+r['gaps'];parts=sorted((F(c['a']),F(c['b'])) for c in cells)
            assert parts[0][0]==0 and parts[-1][1]==1 and all(a<b for a,b in parts)
            assert all(b==c for (a,b),(c,d) in zip(parts,parts[1:]))
            assert all(F(c['upper_exact'])>=0 for c in cells)
            assert sum((F(c['upper_exact']) for c in cells),F(0))==F(r['head_upper_exact'])
            if d['spin']:
                assert dyad(r['threshold']['new_V_over_p_upper_dyadic'])<0
    tail=read('tail_family_4b.json');finite=read('finite_block_family.json');replay=read('tail_replay.json')
    for name in ('head_family_independent_audit.json','refined_heads_independent_audit.json','finite_vectors_independent_audit.json'):
        audit=read(name);assert audit['status']=='PASS'
    assert tail['L']==42 and finite['status']=='PASS' and finite['spins']==list(range(8,42,2))
    assert tail['witness_sha256']==finite['witness_sha256']==witness
    assert replay['status']=='PASS' and replay['input_sha256']==sha(OUT/'tail_family_4b.json')
    assert replay['source_sha256']==sha(Path(__file__).parent/'review_tail_family.py')
    records=[];failed=[]
    def row(data,b,d,field='directions'):
        matches=[r for r in data[field] if F(r['slope_exact'] if field=='directions' else r['b'])==b and F(r['crossing_increment_exact'] if field=='directions' else r['delta_X4'])==d]
        assert len(matches)==1
        return matches[0]
    for t in tail['results']:
        b,d=F(t['b']),F(t['delta_X4']);assert d==(4*b if b>0 else 0)
        if t['status']!='PASS':failed.append({'slope_exact':str(b),'reason':t['status']});continue
        f=row(finite,b,d,'results');assert f['status']=='PASS'
        hs={}
        for l,pair in heads.items():
            vals=[F(row(data,b,d)['head_upper_exact']) for data in pair]
            # Retain the larger total from both independently invoked precisions.
            value=max(vals)
            if b<=0 and l in refined:
                # For d=0,b<=0, positivity of p implies H(U,V+b*p)<=H(U,V).
                # The refined g0 head is therefore also an upper bound here.
                zero=max(F(row(data,F(0),F(0))['head_upper_exact']) for data in refined[l])
                value=min(value,zero)
            hs[str(l)]=value
        k=F(t['tail_kappa_exact'])
        if b>0:
            r=next(r for r in replay['results'] if F(r['b'])==b and F(r['delta_X4'])==d)
            k=min(k,F(r['kappa_exact']))
        assert k>0
        ta=F(t['A_tail_upper_exact'])/(2*k);fi=F(f['finite_support_upper_exact'])
        total=sum(hs.values(),F(0))+fi+ta
        records.append(dict(slope_exact=str(b),crossing_increment_exact=str(d),angle_degrees=math.degrees(math.atan(float(b))),
            inequality=f'g0 + ({b}) g2 <= {outward(total,6)}',certified_upper_exact=str(total),
            upper_outward_6dp=outward(total,6),upper_outward_12dp=outward(total),
            heads_upper_exact={l:str(q) for l,q in hs.items()},finite8_40_upper_exact=str(fi),tail42_infinity_upper_exact=str(ta),
            normalized_support_upper_float=float(total)/math.sqrt(1+float(b)**2),
            head_method='Maximum of256/320bit runs; for b<=0 also use the replayed refined g0 low-spin caps by monotonicity.'))
    assert len(records)==7
    # The strongest positive combination of two certified slopes is again a
    # valid inequality. Exact pairwise interpolation also removes loose points
    # from the lower convex envelope of the certified constants.
    direct=[(F(r['slope_exact']),F(r['certified_upper_exact'])) for r in records]
    for r in records:
        b=F(r['slope_exact']);raw=F(r['certified_upper_exact']);best=raw;mix=None
        for left,cl in direct:
            for right,cr in direct:
                if left<b<right:
                    wl=(right-b)/(right-left);value=wl*cl+(1-wl)*cr
                    if value<best:best=value;mix=dict(left_slope=str(left),right_slope=str(right),left_weight_exact=str(wl),right_weight_exact=str(1-wl))
        r['direct_witness_upper_exact']=str(raw)
        r['positive_interpolation_improvement']=mix
        r['certified_upper_exact']=str(best);r['upper_outward_6dp']=outward(best,6);r['upper_outward_12dp']=outward(best)
        r['inequality']=f'g0 + ({b}) g2 <= {outward(best,6)}'
        r['normalized_support_upper_float']=float(best)/math.sqrt(1+float(b)**2)
    here=Path(__file__);hashes[str(here.relative_to(ROOT))]=sha(here)
    for name in ('head_family.py','review_tail_family.py'):
        p=Path(__file__).parent/name;hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(status='PASS',claim='Seven continuum upper bounds g0+b*g2<=C, all even spins and allmu>4, under the manuscript amplitude assumptions.',
        old_certified_g0_upper=old['certified_g0_upper_outward_6dp'],witness_sha256=witness,results=records,uncertified_directions=failed,
        normalization='Exact rational normals(1,b). Angle labels and normalized plotting supports are floating visualizations; exact claims are the rational inequalities.',
        interpolation='Forb between adjacent certified slopes, linear interpolation of the corresponding exact upper constants is a certified bound. Outward6dp constants can also be interpolated safely.',
        sha256_manifest=hashes)
    (OUT/'complete_directional_certificate.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'status':'PASS','bounds':[r['inequality'] for r in records],'uncertified':failed},indent=2))
if __name__=='__main__':main()
