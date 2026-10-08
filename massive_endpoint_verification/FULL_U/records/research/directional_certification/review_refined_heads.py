"""Check replayed refined g0 head partitions, totals and frozen provenance."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dyad(v):return F(int(v[0]))*F(2)**int(v[1])
def main():
    manifest={};results=[]
    for l in (2,4,6):
        paths=[f'output/directional_certification/refined{suffix}_spin{l}.json' for suffix in ('','_replay')]
        values=[];parts=[]
        for path in paths:
            manifest[path]=sha(path);data=json.loads((ROOT/path).read_text())
            assert data['status']=='PASS' and data['spin']==l and len(data['directions'])==1
            assert data['source_sha256']==sha('research/directional_certification/head_family.py')
            assert data['witness_sha256']==sha('output/angles_resummed/angle_000p000.npz')
            for p,h in data['dependency_sha256'].items():assert sha(p)==h;manifest[p]=h
            row=data['directions'][0];assert F(row['slope_exact'])==F(row['crossing_increment_exact'])==0
            cells=row['cells']+row['gaps'];ranges=sorted((F(c['a']),F(c['b'])) for c in cells)
            assert ranges[0][0]==0 and ranges[-1][1]==1
            assert all(a<b for a,b in ranges) and all(b==c for (a,b),(c,d) in zip(ranges,ranges[1:]))
            assert all(F(c['upper_exact'])>=0 for c in cells)
            value=sum((F(c['upper_exact']) for c in cells),F(0));assert value==F(row['head_upper_exact'])
            assert dyad(row['threshold']['new_V_over_p_upper_dyadic'])<0
            values.append(value);parts.append(ranges)
        assert parts[0]==parts[1]
        assert abs(values[0]-values[1])<F(1,10**12)
        cap=(max(values)*10**12).__ceil__()*F(1,10**12)+F(1,10**12)
        results.append(dict(spin=l,cells_with_gaps=len(parts[0]),producer_upper_exact=str(values[0]),replay_upper_exact=str(values[1]),difference_exact=str(values[1]-values[0]),common_decimal_cap_exact=str(cap),common_decimal_cap_float=float(cap)))
    here=str(Path(__file__).relative_to(ROOT));manifest[here]=sha(here)
    manifest['research/directional_certification/head_family.py']=sha('research/directional_certification/head_family.py')
    out=dict(status='PASS',replay_precision_bits=320,claim='Three full fresh refined head replays; exact total/partition/hash checks.',
        negative_slope_lemma='For b<=0,d=0 and p>0, H(U,V+b*p)<=H(U,V); hence each refined g0 spin bound is also valid for that spin at every negative b.',
        results=results,sha256_manifest=manifest)
    (ROOT/'output/directional_certification/refined_heads_independent_audit.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
