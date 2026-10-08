"""Refine only saved spectral cells failing the new affine sign tests."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/directional_certification')]
from certify_finite_vectors import FiniteV,arb,ctx,ar,tight,rev,C,EXPECTED,dyad,frac
OUT=ROOT/'output/full_direction_20261007/certification'
CASES=[(F(b),F(18,5)*F(b)) for b in ['1','2']]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def good(row):return all(frac(q['V_upper'])+b+d*frac(q['X4_upper'])<0 for q in row.values() for b,d in CASES)
def evaluate(bank,kind,a,b):
    x=arb(ar((a+b)/2),ar((b-a)/2));row={}
    values=bank.fields(x) if kind=='mu' else bank.uv_fields(x)
    for l,v in values.items():
        if kind=='mu':
            delta=x-4;rs=tight(bank.roy.polys[l][1],delta)
            cross=tight(bank.cross[l],delta)/((x-ar(C))*rs)
        else:
            rs=tight(rev(bank.roy.polys[l][1]),x)
            cross=x*x*tight(rev(bank.cross[l],l-1),x)/((1+2*ar(C)*x)*rs)
        if not(v.is_finite() and cross.is_finite()):return None
        row[str(l)]=dict(V_upper=dyad(v),X4_upper=dyad(cross))
    return row
def main():
    tic=time.time();ctx.prec=320;path=ROOT/'output/directional_certification/finite_v_vectors.json'
    old=json.loads(path.read_text());assert old['status']=='PASS' and old['witness_sha256']==EXPECTED
    hashes=dict(old['sha256_manifest']);hashes[str(path.relative_to(ROOT))]=sha(path)
    for p,h in hashes.items():assert sha(ROOT/p)==h,p
    records=[];stack=[]
    for c in old['partition']:
        if good(c['upper_by_spin']):records.append(c)
        else:stack.append((c['variable'],F(c['a']),F(c['b']),0))
    original_bad=len(stack);print(json.dumps(dict(reused=len(records),cells_to_refine=original_bad)),flush=True)
    bank=FiniteV();trials=0
    while stack:
        kind,a,b,depth=stack.pop();trials+=1;assert trials<10000
        row=evaluate(bank,kind,a,b)
        if row is not None and good(row):records.append(dict(variable=kind,a=str(a),b=str(b),upper_by_spin=row));continue
        assert depth<15,('Could not establish sign',kind,a,b,depth,row)
        mid=(a+b)/2
        if depth==0:
            atmid=evaluate(bank,kind,mid,mid)
            assert atmid is not None and good(atmid),('Actual midpoint sign obstruction',kind,mid,atmid)
        stack.extend([(kind,a,mid,depth+1),(kind,mid,b,depth+1)])
    for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
        p=sorted((F(c['a']),F(c['b'])) for c in records if c['variable']==kind)
        assert p[0][0]==lo and p[-1][1]==hi and all(a<b for a,b in p) and all(b==c for (a,b),(c,d) in zip(p,p[1:]))
    hashes[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    result=dict(status='PASS',witness_sha256=EXPECTED,spins=old['spins'],precision_bits=320,partition=records,
        certified_cases=[dict(b=str(b),delta_X4=str(d)) for b,d in CASES],reused_cells=len(old['partition'])-original_bad,
        initial_failed_cells=original_bad,refined_trials=trials,seconds=time.time()-tic,
        acceptance='Every stored old or fresh directed V_upper+b+d*X4_upper is strictly negative for all listed cases and every spin.',sha256_manifest=hashes)
    (OUT/'finite_v_positive_refined.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['partition','sha256_manifest']}),flush=True)
if __name__=='__main__':main()
