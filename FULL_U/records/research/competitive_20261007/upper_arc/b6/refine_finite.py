"""Validated V/p negativity for b6,d107/5, even8..40, every mu>4.

Unchanged b4 enclosures are reused only when their exact affine recombination
is strictly negative. Other cells are bisected and reevaluated with Arb320.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,time,json,hashlib
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[4];sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'research/full_direction_20261007'))
from refine_positive_finite import FiniteV,ctx,evaluate,frac,EXPECTED
OUT=ROOT/'output/competitive_20261007/upper_arc/b6'
B=F(6);D=F(107,5)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def value(q):return frac(q['V_upper'])+B+D*frac(q['X4_upper'])
def good(row):return row is not None and all(value(q)<-F(1,1000000) for q in row.values())
def main():
    tic=time.time();ctx.prec=320
    src=ROOT/'output/competitive_20261007/robust/b4/finite_v_b4.json';old=json.loads(src.read_text())
    assert old['status']=='PASS' and old['witness_sha256']==EXPECTED
    hashes=dict(old['sha256_manifest']);hashes[str(src.relative_to(ROOT))]=sha(src)
    for p,h in hashes.items():assert sha(ROOT/p)==h,p
    records=[];stack=[]
    for c in old['partition']:
        if good(c['upper_by_spin']):records.append(c)
        else:stack.append((c['variable'],F(c['a']),F(c['b']),0))
    bad=len(stack);print('BUILD',len(records),bad,flush=True);bank=FiniteV();trials=0;obstruction=None
    while stack:
        kind,a,b,depth=stack.pop();trials+=1
        assert trials<50000 and depth<24
        row=evaluate(bank,kind,a,b)
        if good(row):records.append(dict(variable=kind,a=str(a),b=str(b),upper_by_spin=row));continue
        mid=(a+b)/2
        if depth==0:
            point=evaluate(bank,kind,mid,mid)
            if not good(point):obstruction=dict(variable=kind,point=str(mid),upper_by_spin=point);break
        stack.extend([(kind,a,mid,depth+1),(kind,mid,b,depth+1)])
        if trials%100==0:print('REFINE',trials,len(stack),time.time()-tic,flush=True)
    hashes[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    result=dict(status='FAILED_SUFFICIENT_COERCIVITY' if obstruction else 'PASS',witness_sha256=EXPECTED,
      spins=old['spins'],precision_bits=320,partition=records if not obstruction else [],
      certified_cases=[dict(b=str(B),delta_X4=str(D))],initial_failed_cells=bad,refined_trials=trials,
      seconds=time.time()-tic,obstruction=obstruction,sha256_manifest=hashes)
    if not obstruction:
        for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
            p=sorted((F(c['a']),F(c['b']))for c in records if c['variable']==kind)
            assert p[0][0]==lo and p[-1][1]==hi and all(a<b for a,b in p) and all(b==c for(a,b),(c,d)in zip(p,p[1:]))
        result['exact_partition_coverage']=True
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'finite_v.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items()if k not in ['partition','sha256_manifest','obstruction']}),flush=True)
    if obstruction:print('OBSTRUCTION',str(obstruction)[:2500],flush=True)
if __name__=='__main__':main()
