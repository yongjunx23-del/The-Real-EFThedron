"""Fresh Arb320/380 direct evaluation of every b6 finite sign cell.

The final upper bound for each spin and cell is the exact maximum of its
two fresh direct enclosures. X4's exact increment is combined with its old
coefficient before interval evaluation. No assumed interval nesting or
floating comparison enters acceptance or the later support merge.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,hashlib,argparse
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parents[4];sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).parent))
from refine_finite import FiniteV,ctx,evaluate,EXPECTED,value,OUT,sha,B,D
from refine_positive_finite import arb,ar
def init():
    global BANKS
    BANKS={}
    for prec in [320,380]:
        ctx.prec=prec;bank=FiniteV();bank.lx+=ar(D);BANKS[prec]=bank
def batch(cells):
    out=[]
    for c in cells:
        a,b=F(c['a']),F(c['b']);pairs=[]
        for prec,bank in BANKS.items():
            ctx.prec=prec;x=arb(ar((a+b)/2),ar((b-a)/2))
            row=bank.fields(x)if c['variable']=='mu'else bank.uv_fields(x)
            assert set(row)==set(map(int,c['upper_by_spin']))
            upper={}
            for l,v in row.items():
                v+=ar(B);assert v.is_finite()
                m,e=v.upper().man_exp();q=F(int(m))*F(2)**int(e)
                assert q<0,('Sign failure',prec,c['variable'],c['a'],c['b'],l,float(q))
                upper[str(l)]=q
            pairs.append(upper)
        bounds={l:str(max(q[l]for q in pairs))for l in pairs[0]}
        assert all(value(v)<0 for v in c['upper_by_spin'].values())
        out.append(dict(variable=c['variable'],a=c['a'],b=c['b'],combined_upper_exact_by_spin=bounds))
    return out
def main(workers):
    tic=time.time();src=OUT/'finite_v.json';q=json.loads(src.read_text())
    assert q['status']=='PASS' and q['witness_sha256']==EXPECTED
    assert q['spins']==list(range(8,42,2)) and q['certified_cases']==[dict(b='6',delta_X4='107/5')]
    hashes=dict(q['sha256_manifest']);hashes[str(src.relative_to(ROOT))]=sha(src)
    for p,h in hashes.items():assert sha(ROOT/p)==h,p
    records=[];cells=q['partition'];chunks=[cells[i:i+100]for i in range(0,len(cells),100)]
    with ProcessPoolExecutor(max_workers=workers,initializer=init)as pool:
        for i,rows in enumerate(pool.map(batch,chunks)):
            records.extend(rows)
            if(i+1)%10==0:print('REPLAY',len(records),time.time()-tic,flush=True)
    assert len(records)==len(cells)
    for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
        ranges=sorted((F(c['a']),F(c['b']))for c in records if c['variable']==kind)
        assert ranges[0][0]==lo and ranges[-1][1]==hi and all(a<b for a,b in ranges)
        assert all(b==c for(a,b),(c,d)in zip(ranges,ranges[1:]))
    for p in [Path(__file__),Path(__file__).with_name('refine_finite.py'),ROOT/'research/full_direction_20261007/refine_positive_finite.py',ROOT/'research/directional_certification/certify_finite_vectors.py',ROOT/'research/tail_certification/certified_finite_v.py']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    out=dict(status='PASS_FRESH_ALL_CELL_REPLAY_AND_OUTWARD_MAXIMA',witness_sha256=EXPECTED,
      precision_bits=[320,380],direction_exact=['1','6'],crossing_increment_exact='107/5',spins=q['spins'],
      cells=len(records),all_spins_finite_and_negative=True,exact_partition_coverage=True,
      producer_sha256=sha(src),partition=records,seconds=time.time()-tic,sha256_manifest=hashes,
      proof='Each saved cell/spin bound is max(fresh320-bit direct upper,fresh380-bit direct upper). Both are strictly negative. Each direct evaluator first adds exact107/5 to the X4 coefficient and then adds exact6 to V/p. Every source cell was freshly reevaluated at both precisions; interval nesting is not assumed. Earlier separately bounded base/X4 cells serve only as valid partition seeds.')
    (OUT/'finite_v_replayed.json').write_text(json.dumps(out,indent=2));print('PASS',len(records),time.time()-tic,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);a=p.parse_args();main(a.workers)
