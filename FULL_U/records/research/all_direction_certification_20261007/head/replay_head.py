"""Replay every rational cell and singular gap of an existing head component.

Stored dyadic bounds are not taken on trust: all interval fields and support
enclosures are recomputed. This is independent of adaptive cell selection,
but intentionally reuses the reviewed field/enclosure implementation.
"""
from finite220 import *
ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--moments',type=Path,required=True)
ap.add_argument('--record',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--ceiling',required=True)
a=ap.parse_args();ctx.prec=256;tic=time.time();r=json.loads(a.record.read_text())
assert hashlib.sha256(a.witness.read_bytes()).hexdigest()==r['witness_sha256']
for name,digest in r['dependencies'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
assert hashlib.sha256((a.moments/'manifest.json').read_bytes()).hexdigest()==r['source_manifest_sha256']
m=FiniteSpin(a.witness,r['spin'],a.moments,N=r.get('N',64))
gaps=[tuple(map(F,x)) for x in r['gap_beta_intervals']]
assert gaps==m.gaps(30)
cells=[(F(c['a']),F(c['b'])) for c in r['cells_data']]
partition=sorted(gaps+cells);assert partition[0][0]==0 and partition[-1][1]==1
assert all(x<y for x,y in partition)
assert all(a1==b0 for (a0,a1),(b0,b1) in zip(partition,partition[1:]))
total=sum((m.gap_bound(a,b) for a,b in gaps),arb(0));per=[]
for k,(aa,bb) in enumerate(cells):
    up,low=m.trapezoid(aa,bb);assert up.is_finite(),(aa,bb)
    total+=up
    if (k+1)%500==0:print('REPLAYED',k+1,'of',len(cells),time.time()-tic,flush=True)
passed=bool(total<arb(a.ceiling))
out=dict(witness_sha256=m.sha,source_record=str(a.record),source_record_sha256=hashlib.sha256(a.record.read_bytes()).hexdigest(),
    spin=m.ell,cells_replayed=len(cells),gaps_replayed=len(gaps),partition_covers_beta_zero_to_one_exactly=True,
    recomputed_upper=str(total.upper()),upper_dyadic=list(map(str,total.upper().man_exp())),ceiling=a.ceiling,passed=passed,
    elapsed_seconds=time.time()-tic,scope='Finite single-spin all-energy head component; no other spins or tail',
    dependency_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('finite220.py'),Path(__file__).with_name('forward220.py')]})
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if not passed:raise SystemExit(1)
