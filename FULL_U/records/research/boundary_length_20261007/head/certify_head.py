"""Checkpointed arbitrary-direction interval head using reviewed220-row fields.

The adaptive rule is the frozen head integrator's rule. Checkpoint writes
are added only after a complete exact partition exists. A quality target
is separate from validity: every finite saved enclosure includes all cells
and singular caps even when that target has not yet been met.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/head'
sys.path.insert(0,str(OLD))
from finite_fast import *
import heapq

MODELS={'centered':FiniteCentered,'quartic':FiniteQuartic,'gap':FiniteGap,'fast':FiniteFast}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyadic(x):return list(map(str,x.upper().man_exp()))
def atomic_json(p,value):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix('.writing.json');tmp.write_text(json.dumps(value,indent=2));tmp.replace(p)


def code_hashes():
    names=['finite220.py','forward220.py','forward_combined.py','finite_combined.py','finite_taylor.py',
           'finite_centered.py','finite_global.py','finite_quartic.py','finite_gap.py','finite_fast.py']
    paths=[Path(__file__)]+[OLD/x for x in names]+[ROOT/'research'/x for x in
        ['head_certification_lowspin.py','head_certification_spin0.py','crossing_kernels.py',
         'certify_real_tail.py','certified_roy_moments.py','certified_strip_moments.py',
         'certified_uv_moments.py','validated_quadrature.py']]
    return {str(p.relative_to(ROOT)):sha(p)for p in paths}


def exact_partition(cells,gaps):
    union=sorted(cells+gaps)
    assert union[0][0]==0 and union[-1][1]==1
    assert all(a<b for a,b in union)
    assert all(x[1]==y[0]for x,y in zip(union,union[1:]))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True)
    ap.add_argument('--moments',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--spin',type=int);ap.add_argument('--target');ap.add_argument('--model',choices=MODELS,default='fast')
    ap.add_argument('--max-cells',type=int,default=10000);ap.add_argument('--bits',type=int,default=256)
    ap.add_argument('--replay',type=Path);ap.add_argument('--checkpoint-every',type=int,default=500)
    a=ap.parse_args();ctx.prec=a.bits;tic=time.time();deps=code_hashes()
    assert not a.out.exists(),'Use a new output path; never overwrite a frozen record.'
    manifest=sha(a.moments/'manifest.json');witness=sha(a.witness)
    if a.replay:
        source=json.loads(a.replay.read_text())
        assert source['witness_sha256']==witness and source['source_manifest_sha256']==manifest
        for name,digest in source['dependencies'].items():assert sha(ROOT/name)==digest,name
        ell=source['spin'];kind=source['model'];target=source['target']
        model=MODELS[kind](a.witness,ell,a.moments,N=128)
        assert [str(x)for x in model.base.direction]==source['direction_exact']
        gaps=[tuple(map(F,x))for x in source['gap_beta_intervals']]
        assert gaps==model.gaps(30)
        cells=[(F(c['a']),F(c['b']))for c in source['cells_data']]
        exact_partition(cells,gaps)
        total=sum((model.gap_bound(lo,hi)for lo,hi in gaps),arb(0))
        for k,(lo,hi)in enumerate(cells):
            up,_=model.trapezoid(lo,hi);assert up.is_finite();total+=up
            if(k+1)%500==0:print('REPLAY',ell,k+1,'of',len(cells),time.time()-tic,flush=True)
        # A replay encloses the whole integral even if it has a wider numerical
        # value than a prior quality target. Its actual directed sum is retained.
        result=dict(status='PASS_FULL_PARTITION_REPLAY',witness_sha256=witness,spin=ell,model=kind,
            source_record=str(a.replay),source_record_sha256=sha(a.replay),source_manifest_sha256=manifest,
            direction_exact=source['direction_exact'],cells_replayed=len(cells),gaps_replayed=len(gaps),
            partition_covers_beta_zero_to_one_exactly=True,passed=True,target=target,
            target_met=bool(total<arb(target)),recomputed_upper=str(total.upper()),upper_dyadic=dyadic(total),
            elapsed_seconds=time.time()-tic,precision_bits=a.bits,dependency_sha256=deps,
            scope='Only this fixed spin, every energy; no omitted-spin bound included.')
        atomic_json(a.out,result);print(json.dumps(result,indent=2),flush=True);return
    assert a.spin is not None and a.target is not None
    ell=a.spin;kind=a.model;model=MODELS[kind](a.witness,ell,a.moments,N=128)
    gaps=model.gaps(30);gapmass=sum((model.gap_bound(lo,hi)for lo,hi in gaps),arb(0))
    upper=gapmass;estimate=arb(0);heap=[];serial=0
    cut=F(math.ceil(math.sqrt(1-4/float(model.uv_start))*2**20),2**20)
    assert F(4)/(1-cut*cut)>=model.uv_start
    points=sorted(set([F(i,256)for i in range(257)]+[cut]+[v for g in gaps for v in g]))
    def insert(lo,hi):
        nonlocal serial,upper,estimate
        up,mid=model.trapezoid(lo,hi)
        if not up.is_finite():
            if hi-lo<F(1,2**45):raise RuntimeError(f'Nonfinite cell {lo} {hi}')
            c=(lo+hi)/2;insert(lo,c);insert(c,hi);return
        upper+=up;estimate+=mid
        heapq.heappush(heap,(-float(up-mid),serial,lo,hi,up,mid));serial+=1
    def save(final=False):
        cells=sorted(heap,key=lambda q:q[2]);exact_partition([(q[2],q[3])for q in cells],gaps)
        result=dict(status='PASS_FINITE_SPIN_ENCLOSURE',final=final,witness_sha256=witness,
            direction_exact=[str(x)for x in model.base.direction],spin=ell,model=kind,N=128,precision_bits=a.bits,
            upper=str(upper.upper()),upper_dyadic=dyadic(upper),midpoint_estimate=str(estimate),
            target=a.target,target_met=bool(upper<arb(a.target)),cells=len(cells),
            gap_upper=str(gapmass.upper()),gap_upper_dyadic=dyadic(gapmass),
            gap_beta_intervals=[[str(x),str(y)]for x,y in gaps],threshold=model.cap_status,
            source_endpoint_caps=getattr(model,'gap_records',[]),source_manifest_sha256=manifest,
            elapsed_seconds=time.time()-tic,dependencies=deps,
            cells_data=[dict(a=str(q[2]),b=str(q[3]),upper_dyadic=dyadic(q[4]))for q in cells],
            claim='Rigorous single-spin all-energy enclosure; target is quality only. Other spins and tails are separate.')
        atomic_json(a.out,result)
    for lo,hi in zip(points,points[1:]):
        if any(lo>=ga and hi<=gb for ga,gb in gaps):continue
        insert(lo,hi)
    save();print('INITIAL',ell,len(heap),str(upper.upper()),'target',a.target,flush=True)
    next_checkpoint=((len(heap)//a.checkpoint_every)+1)*a.checkpoint_every
    while len(heap)<a.max_cells and not upper<arb(a.target):
        _,_,lo,hi,up,mid=heapq.heappop(heap);upper-=up;estimate-=mid
        center=(lo+hi)/2;insert(lo,center);insert(center,hi)
        if len(heap)>=next_checkpoint:
            save();print('CHECKPOINT',ell,len(heap),float(upper.upper()),time.time()-tic,flush=True)
            next_checkpoint+=a.checkpoint_every
    save(final=True);print('COMPLETE',ell,str(upper.upper()),bool(upper<arb(a.target)),time.time()-tic,flush=True)


if __name__=='__main__':main()
