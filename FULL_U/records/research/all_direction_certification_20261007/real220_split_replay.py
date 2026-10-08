"""Directed real-norm enclosure with an exact compact/UV partition and replay.

The output contains real-norm components, not a complete disk certificate.
Coercivity and finite heads must be established separately for the same witness.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import argparse, hashlib, heapq, json, time
from real220 import Norm, ctx, endpoint
from real220_finite import Finite


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--witness', type=Path, required=True)
    ap.add_argument('--mode', choices=['finite', 'tail'], required=True)
    ap.add_argument('--L', type=int, default=448)
    ap.add_argument('--cells', type=int, default=4096)
    ap.add_argument('--bits', type=int, default=256)
    ap.add_argument('--replay-bits', type=int, default=384)
    ap.add_argument('--uv-weight', type=int, default=1000)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    tic = time.time()
    assert args.replay_bits > args.bits and args.uv_weight >= 1
    # Rational lower bracket for sqrt(45/49); compact cells lie within mu<=49.
    denominator = 2**48
    cutoff = F(isqrt(45*denominator**2//49), denominator)
    assert cutoff**2 <= F(45,49) < (cutoff + F(1,denominator))**2
    def build():
        return Finite(args.witness,42,args.L,20) if args.mode=='finite' else Norm(args.witness,args.L,16)
    ctx.prec = args.bits
    model = build()
    heap = []
    totals = {'compact': F(0), 'remainder': F(0)}
    ident = 0
    def insert(lo,hi,region):
        nonlocal ident
        value = model.density(lo,hi)
        if not value.is_finite():
            assert hi-lo > F(1,2**30)
            mid=(lo+hi)/2
            insert(lo,mid,region)
            insert(mid,hi,region)
            return
        mass = max(F(0),endpoint(value))*(hi-lo)
        weight = 1 if region=='compact' else args.uv_weight
        totals[region] += mass
        heapq.heappush(heap,(-float(weight*mass),ident,lo,hi,mass,region))
        ident += 1
    knots = sorted(set([F(i,64) for i in range(65)]+[cutoff]))
    for lo,hi in zip(knots,knots[1:]):
        insert(lo,hi,'compact' if hi<=cutoff else 'remainder')
    while len(heap)<args.cells:
        _,_,lo,hi,mass,region = heapq.heappop(heap)
        totals[region]-=mass
        mid=(lo+hi)/2
        insert(lo,mid,region)
        insert(mid,hi,region)
        if len(heap)%512==0:
            print('cells',len(heap),'A', {k:float(v) for k,v in totals.items()},'seconds',time.time()-tic,flush=True)
    ctx.prec=args.replay_bits
    model=build()
    replay_totals={'compact':F(0),'remainder':F(0)}
    retained_totals={'compact':F(0),'remainder':F(0)}
    records=[]
    cursor=F(0)
    for i,(_,_,lo,hi,mass,region) in enumerate(sorted(heap,key=lambda r:r[2])):
        assert cursor==lo<hi<=1
        assert (hi<=cutoff) if region=='compact' else (lo>=cutoff)
        value=model.density(lo,hi)
        assert value.is_finite()
        replay=max(F(0),endpoint(value))*(hi-lo)
        retained=max(mass,replay)
        replay_totals[region]+=replay
        retained_totals[region]+=retained
        records.append(dict(lo=str(lo),hi=str(hi),region=region,
            first_upper_exact=str(mass),replay_upper_exact=str(replay),retained_upper_exact=str(retained)))
        cursor=hi
        if (i+1)%1024==0:
            print('replayed',i+1,'seconds',time.time()-tic,flush=True)
    assert cursor==1
    out=dict(status='PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY',
        witness=str(args.witness.resolve()),witness_sha256=sha(args.witness),mode=args.mode,
        first_even_spin=42 if args.mode=='finite' else args.L+2,
        last_even_spin=args.L if args.mode=='finite' else 'infinity',
        reference_order=20 if args.mode=='finite' else 16,
        precision_bits=args.bits,replay_precision_bits=args.replay_bits,
        compact_beta_cutoff_exact=str(cutoff),compact_mu_cutoff_exact=str(4/(1-cutoff**2)),
        A_upper_exact={k:str(v) for k,v in retained_totals.items()},
        A_upper_float={k:float(v) for k,v in retained_totals.items()},
        first_upper_exact={k:str(v) for k,v in totals.items()},
        replay_upper_exact={k:str(v) for k,v in replay_totals.items()},
        sha256_manifest={str(p.resolve()):sha(p) for p in [Path(__file__),Path(__file__).with_name('real220.py'),Path(__file__).with_name('real220_finite.py')]},
        cells=records,seconds=time.time()-tic,
        scope='All energy real-norm component. Compact is beta<=rational cutoff strictly below mu49; remainder includes a tiny compact sliver and all UV. Requires separate V/p coercivity over both sets and finite-head bound. Precision replay validates the enclosure arithmetic, not an independent derivation of the mathematics.')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k not in ['cells','sha256_manifest','first_upper_exact','replay_upper_exact','A_upper_exact']},indent=2),flush=True)


if __name__=='__main__':
    main()
