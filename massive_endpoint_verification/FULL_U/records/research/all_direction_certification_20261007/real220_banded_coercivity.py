"""Combine existing exact real-norm cells with proved energy-dependent margins.

Uses the minimum local kappa over all finite spins on each beta cell. This
is cheaper but weaker than the per-spin weighted integration.
"""
from pathlib import Path
from fractions import Fraction as F
import json,argparse,time
from real220_local_coercivity import Coercivity,ROOT,SHA,sha

def dyadic_ceiling(x,bits=40):
    assert x>0
    e=x.numerator.bit_length()-x.denominator.bit_length()
    if F(2)**e>x:e-=1
    step=F(2)**(e-bits+1);n=-((-x/step).numerator//(-x/step).denominator)
    out=n*step;assert x<=out<x+step
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--descriptor',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();tic=time.time();co=Coercivity(args.descriptor.resolve())
    source=ROOT/'output/competitive_20261007/left/real/finite42_448_split_replay.json';d=co.read(source)
    assert d['status']=='PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY' and d['witness_sha256']==SHA
    assert d['mode']=='finite'and d['first_even_spin']==42 and d['last_even_spin']==448
    assert d['replay_precision_bits']>d['precision_bits'];cursor=F(0);total=F(0);records=[];sums={'compact':F(0),'remainder':F(0)}
    for c in d['cells']:
        lo,hi=F(c['lo']),F(c['hi']);assert cursor==lo<hi<=1
        first,replay,retained=map(F,[c['first_upper_exact'],c['replay_upper_exact'],c['retained_upper_exact']])
        assert 0<=first and 0<=replay and retained==max(first,replay)
        kappas=[co.kappa(ell,lo,hi)for ell in range(42,449,2)];kap=min(kappas)
        weight=dyadic_ceiling(1/(2*kap));cost=retained*weight
        total+=cost;sums[c['region']]+=retained;cursor=hi
        records.append(dict(lo=str(lo),hi=str(hi),A_upper_exact=str(retained),kappa_exact=str(kap),reciprocal_weight_upper_exact=str(weight),cost_upper_exact=str(cost)))
    assert cursor==1 and all(sums[k]==F(d['A_upper_exact'][k])for k in sums)
    for p in [Path(__file__),Path(__file__).with_name('real220_local_coercivity.py')]:co.manifest[str(p.relative_to(ROOT))]=sha(p)
    out=dict(status='PASS_BANDED_COERCIVITY_REAL_COST',witness_sha256=SHA,objective_exact=['-1','0'],first_even_spin=42,last_even_spin=448,
        energy_domain=['4','infinity'],cost_upper_exact=str(total),cost_upper_float=float(total),cells=records,sha256_manifest=co.manifest,seconds=time.time()-tic,
        reciprocal_rounding='1/(2*kappa) is rounded upward to40 significant binary bits, with an exact rational inequality check; products and sums then remain dyadic.',
        scope='Complete finite42..448 allenergy support upper; minimum certified kappa across all spins on every existing beta cell. Conservative relative to separate per-spin integration. Requires separate finitehead and infinite-spin cost.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items()if k not in ['cells','sha256_manifest','cost_upper_exact']},indent=2))
if __name__=='__main__':main()
