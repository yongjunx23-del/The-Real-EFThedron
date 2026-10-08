"""Replay frozen fixed-witness components for additional exact rational slopes.

This file never writes frozen 2026-10-06 records.  The imported implementations
are bound in every emitted record; exact slopes and repair ratios are explicit.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, sys, time
ROOT=Path(__file__).resolve().parents[2]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/directional_certification'),str(ROOT/'research')]
OUT=ROOT/'output/full_direction_20261007/certification'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cases(slopes,ratio):return [(F(b),F(0) if F(b)<=0 else F(ratio)*F(b)) for b in slopes.split(',')]

def tail(args):
    import certify_directional_tail as mod
    mod.CASES=cases(args.slopes,args.ratio)
    out=OUT/args.name
    old=sys.argv;sys.argv=[__file__,'--output',str(out)]
    try:mod.main()
    finally:sys.argv=old
    d=json.loads(out.read_text());d['repair_ratio_exact']=str(F(args.ratio))
    d['actual_cases_exact']=[[str(b),str(c)] for b,c in mod.CASES]
    d['sha256_manifest'][str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    out.write_text(json.dumps(d,indent=2))

def finite(args):
    # The reviewed exact partition merge accepts arbitrary nonnegative X4
    # increments. Intercept only its single output pathname, never the inputs.
    import combine_finite_family as mod
    mod.CASES=cases(args.slopes,args.ratio)
    assert all(d>=0 for b,d in mod.CASES)
    original=Path.write_text;target=OUT/args.name
    def redirected(path,text,*a,**k):
        if path==ROOT/'output/directional_certification/finite_block_family.json':
            data=json.loads(text)
            data['repair_ratio_exact']=str(F(args.ratio))
            data['sha256_manifest'][str(Path(__file__).relative_to(ROOT))]=sha(__file__)
            return original(target,json.dumps(data,indent=2),*a,**k)
        raise RuntimeError('Unexpected write attempted: '+str(path))
    Path.write_text=redirected
    try:mod.main()
    finally:Path.write_text=original

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('component',choices=['tail','finite'])
    p.add_argument('--slopes',default='-2,-4,-8,-16,3/10,7/20,2/5,1/2')
    p.add_argument('--ratio',default='18/5');p.add_argument('--name',required=True)
    a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    (tail if a.component=='tail' else finite)(a)
