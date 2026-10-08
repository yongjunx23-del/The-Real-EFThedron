"""Replay new b4 finite-spin sign cells; bind reused old cells exactly.

The finite support merge uses the original rigorous producer enclosures.
This independent precision replay checks sign and enclosure agreement, not
the nesting of all individual interval endpoints.
"""
import json,sys,time,hashlib
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/full_direction_20261007'))
from refine_positive_finite import FiniteV,ctx,evaluate,frac,EXPECTED
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def key(c):return c['variable'],F(c['a']),F(c['b'])
def value(v):return frac(v['V_upper'])+4+F(72,5)*frac(v['X4_upper'])
def main():
    tic=time.time();ctx.prec=380;out=ROOT/'output/competitive_20261007/robust/b4'
    oldpath=ROOT/'output/full_direction_20261007/certification/finite_v_positive_refined.json';newpath=out/'finite_v_b4.json'
    old,new=[json.loads(p.read_text())for p in [oldpath,newpath]]
    assert old['status']==new['status']=='PASS' and old['witness_sha256']==new['witness_sha256']==EXPECTED
    oldrows={key(c):c['upper_by_spin']for c in old['partition']};bank=FiniteV();reused=fresh=0;worst=None;maxdiff=F(0);maxovershoot=F(0);maxrelative=F(0)
    for c in new['partition']:
        k=key(c)
        if k in oldrows and oldrows[k]==c['upper_by_spin']:
            assert all(value(v)<0 for v in c['upper_by_spin'].values());reused+=1;continue
        row=evaluate(bank,*k);assert row is not None and set(row)==set(c['upper_by_spin']);fresh+=1
        for l,v in row.items():
            n=value(v);o=value(c['upper_by_spin'][l]);assert n<0 and o<0
            worst=n if worst is None else max(worst,n);maxdiff=max(maxdiff,abs(n-o));maxovershoot=max(maxovershoot,n-o)
            maxrelative=max(maxrelative,abs(n-o)/max(F(1),abs(n),abs(o)))
    # Interval upper endpoints need not nest under a precision change: the
    # Arb radius precision and interval expression cancellation are separate
    # from midpoint precision. Sign is the replay assertion; differences are
    # diagnostic only and are not used to strengthen the original bound.
    paths=[oldpath,newpath,Path(__file__),ROOT/'research/full_direction_20261007/refine_positive_finite.py',ROOT/'research/tail_certification/certified_finite_v.py']
    result=dict(status='PASS_FRESH_SIGN_REPLAY_AND_EXACT_REUSE',precision_bits=ctx.prec,witness_sha256=EXPECTED,
      exact_reused_cells=reused,fresh_cells=fresh,all_fresh_spins_finite_and_negative=True,maximum_fresh_V_over_p_upper_exact=str(worst),
      maximum_combined_upper_difference_exact=str(maxdiff),maximum_combined_upper_overshoot_exact=str(maxovershoot),
      maximum_combined_upper_relative_difference_exact=str(maxrelative),
      endpoint_agreement_required=False,
      certificate_sha256=sha(newpath),seconds=time.time()-tic,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in paths},
      scope='Fresh380-bit sign replay on all newly refined cells; unchanged cells match a previously certified input exactly. The support merge retains the original320-bit directed enclosures; this record does not assert all interval endpoint nesting.')
    (out/'finite_v_replay.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
