"""Adaptive interval proof of V/p<=-1/4, even42..110, compact mu19/4..49."""
from compact_finite220_combined import *
import argparse,json,time
OUT=ROOT/'output/boundary_length_20261007/compact/angle85'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def up(v):
    assert v.is_finite();m,e=v.upper().man_exp();return F(int(m))*F(2)**int(e)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bits',type=int,default=512);ap.add_argument('--lmax',type=int,default=110);ap.add_argument('--lmin',type=int,default=42)
    ap.add_argument('--mu-min',default='19/4');ap.add_argument('--mu-max',default='49');ap.add_argument('--maxevals',type=int,default=100000);ap.add_argument('--replay',type=Path);ap.add_argument('--partition-seed',type=Path);ap.add_argument('--out',type=Path);a=ap.parse_args()
    ctx.prec=a.bits;tic=time.time();spins=list(range(a.lmin,a.lmax+1,2));mom=OUT/'moments';bank=CompactFinite220(mom,spins)
    print('BUILT',time.time()-tic,flush=True);records={l:[]for l in spins}
    if a.replay:
        saved=json.loads(a.replay.read_text());assert saved['witness_sha256']==EXPECTED
        stack=[(F(c['a']),F(c['b']),[int(l)],0)for l,cells in saved['cells_by_spin'].items()for c in cells]
    elif a.partition_seed:
        # Only the exact interval geometry is reused. Every angle85-witness
        # value is freshly computed, and failing seed cells are subdivided.
        saved=json.loads(a.partition_seed.read_text());groups={}
        assert F(19,4)<=F(a.mu_min)<F(a.mu_max)<=49
        for ell in spins:
            for c in saved['cells_by_spin'][str(min(ell,max(map(int,saved['cells_by_spin']))))]:
                lo=max(F(a.mu_min),F(c['a']));hi=min(F(a.mu_max),F(c['b']))
                if lo<hi:groups.setdefault((lo,hi),[]).append(ell)
        stack=[(lo,hi,active,0)for(lo,hi),active in groups.items()]
    else:
        assert F(19,4)<=F(a.mu_min)<F(a.mu_max)<=49
        cuts=sorted(set([F(a.mu_min),F(a.mu_max)]+[F(x)for x in [5,8,12,24,48,49]if F(a.mu_min)<F(x)<F(a.mu_max)]));stack=[(x,y,spins.copy(),0)for x,y in zip(cuts,cuts[1:])]
    trials=0;worst=F(-10**1000)
    while stack:
        lo,hi,active,depth=stack.pop();trials+=1;assert trials<=a.maxevals
        mu=ar((lo+hi)/2)+arb(0,ar((hi-lo)/2));failed=[]
        try:values=bank.fields(mu,active)
        except (AssertionError,ValueError,ZeroDivisionError):values={}
        for l in active:
            v=values.get(l)
            if v is not None and v.is_finite() and v.upper()<-ar(F(1,4)):
                q=up(v);records[l].append(dict(a=str(lo),b=str(hi),V_over_p_upper_exact=str(q)));worst=max(worst,q)
            else:failed.append(l)
        if failed:
            assert depth<32,('UNRESOLVED',lo,hi,failed,{l:str(values.get(l))for l in failed})
            mid=(lo+hi)/2;stack.extend([(lo,mid,failed,depth+1),(mid,hi,failed,depth+1)])
        if trials%100==0:print(json.dumps(dict(trials=trials,pending=len(stack),accepted=sum(map(len,records.values())),last_mu=[str(lo),str(hi)],seconds=time.time()-tic)),flush=True)
    for l,cells in records.items():
        cells.sort(key=lambda c:F(c['a']));assert F(cells[0]['a'])==F(a.mu_min) and F(cells[-1]['b'])==F(a.mu_max)
        assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
    inputs=[WITNESS,Path(__file__),Path(__file__).with_name('compact_finite220_combined.py'),Path(old.__file__),ROOT/'research/all_direction_certification_20261007/head/forward220.py',ROOT/'research/all_direction_certification_20261007/head/forward_combined.py',ROOT/'research/crossing_kernels.py',ROOT/'output/boundary_length_20261007/tail/angle85/compact_norms.pkl']+list(mom.glob('*.json'))
    if a.partition_seed:inputs.append(a.partition_seed.resolve())
    result=dict(status='PASS_FINITE_COMPACT_BRIDGE',witness_sha256=EXPECTED,spins=spins,mu_domain=[str(F(a.mu_min)),str(F(a.mu_max))],kappa_exact='1/4',precision_bits=a.bits,
      cells_by_spin={str(l):v for l,v in records.items()},worst_upper_exact=str(worst),evaluations=trials,seconds=time.time()-tic,
      replay_input_sha256=sha(a.replay)if a.replay else None,partition_seed_sha256=sha(a.partition_seed)if a.partition_seed else None,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in inputs},
      scope='Uniform V/p<-1/4 for each listed spin and the exact recorded mu domain, using exact moment geometric remainders and tapered-PV endpoint formulas. No other spin/energy or real norm is included.')
    out=a.out or OUT/f'compact_bridge_combined_{a.bits}.json';out.write_text(json.dumps(result,indent=2));print('PASS',float(worst),trials,time.time()-tic,flush=True)
if __name__=='__main__':main()
