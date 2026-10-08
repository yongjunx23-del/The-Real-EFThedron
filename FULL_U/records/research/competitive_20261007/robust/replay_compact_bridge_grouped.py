"""Fresh higher-precision replay of every exact compact bridge interval.

Group identical spectral intervals to share forward evaluations. This does
not merge, omit, or enlarge any spin's original partition.
"""
from compact_finite220_combined import *
import argparse, json, time

OUT=ROOT/'output/competitive_20261007/compact'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upper(v):
    assert v.is_finite()
    m,e=v.upper().man_exp()
    return F(int(m))*F(2)**int(e)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,default=OUT/'compact_bridge_combined_512.json')
    ap.add_argument('--bits',type=int,default=640)
    ap.add_argument('--out',type=Path,default=OUT/'compact_bridge_combined_640.json')
    a=ap.parse_args();ctx.prec=a.bits;tic=time.time()
    saved=json.loads(a.input.read_text())
    assert saved['status']=='PASS_FINITE_COMPACT_BRIDGE'
    assert saved['witness_sha256']==EXPECTED and F(saved['kappa_exact'])==1
    assert saved['spins']==list(range(42,110,2))
    left,right=map(F,saved['mu_domain']);assert 5<=left<right<=49
    for p,h in saved['sha256_manifest'].items(): assert sha(ROOT/p)==h,p
    groups={};records={l:[] for l in saved['spins']}
    for l,rows in saved['cells_by_spin'].items():
        l=int(l);ordered=sorted(rows,key=lambda x:F(x['a']))
        assert F(ordered[0]['a'])==left and F(ordered[-1]['b'])==right
        assert all(F(x['b'])==F(y['a']) for x,y in zip(ordered,ordered[1:]))
        for c in ordered:
            lo,hi=F(c['a']),F(c['b'])
            assert lo<hi and F(c['V_over_p_upper_exact'])<-1
            groups.setdefault((lo,hi),[]).append(l)
    bank=CompactFinite220(OUT/'moments',saved['spins'])
    print('BUILT',time.time()-tic,'groups',len(groups),flush=True)
    worst=F(-10**1000)
    for i,((lo,hi),spins) in enumerate(sorted(groups.items()),1):
        assert len(set(spins))==len(spins)
        mu=ar((lo+hi)/2)+arb(0,ar((hi-lo)/2))
        values=bank.fields(mu,spins)
        for l in spins:
            v=values[l]
            assert v.is_finite() and v.upper()<-1,(l,lo,hi,str(v))
            q=upper(v);worst=max(worst,q)
            records[l].append(dict(a=str(lo),b=str(hi),V_over_p_upper_exact=str(q)))
        if i%100==0:print(json.dumps(dict(groups=i,total=len(groups),seconds=time.time()-tic)),flush=True)
    manifest=dict(saved['sha256_manifest'])
    for p in [a.input,Path(__file__)]:manifest[str(p.relative_to(ROOT))]=sha(p)
    out=dict(status='PASS_FINITE_COMPACT_BRIDGE',witness_sha256=EXPECTED,
        spins=saved['spins'],mu_domain=list(map(str,[left,right])),kappa_exact='1',precision_bits=a.bits,
        cells_by_spin={str(l):sorted(v,key=lambda x:F(x['a'])) for l,v in records.items()},
        worst_upper_exact=str(worst),evaluations=len(groups),seconds=time.time()-tic,
        replay_input_sha256=sha(a.input),sha256_manifest=manifest,
        every_original_spin_interval_replayed=True,
        scope='Fresh interval evaluation of all original exact spin/energy cells at higher precision; same mathematical implementation, no independent implementation claim.')
    a.out.write_text(json.dumps(out,indent=2));print('PASS',float(worst),time.time()-tic,flush=True)

if __name__=='__main__':main()
