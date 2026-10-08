"""Fresh replay on disjoint exact energy groups in independent processes."""
from compact_finite220_combined import *
import argparse,json,time,concurrent.futures

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upper(v):
    m,e=v.upper().man_exp();return F(int(m))*F(2)**int(e)
def worker(arg):
    bits,groups=arg;ctx.prec=bits;tic=time.time()
    bank=CompactFinite220(ROOT/'output/competitive_20261007/compact/moments',range(42,110,2))
    result=[]
    for aa,bb,spins in groups:
        a,b=F(aa),F(bb);mu=ar((a+b)/2)+arb(0,ar((b-a)/2))
        values=bank.fields(mu,spins)
        for ell in spins:
            v=values[ell];assert v.is_finite() and v.upper()<-1,(ell,a,b,str(v))
            result.append((ell,dict(a=aa,b=bb,V_over_p_upper_exact=str(upper(v)))))
    return result,time.time()-tic

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--bits',type=int,default=640)
    ap.add_argument('--jobs',type=int,default=4);a=ap.parse_args();tic=time.time()
    a.input=a.input.resolve();a.out=a.out.resolve()
    while True:
        try:j=json.loads(a.input.read_text());break
        except(FileNotFoundError,json.JSONDecodeError):time.sleep(2)
    assert j['status']=='PASS_FINITE_COMPACT_BRIDGE' and j['witness_sha256']==EXPECTED
    assert j['spins']==list(range(42,110,2)) and F(j['kappa_exact'])==1
    for p,h in j['sha256_manifest'].items():assert sha(ROOT/p)==h,p
    left,right=map(F,j['mu_domain']);groups={};original={}
    for ell,rows in j['cells_by_spin'].items():
        ell=int(ell);rows=sorted(rows,key=lambda c:F(c['a']))
        original[ell]=[(F(c['a']),F(c['b']))for c in rows]
        assert original[ell][0][0]==left and original[ell][-1][1]==right
        assert all(x[1]==y[0]for x,y in zip(original[ell],original[ell][1:]))
        for c in rows:
            lo,hi=F(c['a']),F(c['b']);assert lo<hi and F(c['V_over_p_upper_exact'])<-1
            groups.setdefault((str(lo),str(hi)),[]).append(ell)
    jobs=[[]for _ in range(a.jobs)]
    for i,((aa,bb),spins)in enumerate(groups.items()):jobs[i%a.jobs].append((aa,bb,spins))
    print('START',len(groups),'groups',a.jobs,'processes',flush=True)
    records={ell:[]for ell in j['spins']};times=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.jobs)as pool:
        for rows,elapsed in pool.map(worker,[(a.bits,g)for g in jobs]):
            for ell,row in rows:records[ell].append(row)
            times.append(elapsed)
    for ell,rows in records.items():
        rows.sort(key=lambda c:F(c['a']))
        assert [(F(c['a']),F(c['b']))for c in rows]==original[ell]
    manifest=dict(j['sha256_manifest'])
    for p in [a.input,Path(__file__)]:manifest[str(p.relative_to(ROOT))]=sha(p)
    out=dict(status='PASS_FINITE_COMPACT_BRIDGE',witness_sha256=EXPECTED,spins=j['spins'],
        mu_domain=j['mu_domain'],kappa_exact='1',precision_bits=a.bits,cells_by_spin={str(k):v for k,v in records.items()},
        worst_upper_exact=str(max(F(c['V_over_p_upper_exact'])for rows in records.values()for c in rows)),
        evaluations=len(groups),seconds=time.time()-tic,worker_seconds=times,
        replay_input_sha256=sha(a.input),every_original_spin_interval_replayed=True,sha256_manifest=manifest,
        scope='Fresh higher-precision replay of every original spin/energy cell, partitioned into disjoint groups across processes. Same mathematical evaluator.')
    a.out.write_text(json.dumps(out,indent=2));print('PASS',float(F(out['worst_upper_exact'])),out['seconds'],flush=True)
if __name__=='__main__':main()
