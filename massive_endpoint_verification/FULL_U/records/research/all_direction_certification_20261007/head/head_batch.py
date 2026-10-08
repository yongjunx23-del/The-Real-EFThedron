"""Checkpointed per-spin interval integration for fixed 220-row witnesses."""
from finite220 import *
ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--moments',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
ap.add_argument('--numerical-costs',type=Path,required=True);ap.add_argument('--start',type=int,default=2);ap.add_argument('--stop',type=int,default=40)
ap.add_argument('--relative-slack',type=float,default=.03);ap.add_argument('--absolute-slack',type=float,default=.0001);ap.add_argument('--max-cells',type=int,default=4000)
a=ap.parse_args();ctx.prec=256;tic=time.time();costs=json.loads(a.numerical_costs.read_text())['per_spin'];a.out.mkdir(parents=True,exist_ok=True)
code=[Path(__file__),Path(__file__).with_name('finite220.py'),Path(__file__).with_name('forward220.py'),ROOT/'research/head_certification_lowspin.py',ROOT/'research/head_certification_spin0.py']
deps={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in code}
records=[]
for ell in range(a.start,a.stop+1,2):
    out=a.out/f'head_l{ell}.json'
    if out.exists():records.append(json.loads(out.read_text()));continue
    ts=time.time();m=FiniteSpin(a.witness,ell,a.moments,N=128);target=str(costs[str(ell)]*(1+a.relative_slack)+a.absolute_slack)
    print('START_SPIN',ell,'target',target,'build_s',time.time()-ts,flush=True)
    r=integrate(m,target,a.max_cells,30)
    r.update(spin=ell,witness_sha256=m.sha,threshold=m.cap_status,N=128,elapsed_seconds=time.time()-ts,dependencies=deps,
        source_manifest_sha256=hashlib.sha256((a.moments/'manifest.json').read_bytes()).hexdigest(),
        claim='Rigorous fixed-witness finite single-spin all-energy upper bound only; remaining spins and tails separate')
    out.write_text(json.dumps(r,indent=2));records.append({k:v for k,v in r.items() if k!='cells_data'})
    (a.out/'batch_summary.json').write_text(json.dumps(dict(records=records,seconds=time.time()-tic,complete=ell==a.stop),indent=2))
    print('DONE_SPIN',ell,r['upper'],r['target_met'],'seconds',r['elapsed_seconds'],flush=True)
