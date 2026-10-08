"""Finite-head producer and independent fixed-partition replay scheduling."""
from prepare_head import *
import concurrent.futures,argparse
def run(command,log):
    tic=time.time()
    with log.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
    entry=dict(command=command,log=str(log.relative_to(ROOT)),exit_code=r.returncode,elapsed_seconds=time.time()-tic)
    assert r.returncode==0,entry
    return entry
def work(row):
    ell=row['spin'];folder=BASE/f'spin{ell}';folder.mkdir(exist_ok=True)
    assert not (folder/'commands.json').exists()
    command=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'head_l{ell}.json'),'--spin',str(ell),'--model',row['model'],'--target',row['target'],'--max-cells','6500','--bits','256']
    runs=[run(command,folder/f'head_l{ell}.log')]
    command=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'replay_l{ell}.json'),'--replay',str(folder/f'head_l{ell}.json'),'--bits','384']
    runs.append(run(command,folder/f'replay_l{ell}.log'))
    (folder/'commands.json').write_text(json.dumps(runs,indent=2));print('REPLAYED',ell,flush=True)
    return dict(spin=ell,commands=runs)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);a=ap.parse_args()
    assert sha(WITNESS)==EXPECTED
    start=time.time();rows=json.loads((BASE/'budgets.json').read_text())['rows']
    assert not (BASE/'campaign.json').exists()
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers)as pool:jobs=list(pool.map(work,rows))
    (BASE/'jobs.json').write_text(json.dumps(jobs,indent=2))
    command=[sys.executable,str(OLD/'assemble.py'),'--witness',str(WITNESS),'--base',str(BASE)]
    assembly=run(command,BASE/'assembly.log')
    (BASE/'campaign.json').write_text(json.dumps(dict(status='COMPLETE_REPLAYED_FINITE_HEAD_ONLY',witness_sha256=EXPECTED,jobs=jobs,assembly=assembly,script_sha256=sha(__file__),elapsed_seconds=time.time()-start,scope='Every energy, even spins0..40 only. A same-witness all-spin tail is still required.'),indent=2))
    print((BASE/'head_summary.json').read_text(),flush=True)
if __name__=='__main__':main()
