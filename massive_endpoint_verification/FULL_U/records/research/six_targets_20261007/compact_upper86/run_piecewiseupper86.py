from pathlib import Path
import sys
_ADAPT_ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(_ADAPT_ROOT/'research/six_targets_20261007/tail_upper86'),str(_ADAPT_ROOT/'research/boundary_length_20261007/compact')]
"""Exact local-spin compactupper86 bridges, maximum four thread-limited workers."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import subprocess,sys,concurrent.futures,time,json,threading
LOCK=threading.Lock()
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/six_targets_20261007/upper86/compact'
def run_logged(cmd,log):
    started=time.time()
    with LOCK:
        with (OUT/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(event='start',command=cmd,log=str(log.relative_to(ROOT)),unix_seconds=started))+'\n')
    with log.open('w')as f:
        done=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    with LOCK:
        with (OUT/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(event='complete',command=cmd,returncode=done.returncode,seconds=time.time()-started))+'\n')
    done.check_returncode()

def work(row):
    aa,bb=row['mu_domain_exact'];L=row['minimum_even_spin'];aa=str(F(aa))
    if L==42:return
    tag='piece_'+aa.replace('/','d')+'_'+bb.replace('/','d')
    producer=OUT/f'{tag}_512.json';replay=OUT/f'{tag}_640.json'
    cmd=[sys.executable,str(Path(__file__).with_name('certify_finiteupper86_slice.py')),'--mu-min',aa,'--mu-max',bb,'--lmax',str(L-2),'--bits','512','--out',str(producer)]
    if F(aa)>=5:cmd+=['--partition-seed',str(ROOT/'output/competitive_20261007/compact/compact_bridge_sliced_512.json')]
    if not producer.exists():
        run_logged(cmd,OUT/f'{tag}_512.log')
    print('PRODUCER PASS',tag,flush=True)
    cmd=[sys.executable,str(Path(__file__).with_name('replay_piecewiseupper86.py')),'--input',str(producer),'--out',str(replay),'--bits','640','--jobs','1']
    if not replay.exists():
        run_logged(cmd,OUT/f'{tag}_640.log')
    print('REPLAY PASS',tag,flush=True)
if __name__=='__main__':
    for kind in['roy','strips']:
        while not (OUT/f'moments/{kind}_manifest.json').exists():time.sleep(2)
    run_logged([sys.executable,str(Path(__file__).with_name('validate_compactupper86_combined.py')),'--spins','42,60'],OUT/'point_check.log')
    print('POINT COMPARISON PASS',flush=True)
    plan=json.loads((OUT/'analytic_piecewise_512.json').read_text())['records']
    with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool:list(pool.map(work,plan))
