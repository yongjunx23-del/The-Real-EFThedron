"""Paired low-spin mixture integrations, with explicit one-thread process records."""
from pathlib import Path
import os,sys,json,subprocess,time,hashlib,concurrent.futures
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/six_targets_20261007/upper95'
def work(ell,bits):
    env=os.environ.copy()
    for k in['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
    env['PYTHONDONTWRITEBYTECODE']='1'
    cmd=[sys.executable,'research/six_targets_20261007/upper/upper95_head.py','--spin',str(ell),'--bits',str(bits),'--out',str(OUT/f'head_l{ell}_{bits}.json')]
    tic=time.time();log=OUT/f'head_l{ell}_{bits}.log'
    with log.open('w')as f:p=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
    r=dict(spin=ell,bits=bits,command=cmd,exit_code=p.returncode,elapsed_seconds=time.time()-tic,log=str(log.relative_to(ROOT)))
    (OUT/f'command_l{ell}_{bits}.json').write_text(json.dumps(r,indent=2));print(r,flush=True);assert p.returncode==0
    return r
if __name__=='__main__':
    tic=time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool:r=list(pool.map(lambda q:work(*q),[(e,b)for e in[0,2]for b in[256,384]]))
    (OUT/'run_manifest.json').write_text(json.dumps(dict(status='COMPLETE_PAIRED_LOW_SPIN_MIXTURE',jobs=r,elapsed_seconds=time.time()-tic,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
