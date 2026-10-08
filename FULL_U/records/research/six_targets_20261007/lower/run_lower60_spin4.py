"""Paired lower60 spin4 union-gap integration with retained commands/timings."""
from pathlib import Path
import os,sys,json,time,subprocess,hashlib,concurrent.futures
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/six_targets_20261007/lower/lower60_spin4'
OUT.mkdir(exist_ok=True);SOURCE=Path(__file__).with_name('mixture_head_lower60_spin4.py')
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def work(bits):
 path=OUT/f'head_l4_{bits}.json';assert not path.exists();log=path.with_suffix('.log')
 cmd=[sys.executable,str(SOURCE),'--case','lower60','--spin','4','--bits',str(bits),'--out',str(path)];tic=time.time()
 with log.open('w')as f:r=subprocess.run(cmd,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
 row=dict(case='lower60',spin=4,bits=bits,command=cmd,exit_code=r.returncode,elapsed_seconds=time.time()-tic,log=str(log.relative_to(ROOT)),output=str(path.relative_to(ROOT)))
 (OUT/f'command_l4_{bits}.json').write_text(json.dumps(row,indent=2));assert r.returncode==0,row;print(bits,'PASS',row['elapsed_seconds'],flush=True);return row
tic=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool:jobs=list(pool.map(work,[256,384]))
r=dict(status='COMPLETE_PAIRED_LOWER60_SPIN4_UNION_GAPS',jobs=jobs,elapsed_seconds=time.time()-tic,worker_count=2,source_sha256=sha(SOURCE),driver_sha256=sha(Path(__file__)),scope='Only spin4; full exact substitution into the existing lower60 certificate is required.')
(OUT/'campaign.json').write_text(json.dumps(r,indent=2));print(r['status'],r['elapsed_seconds'],flush=True)
