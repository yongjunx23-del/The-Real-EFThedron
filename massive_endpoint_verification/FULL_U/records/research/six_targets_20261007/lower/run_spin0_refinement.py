"""Four paired-precision spin0 refinement jobs with identical fixed geometry."""
from pathlib import Path
import os,sys,json,time,subprocess,hashlib,concurrent.futures
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/six_targets_20261007/lower/spin0_refined';OUT.mkdir(exist_ok=True)
SOURCE=Path(__file__).with_name('mixture_head_spin0_refined.py');ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def work(job):
 case,bits=job;folder=OUT/case;folder.mkdir(exist_ok=True);path=folder/f'head_l0_{bits}.json';assert not path.exists()
 command=[sys.executable,str(SOURCE),'--case',case,'--spin','0','--bits',str(bits),'--out',str(path)];log=path.with_suffix('.log');tic=time.time()
 with log.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
 row=dict(case=case,spin=0,bits=bits,command=command,exit_code=r.returncode,elapsed_seconds=time.time()-tic,log=str(log.relative_to(ROOT)),output=str(path.relative_to(ROOT)))
 (folder/f'command_l0_{bits}.json').write_text(json.dumps(row,indent=2));assert r.returncode==0,row;print(case,bits,row['elapsed_seconds'],flush=True);return row
tic=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=4)as pool:jobs=list(pool.map(work,[(c,b)for c in['lower15','lower45']for b in[256,384]]))
r=dict(status='COMPLETE_PAIRED_SPIN0_SUBDIVISION',jobs=jobs,elapsed_seconds=time.time()-tic,worker_count=4,source_sha256=sha(SOURCE),driver_sha256=sha(Path(__file__)),scope='Same exact positive mixture and gap bounds; every ordinary spin0 cell bisected, at two precisions. Complete assembly remains required.')
(OUT/'campaign.json').write_text(json.dumps(r,indent=2));print(r['status'],r['elapsed_seconds'],flush=True)
