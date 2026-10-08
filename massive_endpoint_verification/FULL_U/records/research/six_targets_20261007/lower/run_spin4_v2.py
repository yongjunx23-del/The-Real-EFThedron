"""Two independent spin4 jobs, each recomputed at 256 and 384 bits."""
from pathlib import Path
import os,sys,json,time,subprocess,hashlib,concurrent.futures
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/six_targets_20261007/lower'
SOURCE=Path(__file__).with_name('mixture_head_spin4_v2.py')
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def work(job):
 case,spin=job;folder=OUT/case;commands=[]
 for bits in[256,384]:
  path=folder/f'head_l{spin}_{bits}.json';assert not path.exists()
  command=[sys.executable,str(SOURCE),'--case',case,'--spin',str(spin),'--bits',str(bits),'--out',str(path)]
  log=path.with_suffix('.v2.log');start=time.time()
  with log.open('w')as f:result=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
  row=dict(command=command,exit_code=result.returncode,elapsed_seconds=time.time()-start,log=str(log.relative_to(ROOT)),output=str(path.relative_to(ROOT)))
  commands.append(row);(folder/f'commands_l{spin}_v2.json').write_text(json.dumps(commands,indent=2))
  assert result.returncode==0,row
  print(case,spin,bits,'PASS',row['elapsed_seconds'],flush=True)
 return dict(case=case,spin=spin,commands=commands)
start=time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool:rows=list(pool.map(work,[(case,ell)for case in['lower15','lower45']for ell in[4]]))
report=dict(status='COMPLETE_PAIRED_PRECISION_MIXTURE_SPIN4',jobs=rows,elapsed_seconds=time.time()-start,worker_count=2,environment={k:ENV[k]for k in['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']},source_sha256=sha(SOURCE),driver_sha256=sha(Path(__file__)),scope='Only the selected spins. Full certificates require additive decomposition and inherited base tails.')
(OUT/'head4_campaign_v2.json').write_text(json.dumps(report,indent=2));print(report['status'],report['elapsed_seconds'],flush=True)
