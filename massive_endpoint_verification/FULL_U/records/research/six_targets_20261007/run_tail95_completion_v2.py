"""Dependency-ordered95 tail assembly with immutable command records."""
from pathlib import Path
import time,sys,os,json,subprocess,hashlib
ROOT=Path(__file__).resolve().parents[2];CODE=ROOT/'research/six_targets_20261007';OUT=ROOT/'output/six_targets_20261007/upper95_direct'
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');runs=[]
def wait(paths):
 while not all(p.exists()for p in paths):time.sleep(3)
def run(label,args):
 command=[sys.executable,*map(str,args)];log=OUT/(label+'_v2.log');assert not log.exists();tic=time.time()
 with log.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
 row=dict(label=label,command=command,log=str(log.relative_to(ROOT)),returncode=r.returncode,seconds=time.time()-tic)
 runs.append(row)
 with(OUT/'tail_completion_commands_v2.jsonl').open('a')as f:f.write(json.dumps(row)+'\n')
 assert r.returncode==0,row;print(label,'PASS',row['seconds'],flush=True)
def main():
 tic=time.time();assert not(OUT/'tail_completion_commands_v2.jsonl').exists()
 wait([OUT/'tail/finite_queue_results.json',OUT/'tail/uniform_queue_results.json',OUT/'tail_auxiliary_campaign.json',CODE/'REVIEW_SIGNED_SOURCE2_UPPER95.md'])
 run('uv_assembly',[CODE/'assemble_uv_upper95_v2.py','--angle','upper95'])
 run('uv_verify',[CODE/'assemble_uv_upper95_v2.py','--angle','upper95','--verify-only'])
 wait([OUT/'compact/complete_compact_coercivity.json'])
 run('banded_cost',[CODE/'banded_upper95.py','--angle','upper95'])
 run('independent_tail_ledger',[CODE/'review_upper95_tail.py'])
 result=dict(status='PASS_ALL_UPPER95_TAIL_COMPONENTS_ASSEMBLED',seconds=time.time()-tic,jobs=runs,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (OUT/'tail_completion_campaign_v2.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
