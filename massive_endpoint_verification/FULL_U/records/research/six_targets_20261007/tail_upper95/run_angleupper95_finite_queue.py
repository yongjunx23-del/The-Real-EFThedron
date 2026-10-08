"""Two single-thread finite-spin proof/replay workers for angleupper95."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os,sys,subprocess,json,time
ROOT=Path(__file__).resolve().parents[3];CODE=Path(__file__).parent;OUT=ROOT/'output/six_targets_20261007/upper95_direct/tail'
GROUPS=[(42,90),(92,142),(144,192),(194,244),(246,294),(296,346),(348,396),(398,448)]
def main():
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
 def run(g):
  lo,hi=g;tic=time.time()
  for bits in[400,480]:
   target=OUT/f'finite_{lo}_{hi}_{bits}_k1_200.json'
   if target.exists():assert json.loads(target.read_text())['status']=='PASS_ANGLEupper95_FINITE_UV';continue
   args=['--precision',str(bits),'--spins',f'{lo}:{hi}']
   if bits==480:args+=['--replay',str(OUT/f'finite_{lo}_{hi}_400_k1_200.json')]
   log=OUT/f'finite_{lo}_{hi}_{bits}.log'
   with log.open('w')as f:r=subprocess.run([sys.executable,str(CODE/'certify_angleupper95_finite.py'),*args],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
   assert r.returncode==0 and target.exists()and json.loads(target.read_text())['status']=='PASS_ANGLEupper95_FINITE_UV',str(log)
  ans=dict(spins=[lo,hi],status='PASS_PRODUCER_AND_REPLAY',seconds=time.time()-tic);print(json.dumps(ans),flush=True);return ans
 with ThreadPoolExecutor(max_workers=2)as pool:results=list(pool.map(run,GROUPS))
 (OUT/'finite_queue_results.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
