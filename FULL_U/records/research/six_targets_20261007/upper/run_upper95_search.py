"""Record a fresh robust220 search at the requested95-percent upper normal."""
from pathlib import Path
import os,sys,json,subprocess,time,hashlib
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/six_targets_20261007/upper95_search'
if __name__=='__main__':
 OUT.mkdir(parents=True,exist_ok=True);assert not (OUT/'command.json').exists()
 guards=OUT/'threshold_guards.json';guards.write_text(json.dumps([dict(energies=[4.00000001,4.000001,4.0001],spins=[2,4],margin=.001)],indent=2))
 target=next(t for t in json.loads((ROOT/'output/six_targets_20261007/target_map.json').read_text())['targets']if t['tag']=='upper95')
 code='research/full_direction_20261007/right/search_upper_g2_cached_guards.py'
 command=[sys.executable,code,'--angle',str(target['normal_angle_degrees']),'--warm','output/weak_regions_20261007/upper/angle83_repaired.npz','--out',str(OUT/'candidate'),'--cache','output/full_direction_20261007/right/dense220_n192k96_cache.npz','--cached-guards','output/full_direction_20261007/left/dense_v_guards_l10_40_mu6p005_20p005_k96.npz','--cached-guard-margin','.002','--local-guards',str(guards),'--tail-margin','.01','--source-margin','.0001','--max-iter','20','--constraint-generation']
 env=os.environ.copy()
 for k in['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[k]='1'
 env['PYTHONDONTWRITEBYTECODE']='1';tic=time.time()
 with (OUT/'search.log').open('w')as f:p=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
 r=dict(command=command,exit_code=p.returncode,elapsed_seconds=time.time()-tic,source_sha256=hashlib.sha256((ROOT/code).read_bytes()).hexdigest(),environment={k:env[k]for k in['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','PYTHONDONTWRITEBYTECODE']})
 (OUT/'command.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2));assert p.returncode==0
