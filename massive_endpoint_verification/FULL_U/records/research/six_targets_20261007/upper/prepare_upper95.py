"""Pin the legacy source cache after a fresh exact rational rebuild, without altering it."""
from pathlib import Path
import os,sys,json,pickle,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'research'))
from upper95_common import *
from head_certification_spin0 import source_cache
import numpy as np
if __name__=='__main__':
    tic=time.time();OUT.mkdir(parents=True,exist_ok=True);assert not (OUT/'preparation.json').exists()
    row=next(t for t in json.loads((ROOT/'output/six_targets_20261007/target_map.json').read_text())['targets']if t['tag']=='upper95')
    np.savez(OUT/'objective.npz',direction=np.array(row['physical_outward_unit_normal']))
    p=ROOT/'output/certification_head/spin0_source_N64.pkl';before=sha(p);saved=pickle.loads(p.read_bytes())
    code=ROOT/'research/head_certification_spin0.py'
    assert saved['witness']==sha(ROOT/WITNESSES[1])and saved['builder']==sha(code)
    rebuilt=source_cache(ROOT/WITNESSES[1],64);assert rebuilt==saved['data'];assert sha(p)==before
    manifest={str(q.relative_to(ROOT)):sha(q)for q in[p,code,Path(__file__),Path(__file__).with_name('upper95_common.py'),ROOT/WITNESSES[1],OUT/'objective.npz',ROOT/'output/six_targets_20261007/target_map.json']}
    result=dict(status='PASS_EXACT_LEGACY_SOURCE_CACHE_REBUILD',witness_sha256=saved['witness'],source_cache_sha256=before,N=64,all_nested_exact_fields_equal=True,elapsed_seconds=time.time()-tic,sha256_manifest=manifest)
    (OUT/'preparation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
