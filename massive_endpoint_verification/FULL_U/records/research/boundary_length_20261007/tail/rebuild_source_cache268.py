"""Read-only exact regeneration of every fixed-witness R1 source cache entry.

This is a deterministic exact-rational rebuild using the reviewed builder,
not an independent derivation of that builder's algebra.
"""
from pathlib import Path
import argparse,importlib,pickle,json,time
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--angle',required=True,choices=['268']);a=ap.parse_args();m=importlib.import_module('angle'+a.angle+'_campaign')
 import finite_source_corrections as f
 tic=time.time();manifest=m.source_manifest();p=m.OUT/'finite_source_N12_K80.pkl';old=pickle.loads(p.read_bytes());f.WITNESS=m.WITNESS;f.OUT=m.OUT;new=f.build(12,80)
 for k in old:
  if k!='seconds':assert old[k]==new[k],k
 assert old['witness_sha256']==m.SHA and len(old['moments'])==12 and all(len(r)==81 for r in old['moments'])
 files=[Path(__file__),p,Path(f.__file__)]
 manifest.update({str(q.relative_to(m.ROOT)):m.sha(q)for q in files})
 for q,h in manifest.items():assert m.sha(m.ROOT/q)==h,q
 out=dict(status='PASS_EXACT_RATIONAL_SOURCE_REBUILD',witness_sha256=m.SHA,objective_exact=list(map(str,m.OBJECTIVE)),moment_count=972,source_remainder_coefficients=sum(len(x['coefficient_abs_upper_exact'])for x in old['remainder_by_source']),all_entries_equal=True,independent_algebraic_derivation=False,cache_sha256=m.sha(p),sha256_manifest=manifest,seconds=time.time()-tic)
 target=m.OUT/'source_exact_rebuild.json';assert not target.exists();target.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
