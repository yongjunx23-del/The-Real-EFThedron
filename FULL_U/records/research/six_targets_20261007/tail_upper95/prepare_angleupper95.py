"""New source caches with the frozen angleupper95 binary64 coefficients as rationals."""
from angleupper95_campaign import *
import time,pickle
from tail_certify_joint_limit import derive
from matrix_source_norms import source_norms
from compact_source_data import compute
import finite_source_corrections as finite

def main():
 tic=time.time();OUT.mkdir(parents=True,exist_ok=True)
 for name in['roy_exact.json','crossing_exact.json']:
  v=json.loads((AUDITS/name).read_text());assert v['sha256']==SHA and v['passed']
 cache=OUT/'joint_moments.pkl'
 if not cache.exists():cache.write_bytes(pickle.dumps(dict(witness_sha256=SHA,objective_exact=list(map(str,OBJECTIVE)),data=derive(WITNESS))))
 assert pickle.loads(cache.read_bytes())['witness_sha256']==SHA
 norm=OUT/'source_norms.json'
 if not norm.exists():norm.write_text(json.dumps(source_norms(WITNESS,16,16),indent=2))
 assert json.loads(norm.read_text())['sha256']==SHA
 compact=OUT/'compact_norms.pkl'
 if not compact.exists():compact.write_bytes(pickle.dumps(compute(WITNESS,depth=10,N=40)))
 assert pickle.loads(compact.read_bytes())['witness_sha256']==SHA
 source=OUT/'finite_source_N12_K80.pkl'
 if not source.exists():
  finite.WITNESS=WITNESS;finite.OUT=OUT;source.write_bytes(pickle.dumps(finite.build(12,80)))
 assert pickle.loads(source.read_bytes())['witness_sha256']==SHA
 files=[WITNESS,Path(__file__),Path(__file__).with_name('angleupper95_campaign.py'),Path(derive.__code__.co_filename),Path(source_norms.__code__.co_filename),Path(compute.__code__.co_filename),Path(finite.__file__),cache,norm,compact,source,AUDITS/'roy_exact.json',AUDITS/'crossing_exact.json',ROOT/'research/lower_g0_upper_g2_20261006/matrix_roy_guards.py',ROOT/'research/tail_certification/certify_uv.py']
 out=dict(status='PASS_ANGLEupper95_PREPARED_COMPONENTS',witness_sha256=SHA,objective_exact=list(map(str,OBJECTIVE)),norm1_approx=json.loads(norm.read_text())['norm1_approx'],norm2_approx=json.loads(norm.read_text())['norm2_approx'],seconds=time.time()-tic,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in files})
 (OUT/'prepared_components.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='sha256_manifest'},indent=2))

if __name__=='__main__':main()
