"""New left-witness source caches; all g2 and historical records stay immutable."""
from pathlib import Path
from fractions import Fraction as F
import json,pickle,hashlib,time
import finite_source_corrections as src
from compact_source_data import compute
ROOT=src.ROOT;OLD=ROOT/'output/all_direction_certification_20261007/tail/left';OUT=OLD/'refined'
WITNESS=ROOT/'output/full_direction_20261007/left/angle180_final.npz'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 tic=time.time();OUT.mkdir(exist_ok=True);base=json.loads((OLD/'prepared_components.json').read_text());assert base['witness_sha256']==sha(WITNESS)
 for p,h in base['sha256_manifest'].items():assert sha(ROOT/p)==h,p
 for name in['joint_moments.pkl','source_norms.json']:
  target=OUT/name
  if target.exists():assert target.read_bytes()==(OLD/name).read_bytes()
  else:target.write_bytes((OLD/name).read_bytes())
 src.WITNESS=WITNESS;src.OUT=OUT;data=src.build(12,80);p=OUT/'finite_source_N12_K80.pkl';p.write_bytes(pickle.dumps(data))
 report={k:v for k,v in data.items()if k not in['moments','remainder_coefficients']};report['remainder_coefficients_exact']=list(map(str,data['remainder_coefficients']))
 report['sha256_manifest']={str(q.relative_to(ROOT)):sha(q)for q in[WITNESS,p,Path(__file__),Path(src.__file__),ROOT/'research/tail_certification/certify_uv.py']}
 p.with_suffix('.json').write_text(json.dumps(report,indent=2))
 cd=compute(WITNESS,depth=10,N=40);(OUT/'compact_norms.pkl').write_bytes(pickle.dumps(cd))
 base['sha256_manifest'].update({str(q.relative_to(ROOT)):sha(q)for q in[Path(__file__),Path(compute.__code__.co_filename),OUT/'joint_moments.pkl',OUT/'source_norms.json',OUT/'compact_norms.pkl',p,p.with_suffix('.json')]})
 base['preparation_seconds']=time.time()-tic;base['scope']='Left witness copied verified base components, new exact finite-source corrections and compact signed moments.'
 (OUT/'prepared_components.json').write_text(json.dumps(base,indent=2));print(json.dumps(dict(status='PASS_LEFT_REFINED_COMPONENTS',seconds=time.time()-tic,witness_sha256=sha(WITNESS),source_cache_sha256=sha(p))))
if __name__=='__main__':main()
