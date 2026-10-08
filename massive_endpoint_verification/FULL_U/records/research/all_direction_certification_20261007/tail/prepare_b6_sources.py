"""New exact source data for the old150 witness; original records unchanged."""
from finite_source_corrections import *
import finite_source_corrections as src
from tail_certify_joint_limit import derive
W=ROOT/'output/angles_resummed/angle_000p000.npz'
FOLDER=ROOT/'output/all_direction_certification_20261007/tail/b6_extension'

def main():
    FOLDER.mkdir(parents=True,exist_ok=True)
    fixed=ROOT/'output/certification_tail/fixed_tail_certificate.json';rec=json.loads(fixed.read_text())
    for p,h in rec['sha256_manifest'].items():assert sha(ROOT/p)==h,p
    old=ROOT/'output/certification_tail/source_cache_v2_depth10.pkl';v=pickle.loads(old.read_bytes());assert v['sha']==sha(W)
    joint=dict(witness_sha256=sha(W),objective_exact=['1','0'],data=derive(W,80))
    (FOLDER/'joint_moments.pkl').write_bytes(pickle.dumps(joint))
    norms=dict(sha256=sha(W),norm1_rational=str(v['data']['norm1']),norm2_rational=str(v['data']['norm2']))
    (FOLDER/'source_norms.json').write_text(json.dumps(norms,indent=2))
    cp=ROOT/'output/certification_tail/compact_source_cache.pkl'
    (FOLDER/'compact_norms.pkl').write_bytes(cp.read_bytes())
    src.WITNESS=W;data=src.build(12,40);p=FOLDER/'finite_source_N12_K40.pkl';p.write_bytes(pickle.dumps(data))
    files=[Path(__file__),W,old,cp,fixed,Path(src.__file__),Path(derive.__code__.co_filename),FOLDER/'joint_moments.pkl',FOLDER/'source_norms.json',FOLDER/'compact_norms.pkl',p]
    hashes=dict(rec['sha256_manifest']);hashes.update({str(q.relative_to(ROOT)):sha(q)for q in files})
    out=dict(status='PASS_B6_EXACT_SOURCE_PREPARATION',witness_sha256=sha(W),objective_exact=['1','6'],delta_X4_exact='107/5',source_unchanged=True,K=40,N=12,sha256_manifest=hashes)
    (FOLDER/'prepared_components.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='sha256_manifest'},indent=2))

if __name__=='__main__':main()
