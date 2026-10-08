"""Select the tighter spin0 certificate without altering the original head."""
import hashlib,json
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1'
SHA='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyadic(x):
    m,e=map(int,x);return F(m)*F(2)**e
old=json.loads((BASE/'head_summary.json').read_text())
assert old['status']=='COMPLETE_RIGOROUS_FINITE_HEAD' and old['witness_sha256']==SHA
ref=BASE/'low0_refined/head_l0.json';r=json.loads(ref.read_text())
assert r['witness_sha256']==SHA and r['spin']==0 and r['target_met']
for name,digest in r['dependencies'].items():assert sha(ROOT/name)==digest
assert sha(BASE/'moments128/manifest.json')==r['source_manifest_sha256']
rows=[];total=F(0)
for entry in old['records']:
    p=ref if entry['spin']==0 else ROOT/entry['path'];d=json.loads(p.read_text())
    assert d['witness_sha256']==SHA
    for name,digest in d['dependencies'].items():assert sha(ROOT/name)==digest
    value=dyadic(d['upper_dyadic']);total+=value
    rows.append(dict(spin=d['spin'],path=str(p.relative_to(ROOT)),sha256=sha(p),
        upper_exact=str(value),upper_float=float(value),cells=d['cells'],seconds=d['elapsed_seconds']))
assert total<F('1.050950')
base_replay=json.loads((BASE/'final_head_audit.json').read_text())
rp=BASE/'low0_refined/replay_l0.json';complete=base_replay['status']=='COMPLETE_REPLAYED_FINITE_HEAD' and rp.exists()
replayed=None
if complete:
    rd=json.loads(rp.read_text())
    assert rd['passed'] and rd['partition_covers_beta_zero_to_one_exactly'] and rd['witness_sha256']==SHA
    assert rd['source_record_sha256']==sha(ref)
    for name,digest in rd['dependency_sha256'].items():assert sha(ROOT/name)==digest
    first=next(x for x in base_replay['records'] if x['spin']==0)
    replayed=F(base_replay['upper_exact'])-F(first['replayed_upper_exact'])+dyadic(rd['upper_dyadic'])
    assert replayed<F('1.050950')
out=dict(status='COMPLETE_REPLAYED_REFINED_FINITE_HEAD' if complete else 'COMPLETE_RIGOROUS_FINITE_HEAD_REPLAY_PENDING',
    witness_sha256=SHA,upper_exact=str(total),upper_float=float(total),outward_decimal_ceiling='1.050950',
    replay_upper_exact=str(replayed) if replayed is not None else None,
    replay_upper_float=float(replayed) if replayed is not None else None,records=rows,
    original_complete_head_sha256=sha(BASE/'head_summary.json'),
    original_replay_audit_sha256=sha(BASE/'final_head_audit.json'),
    refined_spin0_replay_sha256=sha(rp) if rp.exists() else None,
    source_rebuild_audit_sha256=sha(BASE/'source_rebuild_audit.json'),
    assembly_code_sha256=sha(__file__),
    scope='Even spins0..40 at all energies, same fixed robust-stage1 witness. The complete coupling bound requires the separately certified spins42..infinity tail.')
(BASE/'head_refined_summary.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ['records','upper_exact','replay_upper_exact']},indent=2))
