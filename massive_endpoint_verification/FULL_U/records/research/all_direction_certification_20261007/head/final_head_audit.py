"""Assemble complete same-witness head records and full-partition replays."""
import hashlib,json,sys
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps')]
import flint,numpy,sympy
BASE=ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1'
SHA='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyadic(x):
    m,e=map(int,x);return F(m)*F(2)**e
summary=json.loads((BASE/'head_summary.json').read_text())
assert summary['witness_sha256']==SHA
sources=json.loads((BASE/'source_rebuild_audit.json').read_text())
assert sources['passed'] and sources['witness_sha256']==SHA
for name,digest in sources['source_code_sha256'].items():assert sha(ROOT/name)==digest,name
for c in sources['checks']:
    assert sha(BASE/'moments128'/c['file'])==sha(BASE/'moments128_replay'/c['file'])==c['sha256']
rows=[];missing=[];total=F(0);original=F(0)
byspin={r['spin']:r for r in summary['records']}
for ell in range(0,41,2):
    if ell not in byspin:missing.append(dict(spin=ell,reason='No completed adaptive certificate'));continue
    old=byspin[ell];p=ROOT/old['path'];d=json.loads(p.read_text())
    assert sha(p)==old['record_sha256'] and d['witness_sha256']==SHA
    for name,digest in d['dependencies'].items():assert sha(ROOT/name)==digest,name
    assert sha(BASE/'moments128/manifest.json')==d['source_manifest_sha256']
    original+=dyadic(d['upper_dyadic'])
    q=BASE/'replays'/f'l{ell}.json'
    if not q.exists():missing.append(dict(spin=ell,reason='No completed partition replay'));continue
    replay=json.loads(q.read_text())
    assert replay['passed'] and replay['witness_sha256']==SHA
    assert replay['source_record_sha256']==sha(p)
    assert replay['partition_covers_beta_zero_to_one_exactly']
    for name,digest in replay['dependency_sha256'].items():assert sha(ROOT/name)==digest,name
    value=dyadic(replay['upper_dyadic']);total+=value
    rows.append(dict(spin=ell,record=old['path'],record_sha256=sha(p),
        replay=str(q.relative_to(ROOT)),replay_sha256=sha(q),
        replayed_upper_exact=str(value),replayed_upper_float=float(value),
        cells=replay['cells_replayed'],gaps=replay['gaps_replayed'],
        adaptive_seconds=d['elapsed_seconds'],replay_seconds=replay['elapsed_seconds']))
code=list((ROOT/'research/all_direction_certification_20261007/head').glob('*.py'))
code += [ROOT/'research'/n for n in ['crossing_kernels.py','certify_real_tail.py',
    'head_certification_lowspin.py','head_certification_spin0.py','validated_quadrature.py',
    'certified_roy_moments.py','certified_strip_moments.py','certified_uv_moments.py']]
complete=not missing
if complete:assert total<F('1.058016')
out=dict(status='COMPLETE_REPLAYED_FINITE_HEAD' if complete else 'INCOMPLETE_FINITE_HEAD_REPLAY',
    witness_sha256=SHA,spins_included=[r['spin'] for r in rows],missing=missing,
    upper_exact=str(total),upper_float=float(total),original_upper_exact=str(original),
    conservative_decimal_ceiling='1.058016' if complete else None,
    records=rows,source_rebuild_audit_sha256=sha(BASE/'source_rebuild_audit.json'),
    all_head_code_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(code))},
    versions=dict(python=sys.version,python_flint=flint.__version__,numpy=numpy.__version__,sympy=sympy.__version__),
    scope='Even spins0..40, every energy4<mu<infinity for the fixed witness. A full coupling bound additionally needs same-witness spins42..infinity.',
    timing_note='Component elapsed times are recorded individually. Parallel-job elapsed times must not be summed as total wall time.')
(BASE/'final_head_audit.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ['records','all_head_code_sha256','upper_exact','original_upper_exact']},indent=2))
