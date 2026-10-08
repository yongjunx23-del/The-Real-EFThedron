"""Exact aggregation and provenance checks for the frozen left180 witness."""
import hashlib,json
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/all_direction_certification_20261007/head/left'
SHA='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyadic(x):
    m,e=map(int,x);return F(m)*F(2)**e
source=json.loads((BASE/'moments128/manifest.json').read_text())
assert source['witness_sha256']==SHA
source_checks=[]
for t in source['tables']:
    label=t['label'];m=t['moment_metadata']
    assert 2*m['order']-1>=max(m['energy_degree_bound'],m['angular_degree_bound'])
    for suffix,key in [('_moments.json','moment_sha256'),('_uv_moments.json','uv_sha256')]:
        path=BASE/'moments128'/(label+suffix);assert sha(path)==t[key]
        again=BASE/'moments128_replay'/(label+suffix)
        replayed=again.exists() and sha(again)==t[key]
        source_checks.append(dict(file=str(path.relative_to(ROOT)),sha256=t[key],bit_for_bit_rebuild_match=replayed))
    assert t['uv_metadata']['source_table_sha256']==t['moment_sha256']
rows=[];missing=[];missing_replay=[];total=F(0);replayed_total=F(0)
for ell in range(0,41,2):
    folder=BASE/(f'low{ell}' if ell<4 else f'spin{ell}')
    p=folder/f'head_l{ell}.json';q=folder/f'replay_l{ell}.json'
    if not p.exists():missing.append(ell);continue
    d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['spin']==ell
    for name,digest in d['dependencies'].items():assert sha(ROOT/name)==digest,name
    assert sha(BASE/'moments128/manifest.json')==d['source_manifest_sha256']
    value=dyadic(d['upper_dyadic']);total+=value
    row=dict(spin=ell,path=str(p.relative_to(ROOT)),sha256=sha(p),upper_exact=str(value),
        upper_float=float(value),cells=d['cells'],seconds=d['elapsed_seconds'],target_met=d['target_met'])
    if q.exists():
        r=json.loads(q.read_text());assert r['passed'] and r['witness_sha256']==SHA and r['spin']==ell
        assert r['partition_covers_beta_zero_to_one_exactly'] and r['source_record_sha256']==sha(p)
        for name,digest in r['dependency_sha256'].items():assert sha(ROOT/name)==digest,name
        value_replay=dyadic(r['upper_dyadic']);replayed_total+=value_replay
        row.update(replay=str(q.relative_to(ROOT)),replay_sha256=sha(q),replayed_upper_exact=str(value_replay),
            replayed_upper_float=float(value_replay),replay_seconds=r['elapsed_seconds'])
    else:missing_replay.append(ell)
    rows.append(row)
complete=not missing;all_replayed=complete and not missing_replay and all(c['bit_for_bit_rebuild_match'] for c in source_checks)
if complete:assert total<F('7.75783')
conservative=sum((max(F(r['upper_exact']),F(r.get('replayed_upper_exact',r['upper_exact']))) for r in rows),F(0))
def upward(value,digits=6):
    scale=10**digits;n=-((-value.numerator*scale)//value.denominator)
    return f'{n//scale}.{n%scale:0{digits}d}'
deps=['build_moments.py','assemble_left180.py']
extra=['certified_roy_moments.py','certified_strip_moments.py','certified_uv_moments.py',
       'validated_quadrature.py','crossing_kernels.py','certify_real_tail.py']
code=[Path(__file__).with_name(x) for x in deps]+[ROOT/'research'/x for x in extra]
out=dict(status='COMPLETE_REPLAYED_FINITE_HEAD' if all_replayed else 'COMPLETE_RIGOROUS_FINITE_HEAD_REPLAY_PENDING' if complete else 'INCOMPLETE_FINITE_HEAD',
    witness_sha256=SHA,direction=[-1,0],included_spins=[r['spin'] for r in rows],
    missing_spins=missing,missing_replays=missing_replay,upper_exact=str(total),upper_float=float(total),
    replayed_partial_upper_exact=str(replayed_total),replayed_partial_upper_float=float(replayed_total),
    outward_decimal_ceiling=upward(conservative) if complete else None,records=rows,source_checks=source_checks,
    conservative_max_of_original_replay_upper_exact=str(conservative),
    conservative_max_of_original_replay_upper_float=float(conservative),
    assembly_and_source_code_sha256={str(p.relative_to(ROOT)):sha(p) for p in code},
    scope='Only listed finite spins at all energies. Full left support requires even0..40 complete and a separately certified same-witness spin42..infinity tail.',
    timing_note='Parallel per-component elapsed times are not total wall time.')
(BASE/'head_summary.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ['records','source_checks','assembly_and_source_code_sha256','upper_exact','replayed_partial_upper_exact']},indent=2))
