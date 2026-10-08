"""Compare freshly rebuilt exact source tables with the certificate inputs.

This checks content hashes and the complete polynomial Gauss degree
conditions. It does not substitute for the analytic remainder proofs.
"""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1'
SHA='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'
old=BASE/'moments128';new=BASE/'moments128_replay'
om=json.loads((old/'manifest.json').read_text());nm=json.loads((new/'manifest.json').read_text())
assert om['witness_sha256']==nm['witness_sha256']==SHA
checks=[]
for a,b in zip(om['tables'],nm['tables']):
    assert a['label']==b['label']
    label=a['label']
    for suffix,key in [('_moments.json','moment_sha256'),('_uv_moments.json','uv_sha256')]:
        x=hashlib.sha256((old/(label+suffix)).read_bytes()).hexdigest()
        y=hashlib.sha256((new/(label+suffix)).read_bytes()).hexdigest()
        assert x==y==a[key]==b[key],label+suffix
        checks.append(dict(file=label+suffix,sha256=x,bit_for_bit_rebuild_match=True))
    m=a['moment_metadata'];assert m['witness_sha256']==SHA
    assert 2*m['order']-1>=max(m['energy_degree_bound'],m['angular_degree_bound'])
    assert a['uv_metadata']['source_table_sha256']==a['moment_sha256']
deps=['research/all_direction_certification_20261007/head/build_moments.py',
      'research/certified_roy_moments.py','research/certified_strip_moments.py',
      'research/certified_uv_moments.py','research/validated_quadrature.py',
      'research/all_direction_certification_20261007/head/audit_stage1_sources.py']
out=dict(passed=True,witness_sha256=SHA,checks=checks,
    original_manifest_sha256=hashlib.sha256((old/'manifest.json').read_bytes()).hexdigest(),
    rebuilt_manifest_sha256=hashlib.sha256((new/'manifest.json').read_bytes()).hexdigest(),
    original_seconds=om['seconds'],rebuild_seconds=nm['seconds'],
    source_code_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps},
    scope='All six polynomial source/UV tables rebuilt and matched exactly; analytic remainder lemmas remain separate.')
(BASE/'source_rebuild_audit.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
