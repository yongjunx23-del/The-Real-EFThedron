"""Exact dyadic sum of completed robust-stage1 finite-head certificates."""
import hashlib,json
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1'
SHA='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'
def folder(ell):
    if ell in (0,2):return 'low'
    if ell in (4,6):return 'quartic_low4'
    if ell>=14:return f'gap_l{ell}'
    return {8:'quartic_high',10:'quartic_10_20',22:'quartic_22_30',32:'quartic_32_40'}.get(ell,f'single_l{ell}')
rows=[];missing=[];total=F(0)
for ell in range(0,41,2):
    path=BASE/folder(ell)/f'head_l{ell}.json'
    if not path.exists():missing.append(ell);continue
    d=json.loads(path.read_text());assert d['witness_sha256']==SHA and d['spin']==ell
    for p,sha in d['dependencies'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha,p
    assert hashlib.sha256((BASE/'moments128/manifest.json').read_bytes()).hexdigest()==d['source_manifest_sha256']
    m,e=map(int,d['upper_dyadic']);value=F(m)*F(2)**e;total+=value
    rows.append(dict(spin=ell,path=str(path.relative_to(ROOT)),record_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        upper_exact=str(value),upper_float=float(value),target_met=d['target_met'],cells=d['cells'],seconds=d['elapsed_seconds']))
out=dict(status='COMPLETE_RIGOROUS_FINITE_HEAD' if not missing else 'INCOMPLETE_FINITE_HEAD',witness_sha256=SHA,
    included_spins=[r['spin'] for r in rows],missing_spins=missing,upper_exact=str(total),upper_float=float(total),records=rows,
    claim='All-energy head for even spins0..40 only when complete; higher-spin tail remains separate',
    replay_status='Per-record replays tracked separately; this script checks input/code hashes and sums their directed dyadic bounds exactly')
BASE.mkdir(parents=True,exist_ok=True);(BASE/'head_summary.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ('records','upper_exact')},indent=2))
