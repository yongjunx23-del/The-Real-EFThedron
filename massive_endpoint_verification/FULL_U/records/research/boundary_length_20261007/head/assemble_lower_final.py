"""Hash-checked finite-head assembly for an arbitrary220-row saved witness."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/head'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyadic(pair):
    m,e=map(int,pair);return F(m)*F(2)**e
def upwards(x,n=9):
    scale=10**n;k=-((-x.numerator*scale)//x.denominator)
    return f'{k//scale}.{k%scale:0{n}d}'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True)
    ap.add_argument('--base',type=Path,required=True);a=ap.parse_args()
    base=a.base.resolve();wp=a.witness.resolve();wsha=sha(wp)
    direction=[str(F(float(x)))for x in np.load(wp)['direction']]
    manifest=base/'moments128/manifest.json';m=json.loads(manifest.read_text());assert m['witness_sha256']==wsha
    checks=[];source_rebuilt=True
    for t in m['tables']:
        meta=t['moment_metadata'];uv=t['uv_metadata']
        assert meta['witness_sha256']==uv['witness_sha256']==wsha
        assert 2*meta['order']-1>=max(meta['energy_degree_bound'],meta['angular_degree_bound'])
        assert meta['jmax']>=74 and meta['kmax']>=128 and uv['N']>=32 and uv['jmax']>=40
        assert uv['source_table_sha256']==t['moment_sha256']
        for suffix,key in [('_moments.json','moment_sha256'),('_uv_moments.json','uv_sha256')]:
            p=base/'moments128'/(t['label']+suffix);q=base/'moments128_replay'/(t['label']+suffix)
            assert sha(p)==t[key]
            rebuilt=q.exists()and sha(q)==t[key];source_rebuilt &= rebuilt
            checks.append(dict(file=str(p.relative_to(ROOT)),sha256=t[key],rebuilt_path=str(q.relative_to(ROOT)),bit_for_bit_rebuild_match=rebuilt))
    rows=[];missing=[];pending=[];replay_missing=[];total=F(0);again=F(0);conservative=F(0)
    deps={str(Path(__file__).relative_to(ROOT)):sha(__file__),str((OLD/'build_moments.py').relative_to(ROOT)):sha(OLD/'build_moments.py')}
    review=Path(__file__).with_name('REVIEW_CERTIFY_HEAD.md')
    assert review.exists();deps[str(review.relative_to(ROOT))]=sha(review)
    for ell in range(0,41,2):
        selected_dir=base/(f'spin{ell}_constant' if base.name=='angle268' and ell==2 else f'spin{ell}')
        p=selected_dir/f'head_l{ell}.json';q=selected_dir/f'replay_l{ell}.json'
        if not p.exists():missing.append(ell);continue
        d=json.loads(p.read_text())
        if not d['final']:pending.append(ell);continue
        assert d['status']=='PASS_FINITE_SPIN_ENCLOSURE'and d['witness_sha256']==wsha and d['spin']==ell
        assert d['direction_exact']==direction and d['source_manifest_sha256']==sha(manifest)
        for name,h in d['dependencies'].items():assert sha(ROOT/name)==h;deps[name]=h
        gaps=[tuple(map(F,x))for x in d['gap_beta_intervals']]
        cells=[(F(c['a']),F(c['b']))for c in d['cells_data']];union=sorted(cells+gaps)
        assert len(cells)==d['cells']and union[0][0]==0 and union[-1][1]==1
        assert all(lo<hi for lo,hi in union)and all(x[1]==y[0]for x,y in zip(union,union[1:]))
        value=dyadic(d['upper_dyadic'])
        assert dyadic(d['gap_upper_dyadic'])+sum((dyadic(c['upper_dyadic'])for c in d['cells_data']),F(0))<=value
        total+=value
        row=dict(spin=ell,path=str(p.relative_to(ROOT)),sha256=sha(p),upper_exact=str(value),upper_float=float(value),
            cells=d['cells'],seconds=d['elapsed_seconds'],target_met=d['target_met'],model=d['model'])
        selected=value
        if q.exists():
            r=json.loads(q.read_text())
            assert r['status']=='PASS_FULL_PARTITION_REPLAY'and r['passed']and r['witness_sha256']==wsha and r['spin']==ell
            assert r['direction_exact']==direction and r['source_record_sha256']==sha(p)and r['source_manifest_sha256']==sha(manifest)
            assert r['precision_bits']>d['precision_bits'],'Replay must use strictly higher precision.'
            assert r['partition_covers_beta_zero_to_one_exactly']and r['cells_replayed']==len(cells)and r['gaps_replayed']==len(gaps)
            for name,h in r['dependency_sha256'].items():assert sha(ROOT/name)==h;deps[name]=h
            rv=dyadic(r['upper_dyadic']);again+=rv;selected=max(selected,rv)
            row.update(replay=str(q.relative_to(ROOT)),replay_sha256=sha(q),replayed_upper_exact=str(rv),replayed_upper_float=float(rv),replay_seconds=r['elapsed_seconds'])
        else:replay_missing.append(ell)
        conservative+=selected;rows.append(row)
    complete=not missing and not pending;replayed=complete and not replay_missing and source_rebuilt
    out=dict(status='COMPLETE_REPLAYED_FINITE_HEAD'if replayed else'COMPLETE_RIGOROUS_FINITE_HEAD_REPLAY_PENDING'if complete else'INCOMPLETE_FINITE_HEAD',
        witness_sha256=wsha,witness_path=str(wp.relative_to(ROOT)),direction_exact=direction,
        included_spins=[r['spin']for r in rows],missing_spins=missing,unfinished_spins=pending,missing_replays=replay_missing,
        upper_exact=str(total),upper_float=float(total),replayed_partial_upper_exact=str(again),replayed_partial_upper_float=float(again),
        conservative_max_of_original_replay_upper_exact=str(conservative),conservative_max_of_original_replay_upper_float=float(conservative),
        outward_decimal_ceiling=upwards(conservative)if complete else None,records=rows,source_checks=checks,
        source_manifest_path=str(manifest.relative_to(ROOT)),source_manifest_sha256=sha(manifest),assembly_and_source_code_sha256=deps,
        scope='Only listed finite even spins0..40, every energy. A complete same-witness spin42..infinity tail is separate.',
        timing_note='Parallel per-job elapsed times are not total wall time; numerical targets are quality budgets only.')
    (base/('head_summary_final.json' if base.name=='angle268' else 'head_summary.json')).write_text(json.dumps(out,indent=2))
    print(json.dumps({k:out[k]for k in ['status','included_spins','missing_spins','unfinished_spins','missing_replays','conservative_max_of_original_replay_upper_float','outward_decimal_ceiling']},indent=2))


if __name__=='__main__':main()
