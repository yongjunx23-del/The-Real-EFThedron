"""Select the separately replayed refined spin0; leave base summary unchanged."""
import copy
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'output/all_direction_certification_20261007/head/left'
SHA = 'e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dyadic(pair):
    m,e = map(int,pair)
    return F(m)*F(2)**e


def upward(x,digits=6):
    scale=10**digits
    n=-((-x.numerator*scale)//x.denominator)
    return f'{n//scale}.{n%scale:0{digits}d}'


def main():
    original=BASE/'head_summary.json'
    out=copy.deepcopy(json.loads(original.read_text()))
    assert out['status']=='COMPLETE_REPLAYED_FINITE_HEAD'
    assert out['witness_sha256']==SHA and out['direction']==[-1,0]
    assert out['included_spins']==list(range(0,41,2))
    assert not out['missing_spins'] and not out['missing_replays']
    p=BASE/'low0_refined/head_l0.json'
    q=BASE/'low0_refined/replay_l0.json'
    d=json.loads(p.read_text());r=json.loads(q.read_text())
    assert d['witness_sha256']==r['witness_sha256']==SHA
    assert d['spin']==r['spin']==0 and d['target_met'] and r['passed']
    assert r['partition_covers_beta_zero_to_one_exactly']
    assert r['source_record_sha256']==sha(p)
    assert d['source_manifest_sha256']==sha(BASE/'moments128/manifest.json')
    for deps in [d['dependencies'],r['dependency_sha256']]:
        for name,digest in deps.items():
            assert sha(ROOT/name)==digest,name
    # Rebind every unchanged selected record/replay and its implementation.
    for row in out['records']:
        for key,hkey in [('path','sha256'),('replay','replay_sha256')]:
            record=ROOT/row[key]
            assert sha(record)==row[hkey]
            data=json.loads(record.read_text())
            assert data['witness_sha256']==SHA and data['spin']==row['spin']
            for name,digest in data.get('dependencies',data.get('dependency_sha256',{})).items():
                assert sha(ROOT/name)==digest,name
    previous=out['records'][0]
    value=dyadic(d['upper_dyadic']);replayed=dyadic(r['upper_dyadic'])
    assert value<F(previous['upper_exact'])
    out['records'][0]=dict(spin=0,path=str(p.relative_to(ROOT)),sha256=sha(p),upper_exact=str(value),
        upper_float=float(value),cells=d['cells'],seconds=d['elapsed_seconds'],target_met=d['target_met'],
        replay=str(q.relative_to(ROOT)),replay_sha256=sha(q),replayed_upper_exact=str(replayed),
        replayed_upper_float=float(replayed),replay_seconds=r['elapsed_seconds'])
    total=sum((F(row['upper_exact']) for row in out['records']),F(0))
    again=sum((F(row['replayed_upper_exact']) for row in out['records']),F(0))
    conservative=sum((max(F(row['upper_exact']),F(row['replayed_upper_exact'])) for row in out['records']),F(0))
    out.update(upper_exact=str(total),upper_float=float(total),replayed_partial_upper_exact=str(again),
        replayed_partial_upper_float=float(again),outward_decimal_ceiling=upward(conservative),
        conservative_max_of_original_replay_upper_exact=str(conservative),
        conservative_max_of_original_replay_upper_float=float(conservative),
        base_summary_path=str(original.relative_to(ROOT)),base_summary_sha256=sha(original),
        refined_spin=0,previous_spin0_record=previous)
    out['assembly_and_source_code_sha256'][str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    target=BASE/'head_refined_summary.json'
    target.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:out[k] for k in ['status','witness_sha256','upper_float','replayed_partial_upper_float',
        'conservative_max_of_original_replay_upper_float','outward_decimal_ceiling']},indent=2))


if __name__=='__main__':
    main()
