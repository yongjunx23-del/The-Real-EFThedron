"""Exact, hash-bound packets for independently replayed left objective shifts.

Only the finite all-energy head is assembled here.  The parent must attach a
complete same-base tail, using H(U,V-b*p)<=H(U,V), and normalize outwards.
Incomplete directions are explicitly marked and have no claimed full bound.
"""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'output/all_direction_certification_20261007/head/left'
SHA = 'e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
WITNESS = ROOT / 'output/full_direction_20261007/left/angle180_final.npz'
SLOPES = [F(1, 4), F(1, 2), F(1), F(2), F(4), F(8), F(16), F(32), F(64), F(128)]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dyadic(pair):
    m, e = map(int, pair)
    return F(m) * F(2) ** e


def relative(p):
    return str(Path(p).relative_to(ROOT))


def dependencies(mapping):
    for name, digest in mapping.items():
        assert sha(ROOT / name) == digest, name


def load_component(ell, folder):
    p = BASE / folder / f'head_l{ell}.json'
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    assert r['status'] == 'PASS_FIXED_PARTITION_FAMILY_HEAD_COMPONENT'
    assert r['base_witness_sha256'] == SHA and r['spin'] == ell
    assert r['source_manifest_sha256'] == sha(BASE / 'moments128/manifest.json')
    dependencies(r['dependencies'])
    slopes = list(map(F, r['slopes']))
    assert slopes == sorted(set(slopes)) and all(x > 0 for x in slopes)
    assert r['objectives'] == [['-1', str(-b)] for b in slopes]
    source = BASE / (f'low{ell}' if ell < 4 else f'spin{ell}') / f'head_l{ell}.json'
    assert sha(source) == r['source_record_sha256']
    s = json.loads(source.read_text())
    assert s['witness_sha256'] == SHA and s['spin'] == ell and s['target_met']
    dependencies(s['dependencies'])
    assert s['source_manifest_sha256'] == r['source_manifest_sha256']
    assert r['gap_intervals'] == s['gap_beta_intervals']
    cells = [(F(c['a']), F(c['b'])) for c in r['cells']]
    assert cells == [(F(c['a']), F(c['b'])) for c in s['cells_data']]
    partition = sorted(cells + [tuple(map(F, x)) for x in r['gap_intervals']])
    assert partition[0][0] == 0 and partition[-1][1] == 1
    assert all(a < b for a, b in partition)
    assert all(x[1] == y[0] for x, y in zip(partition, partition[1:]))
    replay_candidates = [
        source.with_name(f'replay_l{ell}.json'),
        BASE / 'family/base_replays' / f'l{ell}.json',
        BASE / 'family10/base_replays' / f'l{ell}.json',
    ]
    replay = next((q for q in replay_candidates if q.exists() and sha(q) == r['original_replay_sha256']), None)
    assert replay is not None, ('missing original replay', str(p))
    q = json.loads(replay.read_text())
    assert q['passed'] and q['witness_sha256'] == SHA and q['spin'] == ell
    assert q['partition_covers_beta_zero_to_one_exactly']
    assert q['source_record_sha256'] == sha(source)
    dependencies(q['dependency_sha256'])
    gap = dyadic(r['gap_bound_dyadic'])
    for k in range(len(slopes)):
        actual_sum = gap + sum((dyadic(c['shifted_upper_dyadic'][k]) for c in r['cells']), F(0))
        assert actual_sum <= dyadic(r['upper_dyadic'][k])
        for c in r['cells']:
            assert dyadic(c['shifted_upper_dyadic'][k]) <= dyadic(c['original_upper_dyadic'])
            if k:
                assert dyadic(c['shifted_upper_dyadic'][k]) <= dyadic(c['shifted_upper_dyadic'][k - 1])
    return dict(path=p, data=r, source=source, replay=replay, slopes=slopes)


def main():
    assert sha(WITNESS) == SHA
    source_manifest = BASE / 'moments128/manifest.json'
    source = json.loads(source_manifest.read_text())
    assert source['witness_sha256'] == SHA
    source_checks = []
    for t in source['tables']:
        m = t['moment_metadata']
        assert 2 * m['order'] - 1 >= max(m['energy_degree_bound'], m['angular_degree_bound'])
        for suffix, key in [('_moments.json', 'moment_sha256'), ('_uv_moments.json', 'uv_sha256')]:
            p = BASE / 'moments128' / (t['label'] + suffix)
            q = BASE / 'moments128_replay' / (t['label'] + suffix)
            assert sha(p) == sha(q) == t[key]
            source_checks.append(dict(path=relative(p), sha256=t[key], rebuilt_path=relative(q), exact_rebuild_match=True))
    components = {}
    for ell in range(0, 41, 2):
        # Each record remains independently available; ten slopes are preferred
        # when both exist, with no arithmetic mixing of different row banks.
        components[ell] = [x for x in (load_component(ell, 'family10'), load_component(ell, 'family')) if x]
    packets = []
    for b in SLOPES:
        selected = []
        total = F(0)
        deps = {}
        missing = []
        for ell in range(0, 41, 2):
            item = next((x for x in components[ell] if b in x['slopes']), None)
            if item is None:
                missing.append(ell)
                continue
            r = item['data']
            k = item['slopes'].index(b)
            value = dyadic(r['upper_dyadic'][k])
            total += value
            for name, digest in r['dependencies'].items():
                assert name not in deps or deps[name] == digest
                deps[name] = digest
            selected.append(dict(spin=ell, component_path=relative(item['path']), component_sha256=sha(item['path']),
                selected_slope_index=k, upper_exact=str(value), upper_float=float(value),
                source_record_path=relative(item['source']), source_record_sha256=sha(item['source']),
                original_replay_path=relative(item['replay']), original_replay_sha256=sha(item['replay']),
                cells=len(r['cells']), gap_intervals=r['gap_intervals'], precision_bits=r['precision_bits']))
        complete = not missing
        packets.append(dict(status='COMPLETE_REPLAYED_FINITE_HEAD' if complete else 'INCOMPLETE_FINITE_HEAD',
            slope=str(b), exact_objective=['-1', str(-b)], base_witness_path=relative(WITNESS), base_witness_sha256=SHA,
            included_spins=[x['spin'] for x in selected], missing_spins=missing,
            head_upper_exact=str(total) if complete else None, head_upper_float=float(total) if complete else None,
            partial_sum_exact=str(total), partial_sum_float=float(total), selected_components=selected,
            dependency_sha256=deps, source_manifest_path=relative(source_manifest), source_manifest_sha256=sha(source_manifest),
            scope='Finite even spins0..40, all energies; no infinite-spin tail or unit-direction normalization included.'))
    out = dict(status='COMPLETE_REPLAYED_FINITE_HEAD_FAMILY' if all(p['status'].startswith('COMPLETE') for p in packets) else 'INCOMPLETE_FAMILY',
        base_witness_sha256=SHA, witness_definition='Exact base multipliers and source rows; objective(-1,-b).',
        identity='U_b=U_base and V_b=V_base-b*p, with p>0 on the physical open domain.',
        tail_requirement='Attach complete same-base spin42..infinity bound via H(U,V-b*p)<=H(U,V); normalize outwards separately.',
        packets=packets, source_checks=source_checks, assembly_sha256=sha(__file__))
    output = BASE / 'family_packets.json'
    output.write_text(json.dumps(out, indent=2))
    print(json.dumps(dict(status=out['status'], packets=[{k:p[k] for k in ['slope','status','head_upper_float','missing_spins']} for p in packets]), indent=2))


if __name__ == '__main__':
    main()
