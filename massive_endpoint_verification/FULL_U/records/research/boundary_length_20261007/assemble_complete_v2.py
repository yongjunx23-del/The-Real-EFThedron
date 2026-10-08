"""Fail-closed assembly of a same-witness all-spin, all-energy certificate.

This does not optimize or sample the functional. Every input must already
contain a completed interval proof and its higher-precision replay.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math, sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'research'))
from angle_benchmarks import support_report


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dyadic(pair):
    m, e = map(int, pair)
    return F(m) * F(2) ** e


def ceiling_decimal(v, places=6):
    scale = 10 ** places
    n = -((-v.numerator * scale) // v.denominator)
    return f'{n // scale}.{n % scale:0{places}d}'


class Inputs:
    def __init__(self):
        self.manifest = {}
        self.seen = set()

    def bind(self, path, expected=None):
        p = Path(path)
        p = (p if p.is_absolute() else ROOT / p).resolve()
        key = str(p.relative_to(ROOT))
        actual = sha(p)
        if expected is not None:
            assert actual == expected, key
        assert key not in self.manifest or self.manifest[key] == actual, key
        self.manifest[key] = actual
        return p

    def read(self, path, expected=None):
        p = self.bind(path, expected)
        d = json.loads(p.read_text())
        if p not in self.seen:
            self.seen.add(p)
            for sub, h in d.get('sha256_manifest', {}).items():
                q = self.bind(sub, h)
                if q.suffix == '.json':
                    self.read(q, h)
        return d


def head_cost(inputs, path, digest, direction):
    d = inputs.read(path)
    assert d['status'] == 'COMPLETE_REPLAYED_FINITE_HEAD'
    assert d['witness_sha256'] == digest and d['direction_exact'] == direction
    assert d['included_spins'] == list(range(0, 41, 2))
    assert not d['missing_spins'] and not d['unfinished_spins'] and not d['missing_replays']
    for p, h in d['assembly_and_source_code_sha256'].items():
        inputs.bind(p, h)
    for row in d['source_checks']:
        assert row['bit_for_bit_rebuild_match']
        inputs.bind(row['file'], row['sha256'])
        inputs.bind(row['rebuilt_path'], row['sha256'])
    sources = inputs.read(d['source_manifest_path'], d['source_manifest_sha256'])
    assert sources['witness_sha256'] == digest
    for table in sources['tables']:
        m = table['moment_metadata']
        assert 2*m['order']-1 >= max(m['energy_degree_bound'], m['angular_degree_bound'])
    first = replay = retained = F(0)
    assert [r['spin'] for r in d['records']] == list(range(0, 41, 2))
    for row in d['records']:
        p = inputs.read(row['path'], row['sha256'])
        q = inputs.read(row['replay'], row['replay_sha256'])
        assert p['status'] == 'PASS_FINITE_SPIN_ENCLOSURE' and p['final']
        assert q['status'] == 'PASS_FULL_PARTITION_REPLAY' and q['passed']
        assert p['witness_sha256'] == q['witness_sha256'] == digest
        assert p['direction_exact'] == q['direction_exact'] == direction
        assert p['spin'] == q['spin'] == row['spin']
        assert p['source_manifest_sha256'] == q['source_manifest_sha256'] == d['source_manifest_sha256']
        assert q['precision_bits'] > p['precision_bits']
        assert q['source_record_sha256'] == row['sha256']
        assert q['partition_covers_beta_zero_to_one_exactly']
        cells = [(F(c['a']), F(c['b'])) for c in p['cells_data']]
        gaps = [tuple(map(F, c)) for c in p['gap_beta_intervals']]
        union = sorted(cells + gaps)
        assert union[0][0] == 0 and union[-1][1] == 1
        assert all(a < b for a, b in union)
        assert all(a[1] == b[0] for a, b in zip(union, union[1:]))
        assert q['cells_replayed'] == len(cells) and q['gaps_replayed'] == len(gaps)
        for group in (p['dependencies'], q['dependency_sha256']):
            for file, h in group.items():
                inputs.bind(file, h)
        x, y = dyadic(p['upper_dyadic']), dyadic(q['upper_dyadic'])
        assert x >= 0 and y >= 0
        if 'combined_chunk_paths' in q:
            chunks = [inputs.read(path) for path in q['combined_chunk_paths']]
            chunks.sort(key=lambda r:r['first_cell'])
            cursor = 0
            chunk_sum = F(0)
            for chunk in chunks:
                assert chunk['status'] == 'PASS_DISJOINT_REPLAY_CHUNK'
                for key in ['witness_sha256', 'source_manifest_sha256', 'source_record_sha256', 'spin', 'model', 'direction_exact', 'precision_bits']:
                    assert chunk[key] == q[key]
                assert chunk['first_cell'] == cursor < chunk['last_cell_exclusive'] <= len(cells)
                assert chunk['cells_replayed'] == chunk['last_cell_exclusive'] - cursor
                assert chunk['gaps_replayed'] == (len(gaps) if cursor == 0 else 0)
                assert chunk['gap_beta_intervals'] == (p['gap_beta_intervals'] if cursor == 0 else [])
                for file, h in chunk['dependency_sha256'].items():
                    inputs.bind(file, h)
                chunk_sum += dyadic(chunk['upper_dyadic'])
                cursor = chunk['last_cell_exclusive']
            assert cursor == len(cells) and chunk_sum == y
        assert dyadic(p['gap_upper_dyadic']) + sum((dyadic(c['upper_dyadic']) for c in p['cells_data']), F(0)) <= x
        assert x == F(row['upper_exact']) and y == F(row['replayed_upper_exact'])
        first += x
        replay += y
        retained += max(x, y)
    assert first == F(d['upper_exact']) and replay == F(d['replayed_partial_upper_exact'])
    assert retained == F(d['conservative_max_of_original_replay_upper_exact'])
    return retained


def real_norms(inputs, path, digest, mode):
    d = inputs.read(path)
    assert d['status'] == 'PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY'
    assert d['witness_sha256'] == digest and d['mode'] == mode
    assert (d['first_even_spin'], d['last_even_spin']) == ((42, 448) if mode == 'finite' else (450, 'infinity'))
    assert d['replay_precision_bits'] > d['precision_bits']
    cut = F(d['compact_beta_cutoff_exact'])
    assert 0 < cut < 1 and cut*cut <= F(45, 49)
    assert F(d['compact_mu_cutoff_exact']) == 4/(1-cut*cut)
    sums = {'compact': F(0), 'remainder': F(0)}
    first = {k: F(0) for k in sums}
    replay = {k: F(0) for k in sums}
    cursor = F(0)
    for cell in d['cells']:
        lo, hi = F(cell['lo']), F(cell['hi'])
        assert cursor == lo < hi <= 1
        cursor = hi
        region = cell['region']
        assert region in sums
        assert hi <= cut if region == 'compact' else lo >= cut
        x, y, z = map(F, [cell['first_upper_exact'], cell['replay_upper_exact'], cell['retained_upper_exact']])
        assert x >= 0 and y >= 0 and z == max(x, y)
        first[region] += x
        replay[region] += y
        sums[region] += z
    assert cursor == 1
    for k in sums:
        assert sums[k] == F(d['A_upper_exact'][k])
        assert first[k] == F(d['first_upper_exact'][k])
        assert replay[k] == F(d['replay_upper_exact'][k])
    return sums


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--angle', choices=['85', '255', '259', '265'], required=True)
    ap.add_argument('--head', type=Path, help='Optional frozen same-witness head summary, e.g. a parallel replay assembly.')
    ap.add_argument('--out', type=Path)
    a = ap.parse_args()
    base = ROOT / 'output/boundary_length_20261007'
    name = 'angle' + a.angle
    witness = ROOT / {
        '85': 'output/full_direction_20261007/right/angle85_220.npz',
        '255': 'output/full_direction_20261007/left/angle255.npz',
        '265': 'output/boundary_length_20261007/search/angle265.npz',
        '259': 'output/boundary_length_20261007/search/angle259p2941_from255.npz',
    }[a.angle]
    digest = sha(witness)
    inputs = Inputs()
    inputs.bind(witness)
    with np.load(witness) as data:
        direction = [str(F(float(x))) for x in data['direction']]
        assert len(data['names']) == len(data['multipliers']) == 220
    head = head_cost(inputs, a.head or base/'head'/name/'head_summary.json', digest, direction)
    if a.angle in ['259', '265']:
        cap = inputs.read(base/'head'/name/'constant_cap_preflight.json')
        assert cap['status'] == 'PASS_CONSTANT_CAP_PREFLIGHT' and cap['witness_sha256'] == digest and cap['spin'] == 2
        assert [r['bits'] for r in cap['rows']] == [256, 384]
        inputs.bind(ROOT/('research/boundary_length_20261007/head/REVIEW_CONSTANT_CAP_259.md' if a.angle == '259' else 'research/boundary_length_20261007/head/REVIEW_CONSTANT_CAP.md'))
        if a.angle == '259':
            endpoints = inputs.read(base/'head'/name/'constant_cap_exact_sign.json')
            assert endpoints['status'] == 'PASS_EXACT_NEGATIVE_CONSTANT_CAP_ENDPOINTS'
            assert endpoints['witness_sha256'] == digest and endpoints['spin'] == 2
            assert endpoints['source_manifest_sha256'] == cap['source_manifest_sha256']
            for file, h in endpoints['dependency_sha256'].items():
                inputs.bind(file, h)
            assert [r['bits'] for r in endpoints['records']] == [256, 384]
            for row in endpoints['records']:
                assert row['beta_domain_exact'] == ['0', '1/256']
                assert dyadic(row['V_over_p']['lower_dyadic']) <= dyadic(row['V_over_p']['upper_dyadic']) < 0
    compact = inputs.read(base/'compact'/name/'complete_compact_coercivity.json')
    assert compact['status'] == 'PASS_COMPLETE_COMPACT_COERCIVITY'
    assert compact['witness_sha256'] == digest and compact['objective_exact'] == direction
    assert compact['minimum_even_spin'] == 42 and compact['mu_domain_exact'] == ['4', '49']
    assert compact['lower_energy_endpoint_open']
    ck = F(compact['kappa_exact'])
    uv = inputs.read(base/'tail'/name/'uv_complete.json')
    assert uv['status'] == 'PASS_COMPLETE_DIRECTIONAL_UV_COERCIVITY'
    assert uv['witness_sha256'] == digest and uv['objective_exact'] == direction
    assert (uv['finite_first_even_spin'], uv['finite_last_even_spin'], uv['infinite_first_even_spin']) == (42, 448, 450)
    assert uv['energy_domain'] == ['49', 'infinity']
    fk, ik = F(uv['finite_kappa_exact']), F(uv['infinite_kappa_exact'])
    assert min(ck, fk, ik) > 0
    algebra = inputs.read(base/'tail'/name/'exact_rebuild_and_algebra_audit.json')
    assert algebra['status'] == 'PASS_UV_EXACT_REBUILD_AND_ALGEBRA'
    assert algebra['witness_sha256'] == digest and algebra['objective_exact'] == direction
    assert algebra['all_entries_equal'] and algebra['moment_count'] == 972 and algebra['remainder_coefficient_count'] == 9
    assert algebra['source_polynomial_equality_spins'] == algebra['crossing_polynomial_equality_spins'] == [450, 452, 600]
    assert algebra['precision_bits'] == 512 and len(algebra['checks']) == 18
    finite = real_norms(inputs, base/'real'/name/'finite_split_replay.json', digest, 'finite')
    infinite = real_norms(inputs, base/'real'/name/'tail_split_replay.json', digest, 'tail')
    # The exact dyadic beta cutoff is rounded below sqrt(45/49). The tiny
    # compact-energy sliver belongs to "remainder", so use min(ck, UV kappa).
    components = {
        'head_0_40': head,
        'finite_42_448_compact': finite['compact']/(2*ck),
        'finite_42_448_UV_and_sliver': finite['remainder']/(2*min(ck, fk)),
        'infinite_450_up_compact': infinite['compact']/(2*ck),
        'infinite_450_up_UV_and_sliver': infinite['remainder']/(2*min(ck, ik)),
    }
    check = inputs.read(base/'real'/name/'independent_field_check.json')
    assert check['status'] == 'PASS_INDEPENDENT_REAL_KERNEL_COMPARISON' and check['witness_sha256'] == digest
    for p in [Path(__file__), ROOT/'research/all_direction_certification_20261007/check_real220_fields.py', ROOT/'research/boundary_length_20261007/REVIEW_NEW_ENVELOPES.md', ROOT/'research/boundary_length_20261007/head/REVIEW_COMPLETE_ASSEMBLY.md', ROOT/'research/all_direction_certification_20261007/head/REVIEW_REAL220.md', ROOT/'research/all_direction_certification_20261007/alternative/rational_head.py', ROOT/'research/all_direction_certification_20261007/alternative/exact_bernstein.py', ROOT/'research/lower_g0_upper_g2_20261006/matrix_real_tail_gram.py']:
        inputs.bind(p)
    total = sum(components.values(), F(0))
    outward = ceiling_decimal(total)
    assert F(outward) >= total
    dx, dy = map(F, direction)
    angle = math.degrees(math.atan2(float(dy), float(dx))) % 360
    norm = math.hypot(float(dx), float(dy))
    dual, primal = support_report(angle), support_report(angle, 'primal')
    result = dict(
        status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND', witness=str(witness.relative_to(ROOT)),
        witness_sha256=digest, objective_exact=direction, nominal_angle_degrees=angle,
        direction_norm_squared_exact=str(dx*dx+dy*dy),
        inequality=f'({direction[0]}) g0 + ({direction[1]}) g2 <= {outward}',
        support_upper_exact=str(total), support_upper_outward_6dp=outward,
        support_normalized_float=float(F(outward))/norm,
        components_exact={k:str(v) for k,v in components.items()},
        components_float={k:float(v) for k,v in components.items()},
        coercivity_constants_exact=dict(compact=str(ck), finite_UV=str(fk), infinite_UV=str(ik)),
        tail_upper_exact=str(total-head), tail_upper_float=float(total-head),
        coverage='Every even spin ell>=0 and every energy mu>4, including threshold and source endpoint caps, infinite energy and infinite spin.',
        graphical_comparison=dict(dual=dual, primal=primal, relative_to_dual_estimate_percent=100*(float(F(outward))/norm/dual['support']-1)),
        comparison_warning='The published curves were extracted graphically. They are not exact author optimizer data, and analytic assumptions must match for a scientific dominance claim.',
        proof_scope='Internal computer-assisted directed interval certificate, exact stored witness, source identities and remainders, rational partitions and higher-precision replays. Not formal verification, external review, or an optimality claim.',
        sha256_manifest=inputs.manifest,
    )
    target = a.out or base/f'complete_{name}_certificate.json'
    assert not target.exists(), 'Choose a new output path; do not overwrite a frozen completed certificate.'
    target.write_text(json.dumps(result, indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['sha256_manifest','components_exact','support_upper_exact','tail_upper_exact','graphical_comparison']}, indent=2))


if __name__ == '__main__':
    main()

