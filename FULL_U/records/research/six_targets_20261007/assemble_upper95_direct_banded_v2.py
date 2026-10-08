"""Replace the complete global packet's finite-spin tail by its reviewed local bound.

This is an exact same-witness component substitution, with no kernel evaluation.
The original complete packet and all of its proof bytes remain unchanged.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, math
from assemble_upper95_direct_v2 import Inputs, ROOT, BASE, SHA, sha, ceiling_decimal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--global-certificate', type=Path, default=ROOT/'output/six_targets_20261007/complete_upper95_direct_global_v2_certificate.json')
    ap.add_argument('--out', type=Path, default=ROOT/'output/six_targets_20261007/complete_upper95_direct_certificate.json')
    args = ap.parse_args()
    assert not args.out.exists()
    inputs = Inputs()
    old = inputs.read(args.global_certificate)
    assert old['status'] == 'PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
    assert old['witness_sha256'] == SHA
    direction = old['objective_exact']
    parts = {k: F(v) for k, v in old['components_exact'].items()}
    expected = {'head_0_40', 'finite_42_448_compact', 'finite_42_448_UV_and_sliver',
                'infinite_450_up_compact', 'infinite_450_up_UV_and_sliver'}
    assert set(parts) == expected and all(v >= 0 for v in parts.values())
    assert sum(parts.values(), F(0)) == F(old['support_upper_exact'])
    assert sum(v for k, v in parts.items() if k != 'head_0_40') == F(old['tail_upper_exact'])

    band_path = BASE/'real/finite_banded_cost.json'
    band = inputs.read(band_path)
    audit = inputs.read(BASE/'tail/independent_coverage_accounting.json')
    real = inputs.read(BASE/'real/finite_split_replay.json')
    assert band['status'] == 'PASS_BANDED_COERCIVITY_REAL_COST'
    assert audit['status'] == 'PASS_UPPER_TAIL_INDEPENDENT_COVERAGE_AND_ACCOUNTING'
    assert real['status'] == 'PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY'
    assert band['witness_sha256'] == audit['witness_sha256'] == real['witness_sha256'] == SHA
    assert band['objective_exact'] == direction
    assert band['first_even_spin'] == real['first_even_spin'] == 42
    assert band['last_even_spin'] == real['last_even_spin'] == 448
    assert band['energy_domain'] == ['4', 'infinity'] and real['mode'] == 'finite'
    assert audit['real_cells_checked'] == len(real['cells']) == len(band['cells']) == 4096
    assert audit['spins_per_real_cell'] == 204 and audit['paired_finite_uv_cells'] > 0
    # The independently rebuilt coverage ledger must bind this precise cost packet.
    audit_manifest = {str((ROOT/Path(p)).resolve()): h for p, h in audit['sha256_manifest'].items()}
    assert audit_manifest[str(band_path.resolve())] == sha(band_path)
    assert audit_manifest[str((BASE/'real/finite_split_replay.json').resolve())] == sha(BASE/'real/finite_split_replay.json')
    total_finite = F(0)
    cursor = F(0)
    region_costs = {'compact': F(0), 'remainder': F(0)}
    region_norms = {'compact': F(0), 'remainder': F(0)}
    for c, r in zip(band['cells'], real['cells']):
        lo, hi = F(c['lo']), F(c['hi'])
        assert cursor == lo == F(r['lo']) < hi == F(r['hi']) <= 1
        assert c['region'] == r['region'] and c['spins_checked'] == 204
        amount = max(F(r['first_upper_exact']), F(r['replay_upper_exact']))
        assert amount == F(r['retained_upper_exact']) == F(c['A_upper_exact']) >= 0
        kappa, weight = F(c['kappa_exact']), F(c['reciprocal_weight_upper_exact'])
        assert kappa > 0 and weight >= 1/(2*kappa)
        cost = amount*weight
        assert cost == F(c['cost_upper_exact'])
        total_finite += cost
        region_costs[c['region']] += cost
        region_norms[c['region']] += amount
        cursor = hi
    assert cursor == 1
    assert all(region_norms[k] == F(real['A_upper_exact'][k]) for k in region_norms)
    assert all(region_costs[k] == F(band['cost_by_real_region_exact'][k]) for k in region_costs)
    assert total_finite == F(band['cost_upper_exact']) == F(audit['cost_upper_exact'])
    removed = parts.pop('finite_42_448_compact') + parts.pop('finite_42_448_UV_and_sliver')
    assert 0 <= total_finite < removed
    parts['finite_42_448_banded_all_energy'] = total_finite
    total = sum(parts.values(), F(0))
    assert total == F(old['support_upper_exact']) - removed + total_finite

    dx, dy = map(F, direction)
    norm2 = dx*dx + dy*dy
    scale = 2**256
    normlo = F(math.isqrt((norm2.numerator*scale*scale)//norm2.denominator), scale)
    assert normlo > 0 and normlo**2 <= norm2 < (normlo + F(1, scale))**2
    assert str(norm2) == old['direction_norm_squared_exact'] and str(normlo) == old['direction_norm_lower_exact']
    outward = ceiling_decimal(total, 9)
    normalized = ceiling_decimal(F(outward)/normlo, 9)
    for p in [Path(__file__), ROOT/'research/six_targets_20261007/assemble_upper95_direct_v2.py',
              ROOT/'research/boundary_length_20261007/assemble_complete_v3.py',
              ROOT/'research/six_targets_20261007/REVIEW_UPPER95_DIRECT_BANDED_ASSEMBLY.md',ROOT/'research/six_targets_20261007/REVIEW_UPPER95_PARALLEL_ASSEMBLY.md']:
        inputs.bind(p)
    result = dict(old)
    result.update(
        inequality=f'({direction[0]}) g0+({direction[1]}) g2 <= {outward}',
        support_upper_exact=str(total), support_upper_outward_9dp=outward,
        support_normalized_upper_outward_9dp=normalized, support_normalized_float=float(normalized),
        components_exact={k: str(v) for k, v in parts.items()},
        components_float={k: float(v) for k, v in parts.items()},
        tail_upper_exact=str(total-parts['head_0_40']), tail_upper_float=float(total-parts['head_0_40']),
        finite_tail_method='Same-witness local coercivity minima across every overlapping closed energy cell and all204 finite spins; exact independently rebuilt cell ledger.',
        strengthening=dict(global_certificate=str(args.global_certificate.resolve().relative_to(ROOT)),
                           global_certificate_sha256=sha(args.global_certificate),
                           replaced_finite_components_exact=str(removed),
                           replacement_finite_component_exact=str(total_finite),
                           improvement_exact=str(removed-total_finite), improvement_float=float(removed-total_finite)),
        sha256_manifest=inputs.manifest)
    args.out.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: result[k] for k in ['status', 'support_normalized_float', 'components_float', 'tail_upper_float', 'strengthening']}, indent=2))


if __name__ == '__main__':
    main()
