"""Cache numerical V/p guards for all even spins 10..40 on a dense energy grid.

The output is a numerical guard bank, not a continuum sign certificate.
Normalized total V/p = a*v0 + b + einsum('i,ijk->jk', lam, rowVnorm).
The row order is exactly that of the frozen lower_enriched_cache.npz.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
from pathlib import Path
import argparse, hashlib, json, sys, time
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'research'))
import numpy as np
from angle_audit import normalized_inventory


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--order', type=int, default=96)
    ap.add_argument('--output', type=Path,
                    default=ROOT/'output/full_direction_20261007/left/dense_v_guards_l10_40_mu6p005_20p005_k96.npz')
    args = ap.parse_args()
    started = time.perf_counter()
    base = ROOT/'output/lower_g0_upper_g2_20261006/lower_enriched_cache.npz'
    with np.load(base) as old:
        names = np.array(old['names'], copy=True)
    assert len(names) == 220
    spins = np.arange(10, 41, 2, dtype=int)
    # Rational-grid specification, each node rounded once to binary64.
    mu = np.array([(1201+10*k)/200 for k in range(281)], dtype=float)
    tic = time.perf_counter()
    u0, v0, p, rowU, rowVnorm = normalized_inventory(
        names, mu, spins, kernel_order=args.order, angular_order=args.order)
    build_seconds = time.perf_counter()-tic
    assert rowVnorm.shape == (220, 16, 281)
    for value in (u0, v0, p, rowU, rowVnorm):
        assert np.isfinite(value).all()
    assert (p > 0).all()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, names=names, spins=spins, mu=mu,
                        u0=u0, v0=v0, p=p, rowVnorm=rowVnorm,
                        kernel_order=np.array(args.order), angular_order=np.array(args.order))
    sources = ['research/angle_audit.py', 'research/roy_kernels.py',
               'research/transfer_kernels.py', 'research/crossing_kernels.py',
               'research/forward_transfer_derivatives.py',
               'research/existing_solver/massive_forward_kernels.py']
    record = dict(status='FINITE_NUMERICAL_GUARD_BANK', rows=len(names),
        guard_points=int(len(spins)*len(mu)), spins=spins.tolist(),
        energy_grid='mu_k=(1201+10*k)/200, k=0,...,280; binary64 evaluations',
        mu_min=float(mu[0]), mu_max=float(mu[-1]),
        rowVnorm_shape=list(rowVnorm.shape), kernel_order=args.order,
        angular_order=args.order, formula="V/p=a*v0+b+einsum('i,ijk->jk',lambda,rowVnorm)",
        u0_normalization='Unnormalized real g0 kernel; not U/p.',
        v0_normalization='Imaginary g0 kernel divided by positive g2 kernel p.',
        flattening='C order: spin major, energy minor. rowVnorm.reshape(220,-1).',
        build_seconds=build_seconds, total_seconds=time.perf_counter()-started,
        cache_sha256=sha(args.output), inventory_cache_sha256=sha(base),
        producer_sha256=sha(__file__),
        source_sha256={path:sha(ROOT/path) for path in sources},
        claims='Numerical pointwise guards only. No optimization, integration, interval enclosure or sign certificate.')
    args.output.with_suffix('.json').write_text(json.dumps(record,indent=2), encoding='utf8')
    print(json.dumps(record,indent=2), flush=True)


if __name__ == '__main__':
    main()
