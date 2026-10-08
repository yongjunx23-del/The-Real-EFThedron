# Twelve-null full-unitarity Gram head certificate

The complete finite head, every even spin from zero through 18 and every energy above threshold, has the exact rational support upper bound

**11.693432984.**

This is a component of the full coupling bound. The independently certified tail starting at spin 20 must be added; it is not set to zero.

## Exact witness

The objective kernels are `(U,V)=(0,p)+sum lambda_a*(U_a,V_a)`. Thus their integral is `g2`. The exact multipliers are:

| Null equation | Multiplier in the objective |
|---|---:|
| `Dt1G_regular` | `-4591719/20000000` |
| `H2` | `-151061389/100000000` |
| `H3` | `-60956417/50000000` |
| `H4` | `7391707/100000000` |
| `H5` | `32663847/100000000` |
| `H6` | `1740043/10000000` |
| `H7` | `15796213/100000000` |
| `H8` | `-9007949/100000000` |
| `Dt1H_j2` | `13791899/100000000` |
| `Dt1H_j3` | `6721549/50000000` |
| `Dt1H_j4` | `-12089969/100000000` |
| `X4` | `-82836871/20000000` |

The signs are the negatives of the `raw_null_weights` in the high-spin sign certificate. The witness file records the hash of that source. Rational rounding is part of the definition of this witness; it is not an unaccounted numerical perturbation.

## Finite-head results

| Spin | Outward support upper bound | Rational velocity cells | Polynomial degree |
|---:|---:|---:|---:|
| 0 | 5.420775896 | 835 | 44 |
| 2 | 2.406588483 | 1314 | 48 |
| 4 | 3.400237797 | 550 | 52 |
| 6 | 0.380467973 | 602 | 56 |
| 8 | 0.060629224 | 466 | 60 |
| 10 | 0.014920805 | 374 | 64 |
| 12 | 0.005464992 | 273 | 68 |
| 14 | 0.002470140 | 204 | 72 |
| 16 | 0.001234513 | 148 | 76 |
| 18 | 0.000643165 | 125 | 80 |

There are **4891** cells and **273411** positive-semidefinite Bernstein coefficient matrices. Their exact integral sum is

`2658244733366035013/227328000000000000`.

The main certificate hash is `42a2ece2500ced7afd629d65f3a62c1521f3d80201c828b333f1632cd464f330`.

## What the matrix condition is

For each energy and spin, the unitarity matrix is `W=[[b,r],[r,2-b]]>=0`. We construct the entire dual matrix `C=[[h-2V,-U],[-U,h]]>=0`. Pairing them gives `h-U*r-V*b=tr(C*W)/2>=0`.

The full `2-b` entry is retained. On each velocity cell, `h=(1-beta^2)^2*H`, where `H` is a recorded nonnegative rational. Clearing positive denominators and applying a positive congruence produces a rational polynomial matrix. Every Bernstein coefficient matrix is checked positive semidefinite exactly, proving positivity over the continuum of energies.

The reference integral of this piecewise function is exact: each cell contributes `4*(2*ell+1)*H*(beta_hi-beta_lo)/pi`. The exact inequality `333/106<pi` gives the final rational outward rounding. No numerical quadrature is used to accept the head bound.

These polynomial matrix certificates have explicit interval-SOS Gram representations. Their largest pair has dimensions **82 by 82 and 80 by 80**. The full dense rational matrices for a representative spin-18 cell are in `certificate_m011_sos_gram_82.json`; their exact PSD decomposition and coefficient identity are included. The compact main certificate reconstructs the corresponding large matrices on every cell.

## Validation

- `certificate_replay_head.py` rebuilds every matrix by a separate homogeneous-Horner/projective Bernstein transformation, checks every diagonal and determinant exactly, checks gapless partitions, and reconstructs the exact support sum. Its completed result is `certificate_m011_head_replay.json`.
- `certificate_check_kernels.py` checks the rational kernel construction against the existing independently differentiated transfer-kernel implementation at 70 spin/energy points. The maximum relative-scaled discrepancy was approximately `3.83e-13`. This is an implementation cross-check, not the positivity proof.
- `certificate_export_sos_gram.py` constructs the actual large Gram matrices, proves their PSD property by exact decomposition, and matches every polynomial coefficient exactly.

## Approximations and scope

The multiplier search was numerical, but the displayed rational multipliers define the final witness exactly. Search convergence or optimality is not required for validity. Reference polynomials were selected numerically and then accepted only by exact rational PSD checks. The Bernstein/SOS ansatz is conservative; its approximation error appears as a larger objective, never as an omitted error term.

The twelve null equations are a finite subset of the manuscript's valid crossing/harmonic/transfer equations. Choosing a subset weakens the bound and does not impose an amplitude ansatz. The same pole-free amplitude class, unitarity normalization, dispersion assumptions and locally uniform convergence required for first transfer derivatives remain necessary.

The original finite-spin search is not itself a physical truncation: the final result must combine this head with the separately proved all-spin tail. The head certificate alone does not establish a bound on the complete `g2`.

This construction is a family of interval matrix-SOS Gram certificates. It should not be described as a single global monomial Hankel moment matrix or as new physical information beyond the retained valid equations.

## Files

- `certificate_m011_witness.json`: exact objective multipliers and source hash.
- `certificate_m011_head.json`: full rational polynomial/cell certificate.
- `certificate_m011_head_replay.json`: separate exact arithmetic replay.
- `certificate_m011_head_kernel_check.json`: kernel normalization cross-check.
- `certificate_m011_sos_gram_82.json`: actual 82 by 82 and 80 by 80 rational Gram matrices.
- `certificate_general_head.py`: producer for the general rational first-derivative/harmonic family.
- `certificate_large_sos_gram_proof.md`: explicit conversion to large interval-SOS Gram matrices.
