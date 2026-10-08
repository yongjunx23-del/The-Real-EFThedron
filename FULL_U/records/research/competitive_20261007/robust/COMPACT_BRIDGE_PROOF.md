# Compact finite-spin coercivity for the robust upper-g2 witness

This component concerns only the frozen witness
`output/competitive_20261007/robust/g2_robust_stage1.npz`, SHA256
`933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92`,
with objective exactly `(0,1)`. Multipliers are interpreted as their exact
binary64 rational values. Transfer endpoints and named parameters have the
same exact rational conventions as the head and tail proofs.

The completed packet proves `V_l(mu)/p_l(mu) < -1` for every even spin at
least 42 and every `4 < mu <= 49`. It is a coercivity component, not a whole
coupling bound. Integrating the real-part norm and covering higher energies
are separate obligations.

## Finite-spin algebra

Put `delta=mu-4`, `c=4/3`, and `Q=mu-c`. Define

`R_l(delta,t)=sum_j binom(l,j) binom(l+j,j) t^j delta^(l-j)`

and `Rstar=R_l(delta,c)=delta^l P_l(1+2c/delta)`. This is positive on the
physical interval. All direct finite-spin cells start at `mu=5`, so no
threshold endpoint division occurs in this evaluator.

The Roy source has angular degree at most six. Its diagonal PV residue is
therefore zero for every spin used here. Polynomial division gives its
remaining PV term exactly: the coefficient of `delta^p` is

`sum_(j=0)^(l-1-p) C_lj M_s[j,l-j-1-p]`.

For every source the divided difference

`(Rstar-R_l(delta,t))/(delta+c+t)`

is a polynomial because even-spin reflection makes the numerator vanish at
`t=-delta-c`. Long division is performed algebraically; integration uses
the exact polynomial source moments `T_j`. No sampled removal of either
pole is used.

The two nonsingular source denominators are expanded geometrically through
order 255, using moment tables through order 256. These tables enclose
degree-exact polynomial integrals with certified Gauss nodes and weights.
Their construction has no quadrature truncation error. The geometric
residual remains explicit. If `r_u`, `r_t` are the source half-widths and
`h`, `t_m` its expansion centers, its absolute normalized residual is at
most

`Q^3 ||W||_1 [rho_l q_u^256/((mu+h)(1-q_u)) + q_t^256/((mu+t_m)(1-q_t))]`,

where `q_u=r_u/(mu+h)` and `q_t=r_t/(mu+t_m)`. Every denominator is positive
on this domain. The physical-angle factor is `rho_l=1`; for the small
nonphysical Roy subregion the evaluator instead bounds the scaled
Legendre polynomial at its extreme angular endpoint and divides by a
strictly positive interval lower bound for `Rstar`. A cell is rejected if
that lower bound cannot be established. The total certified strip norm is
used for each strip separately, which is conservative.

Roy polynomial coefficients above degree `l+2N-2` vanish analytically by
the source orthogonality identities. The corresponding computed Arb
coefficients must contain zero before those known-zero coefficients are
removed. Tapered strip PV terms retain their exact Legendre-Q expressions,
including the continuous endpoint formulas in terms of `x log|x|`.

The source evaluator inherits these formulas and error terms from the
reviewed `certified_finite_v.Source.eval`; the new constructor supports the
larger spin range and new witness-specific moment tables.

## Forward, derivative, and crossing rows

The unchanged reviewed `ForwardCombined` evaluator supplies all forward
harmonic, rational, and transfer-derivative rows. No row is discarded.
Derivative expressions use their factored threshold form, although this
finite bridge stays strictly above threshold.

The exact crossing numerator is combined before interval evaluation. Let
`Rstar_c` and `Rstar_cc` denote the derivatives with respect to `c` and
`cross=(Q Rstar_cc-3 Rstar_c)/2`. The X4 term is
`lambda_X4 cross/(Q Rstar)`. Its contribution to the common denominator
`Q^6 Rstar` is consequently `lambda_X4 cross Q^5`.

The objective contributes `Rstar Q^6`. Each higher crossing row of degree
`D`, with derivative index `j` and coefficient `beta_j`, contributes

`lambda beta_j [sum_(k>=j) C_lk binom(k,j)c^(k-j)delta^(l-k)] Q^(6-(D-2-j))`.

Every exponent is a nonnegative integer for this bank. Combining these
polynomials exactly substantially reduces interval dependency; it does not
alter the functional.

An early, unpublished diagnostic version omitted a factor `k` in
`Rstar_cc`. That version produced a large point mismatch and was rejected
before any certificate was generated. The chosen source uses the complete
factor `k(k-1)`. The corrected combined model agrees with independent
floating kernels at spins 42 and 108 across `mu=4.5..49` to better than
`8.5e-13` relative in the saved point check. This comparison diagnoses
implementation mistakes; interval enclosures, not that agreement, certify
the spectral continuum.

## Exact interval coverage and fresh replay

The independently computed energy pieces are six half-unit intervals
covering `[5,8]`, followed by `[8,12]`, `[12,24]`, `[24,36]`, `[36,48]`,
and `[48,49]`. A separate unsplit `[5,8]` run is retained as an unused
diagnostic. For every even spin from 42
through 108, bisection saves exact rational interval endpoints. Acceptance
requires both a finite Arb value and an exact directed upper bound below
`-1`. Float values affect neither acceptance nor the final coverage test.

The producer works at 512 bits. A fresh 640-bit replay evaluates every
original spin/interval pair. Identical energy intervals are grouped only
to share calculations across spins; none is omitted or widened. The exact
partition of each spin must remain unchanged. The union checker enforces
ordered positive-width cells, identical neighboring boundaries, endpoints
5 and 49, the full even-spin list, witness identity, and input hashes.

Both runs use the same mathematical implementation. This is a fresh
higher-precision replay, not an independent implementation or a formal
proof-assistant check. Formula review, separate point comparisons, exact
moment audits, and the inherited component reviews complement it.

## Completing all compact spins

Two separate analytic envelopes close the other regions:

* `low_allspin42_400.json` and `low_allspin42_512.json` cover every even
  spin at least 42 on `4<mu<=5`. The sharp compact envelope uses
  `Qmax=11/3`, the positive Legendre mean at `delta=1`, and its negative
  spin derivative. Reusing the larger conservative source-resolvent
  constants from the original envelope only weakens this bound.
* `compact_highspin_sharp_L110_320.json` and its 400-bit replay cover every
  even spin at least 110 on `4<mu<=49`. The sharp complex-transfer bound
  and fifth-power crossing envelope are reviewed in the associated
  all-spin proof.

The three regions have overlapping closed seams at `mu=5` and no missing
even spin between 108 and 110. The assembly hashes all component packets
and their input manifests and records the exact coercivity constant one.
The endpoint `mu=4` is not asserted as an ordinary finite value; the
analytic low-energy theorem covers every point strictly above it.

## Reproduction

From the workspace root, with the project's Python/FLINT environment:

```text
python research/competitive_20261007/robust/run_compact_slices.py
python research/competitive_20261007/robust/replay_compact_slices.py
python research/competitive_20261007/robust/run_compact_low_slices.py
python research/competitive_20261007/robust/replay_compact_parallel.py --input output/competitive_20261007/compact/slice_5_5.5_512.json --out output/competitive_20261007/compact/slice_5_5.5_640_parallel.json --jobs 6
python research/competitive_20261007/robust/merge_compact_slices.py
python research/competitive_20261007/robust/assemble_compact_coercivity.py
```

The first four commands can overlap: replay waits until each producer
packet exists, and the low-energy half-unit jobs replay immediately. The
selected replay of `[5,5.5]` uses six disjoint groups of original intervals
in separate processes; exact partition equality is checked after joining.
The serial replay of that piece is redundant and not selected. The final two commands must follow successful completion
of all producers and replays. Source moments and the paired all-spin
packets are prerequisites, with their exact hashes checked during
assembly. The producer may require adaptive subdivision even when a
pointwise plot appears comfortably negative.
