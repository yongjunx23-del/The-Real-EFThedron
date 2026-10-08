# Scoped tail certification, 7 October 2026

## Scope and current status

This directory contains additional **components**, not a complete certificate
for either competitive 220-row witness. Frozen earlier files are unchanged.

For each of the following exact binary64 witnesses, the new calculation proves
`V_l(mu)/p_l(mu) < 0` for every even `l >= 10000` and

```
mu >= 4 + l(l+1)/256,
```

including the endpoint at infinity by continuous compactification.

| Witness | SHA256 | Objective | Accepted cells |
|---|---|---|---:|
| `output/full_direction_20261007/right/g2_guarded192_refined.npz` | `825dcff5caca223ae757bc027b2e1438d80e9878f4c2801db80e2a5790452ca7` | `g2` | 1342 |
| `output/full_direction_20261007/left/angle180_final.npz` | `e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b` | `-g0` | 1472 |

Both adaptive calculations were performed afresh at Arb280 and Arb360. The
partitions, exact rational outward bounds, source data, and dependency hashes
are checked by `verify_scoped_packet.py`. The resulting ledger is
`output/all_direction_certification_20261007/tail/scoped_uv_packet_verified.json`.
The largest directed cell upper bounds are approximately `-1.16867e-6` and
`-1.81830e-5`, respectively. These are partition bounds, not the actual maxima.

The following are **not covered** by this packet:

* spins below10000;
* the complementary region `tau=l(l+1)/(mu-4)>256`;
* any bound on the integrated real-kernel norm;
* a complete dual support integral or coupling inequality.

The source assumption, dispersion conventions and disk relaxation are exactly
those of the existing massive4D construction. In particular,
`c=4/3`, `Q=mu-c`, `Z=Q/(mu-4)`, `p=P_l(Z)/Q^3` and
`dnu=(2l+1)dmu/(2pi sqrt(1-4/mu))`.

## Finite-spin proof, rather than a joint-limit sample

`prepare_components.py` rebuilds exact rational joint-source moments,
16×16 Bernstein-subdivided combined-source absolute norms, and the crossing
polynomial certificate. Source weights are the combined signed polynomial
densities, not a triangle inequality on individual basis coefficients.

The combined norms are approximately

| Witness | `integral abs(W K1_infinity)` | `integral abs(W K2_infinity)` |
|---|---:|---:|
| g2 | 309.13985148 | 71.43483956 |
| left | 1654.88433604 | 426.56608838 |

All proof calculations use the exact rational fields in the saved JSON.
`finite_uv_component.py` applies the reviewed finite-energy enclosure in
`research/lower_g0_upper_g2_20261006/tail_finite_uv_region.py`, with
`T=l(l+1)`, `eta=1/T`, `w=1/(mu-4)=tau eta`. It retains the exact finite-energy
forward and derivative jets, the first64 positive Legendre terms, the exact
crossing polynomial, the source Mehler error and the rational kernel errors.
The exact polynomial crossing certificate permits omission of terms `k>=64`
because every such term is nonpositive. The interval `eta in[0,1/(10000*10001)]`
covers all the required integer spins; it is not a sampled-spin interpolation.

An important adaptation is explicit addition of the objective coefficient
`b` to `V/p`. The historical finite-UV helper was written for `(-1,0)` and did
not itself add `b`. Here the repaired g2 objective contributes **exactly +1**.
The two prepared objective arrays are checked in the final packet verifier.

The failed L1000 diagnostic envelopes are saved as `*_screen_*.json`; they
are only pointwise screens and have no certification status. Near tau2.25
the finite interval upper bounds are positive, while the L10000 envelopes
pass. The source norms are already combined-density Bernstein bounds. The
main residual difficulty is cancellation in large derivative/harmonic
coefficients under a finite eta interval, together with shallow joint-limit
sign margins. A lower cutoff therefore needs a materially sharper finite
parameter enclosure, not just more angular samples.

## Reproduction

From the workspace root, using the configured Python environment with
python-flint under `tmp/certification_deps`:

```powershell
python research/all_direction_certification_20261007/tail/prepare_components.py output/full_direction_20261007/right/g2_guarded192_refined.npz --label g2
python research/all_direction_certification_20261007/tail/prepare_components.py output/full_direction_20261007/left/angle180_final.npz --label left
python research/all_direction_certification_20261007/tail/finite_uv_component.py --label g2 --lmin 10000 --precision 280
python research/all_direction_certification_20261007/tail/finite_uv_component.py --label g2 --lmin 10000 --precision 360
python research/all_direction_certification_20261007/tail/finite_uv_component.py --label left --lmin 10000 --precision 280
python research/all_direction_certification_20261007/tail/finite_uv_component.py --label left --lmin 10000 --precision 360
python research/all_direction_certification_20261007/tail/verify_scoped_packet.py
```

Preparation takes about3.4seconds per witness and each complete adaptive UV
run takes about8–10seconds on this machine. Rebuilding the pickle files can
change their byte hashes even if the mathematical content is unchanged, so
rebuild the downstream records and packet together when preparing afresh.

## Additional old150 components for slopes4 and6

For the already certified old150 witness with SHA
`bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c`, take
objective `g0+b*g2` and change only its X4 multiplier by `(18/5)b`.
The source polynomials and U kernel remain unchanged. Two new analytic
Legendre estimates improve the compact and UV envelopes.

Write `x=Z-1=2c/(mu-4)` and let the normalized positive Taylor coefficients
of `P_l(1+x)` define a distribution on `k`. Set
`m=E[k]=x P_l'(Z)/P_l(Z)`. Exactly,

```
chi = X4/p
    = E[k(k-1)]/(2c^2) - 3 E[k]/(2c Q)
   >= [m^2 - m(1+4/Q)]/(32/9).
```

This is Jensen's inequality `E[k^2]>=E[k]^2`, not a numerical fit. The
quadratic is increasing in `m` once `m>(1+4/Q)/2`; code checks that condition
before using it. Otherwise the known nonnegative lower bound is retained.

For a compact-energy cell and all `l>=L`, the moment representation
`P_l(Z)=pi^-1 integral_0^pi (Z+sqrt(Z^2-1)cos(theta))^l dtheta` implies
log-convexity in `l`. The standard derivative recurrence therefore gives

```
m_l/l = [Z-P_(l-1)(Z)/P_l(Z)]/(Z+1) >= m_L/L.
```

Also `dm_L/dlog(x)=Var(k)>=0`, so the endpoint with largest mu gives a valid
lower bound throughout the cell. The resulting quadratic lower bound in
`l` supplies both the value and the spin derivative used to close the
unbounded-spin domain in `compact_jensen.py`.

At fixed `tau=l(l+1)/(mu-4)`, the positive coefficients are

```
(c*tau)^k/(k!)^2 * product_(j<k) [1-j(j+1)/(l(l+1))].
```

Increasing spin increases the ratio of coefficients monotonically with k,
and adds only higher-index terms. Hence the resulting distributions are
ordered by monotone likelihood ratio and `m_l(tau)>=m_L(tau)`. Combined
with monotonicity in tau, this supplies the UV lower bound for chi, including
the entire unbounded region `tau>=256`. `positive_uv_jensen.py` retains the
negative X4 contribution there instead of discarding it.

The support estimate is the already reviewed disk inequality

```
h(U,V) <= U^2/(2 eta p) + 2p(V/p+eta)_+,
```

and the old fully certified `A40` real norm plus the all-spin positive-band
mass estimate. The new source and real norms are not inferred from the220
witness; this subsection uses only the old150 witness and its hash-bound
component packet.

Fresh Arb280/360 runs yield the following **scoped** results:

| b | X4 increment | Compact sign domain | Compact V/p upper | UV support, all even l>=42 and mu>=49 |
|---:|---:|---|---:|---:|
| 4 | 72/5 | every even l>=46, 4<mu<=49 | -31.70326549 | <0.011315489 |
| 6 | 108/5 | every even l>=108, 4<mu<=49 | -63.29990063 | <0.061054370 |

The exact records are in
`output/all_direction_certification_20261007/tail/old150_extension/partial_components_{280,360}.json`.
The compact bridges `l=42,44` for b4 and `l=42,44,...,106` for b6 are
**not covered**. Neither result is a complete directional bound; the finite
head components and full assembly must also be supplied. There is no
inference that a positive loose envelope represents an actual violation.

```powershell
python research/all_direction_certification_20261007/tail/old150_partial_components.py --precision 280
python research/all_direction_certification_20261007/tail/old150_partial_components.py --precision 360
```

Each run takes about20seconds. These new Jensen extensions have not received
an independent implementation or external mathematical audit. The code
uses exact rational endpoints and directed Arb arithmetic. The higher-
precision replays are checks of the same implementation, not formal proof
verification.
