# Internal review of the new boundary-length campaign envelopes

Reviewed on 2026-10-07. This is an internal mathematical and implementation
review, not an independent external review or formal verification. A review
of a formula does not substitute for the saved interval cover of its domain.

## Angle 255: retaining every negative shifted crossing power

Reviewed `tail/certify_angle255_large_tau.py`, including its 320- and 400-bit
records. The witness SHA256 is
`5b7febef7cb467ff224cd7162a82d76a3534fd3439102f459733e513b4f31615`.

Write the crossing polynomial as `F(k,z)` and shift `k=y+10`. On the entire
exact `z` interval `[0,3/8]`, each coefficient of `y^j` is bounded above by
the maximum of its exact rational Bernstein coefficients, `c_j<=0`. Set
`a_j=-c_j`. For integer `k>=10`, this gives

    F(k,z) <= -sum_{j>=1} a_j (k-10)^j.

The constant term is nonpositive and may be discarded in an upper bound.
For `k=0,...,9`, the separately computed exact Bernstein cap `M>=0` gives
`F(k,z)<=M`. Combining both domains therefore gives

    F(k,z) <= M - sum_{j>=1} a_j (k-10)_+^j.

The normalized positive Legendre coefficients define a probability measure
on the integer index `k`. For each integer `j>=1`, `(k-10)_+^j` is convex,
including `j=1`. Jensen's inequality can be applied to each term separately:

    E F <= M - sum_{j>=1} a_j (E k-10)_+^j.

There is no assumption that the crossing polynomial itself is convex.
There is also no Jensen step for `j=0`.

The reviewed positive-coefficient stochastic ordering gives
`E k >= m_108(256)>10` for every physical spin `ell>=108` and joint coordinate
`tau>=256`. Every retained negative-power contribution therefore becomes at
least as negative when the mean increases. The previous positive monomial
degree-20 bounds on forward and source terms still apply; the forward
numerator has degree at most three, so all powers `tau^(j-20)` decrease on
this domain. The exact objective's second component is added to the forward
constant, and its first component is already included in the source-prepared
forward coefficient. Neither component is replaced by a nominal angle.

The saved interval upper bound is below -816, comfortably below the claimed
`-1`. Both precisions agree on this conclusion. No mathematical defect was
identified in this strengthening of the old fifth-power-only estimate.
It covers the unbounded-tau UV component only, not the full functional cost.

## Angle 85: compact finite-spin evaluator

Reviewed `compact/compact_finite220_combined.py` against the previously
reviewed left-witness evaluator. The witness SHA256 is
`58fafda2eafa3a741610175f416a009c84f9a18f309d36e74d0e4007f65fd6d8`.

The source polynomial division, source geometric remainders, removable pole
treatment, and derivative terms are unchanged. The replacement source tables
have an independent exact-beta audit: all 333 beta identities and all 56,934
overlapping source moment entries agree exactly with the independently
rebuilt finite-head tables. These checks are saved separately and must be
bound into the final certificate manifest.

The exact objective `b` enters the combined crossing numerator as
`b*Rstar*Q^6`; division by `Q^6*Rstar` contributes precisely `b` to `V/p`.
The objective `a` enters `ForwardCombined.alpha` through the unchanged exact
source/subtraction identities. The higher crossing rows remain present.

Moving `ForwardCombined.core(beta)` outside the spin loop is legitimate:
inspection of `Forward220.core` shows it uses only beta, alpha, and the
forward polynomial/rational coefficients, not `ell` or its Legendre
coefficients. The derivative and Legendre terms remain inside the spin loop.
It evaluates the identical interval expression fewer times and does not
alter the enclosure. No formula defect was identified in this adaptation.

The analytic compact driver uses the reviewed sharp envelope with the exact
new objective and witness. Its original saved coverage is every even spin
at least 42 on `4<mu<=19/4`, and every even spin at least 112 on `4<mu<=49`.
The finite bridge and any later piecewise improvements require their own
complete coverage and higher-precision replay records. This review does not
assert that those numerical jobs are finished.
