# Review of exact objective shifts on the saved 268-degree witness

This internal code and mathematical review concerns
`research/weak_regions_20261007/head/shifted_head.py` and
`research/weak_regions_20261007/certify_tail_transfer_v2.py`.
It is a separate agent review, not an external review or a formal proof.

## Exact functional and transfer

The new NPZ keeps all 220 binary64 multipliers and the first objective
coordinate bit-for-bit identical to the saved 268-degree witness. The
second objective coordinate is a new binary64 number interpreted as an
exact dyadic rational. If its change is delta < 0, the exact fields obey

    U_new = U_base
    V_new = V_base + delta p
    p = P_l(1 + 8/[3(mu-4)]) / (mu-4/3)^3 > 0.

For even nonnegative spins and mu > 4, the argument of P_l is larger
than 1 and P_l is positive. The disk support H(U,V)=sqrt(U^2+V^2)+V
is nondecreasing in V. Thus every previously certified cell integral
and every previously certified singular cap upper remains valid.

No multiplication or rounding of the 220 coefficients is performed.
The exact objective is generally not a unit vector. The exact reported
inequality uses its unnormalized coordinates; Euclidean normalization
is only needed for plotting and comparison.

## Shifted finite-head bounds

The interval fields in `fields(a,b)` include the source remainder in the
V ball. The third return value is a list of logarithmic factors, not an
unaccounted source error. The ordinary cells avoid source endpoints.

`pscaled` evaluates p/e^2, with e=1-beta^2, matching the scaling of
the unchanged model's U and V fields. Adding delta times this quantity
is therefore the correct objective change in the ordinary fields.

Three new estimates are independently valid on each cell:

1. An interval Riemann bound of the shifted disk support.
2. Endpoint convex interpolation plus the unchanged Lipschitz/curvature
   remainder, with the shifted V jet and its second Taylor coefficient.
   The existing separate approximation error retains its factor 4h.
3. The existing normalized negative-V estimate, after replacing V/p by
   V/p+delta, only where the original method permits this branch.

The normalized branch below beta=1/2, or for positive spin between
beta=1/2 and beta=0.99, preserves the original applicability conditions.
`normalized_interval` selects its principal-value approximant from the
whole cell before evaluating range and midpoint jets. This prevents a
stale or mismatched midpoint source branch.

Each cell takes the smallest exact directed upper among these estimates
and the saved monotonicity fallback. Taking this minimum is valid because
every option bounds the same complete cell integral. Skipping evaluation
when a saved cell upper is below a specified cutoff only retains the
rigorous fallback; it never drops that cell's contribution.

The five singular gaps are reused by monotonicity, including the selected
constant threshold cap for spin 2. They remain included in the exact sum.
They are not newly reevaluated at each precision.

The assembly must check the complete base certificate and its dependency
hashes, exact cell/gap coverage, exact equality to the saved base partition,
the base source manifest and record hashes, exact shifted objective and
coefficient identities, exact per-cell sums, and both 256/384-bit results.
For each spin it retains the larger of the two complete integral upper
bounds. A maximum of per-cell values would be at least as large and is
unnecessary. All even spins 0 through 40 must be present.

## All remaining spins

The old base tail was reassembled with the unchanged version-3 verifier,
including reselection of every local coercivity margin. The transfer
reuses the exact bound A >= integral U^2/p on each finite-tail energy cell.
It replaces each positive margin kappa by kappa-delta and charges

    A / [2(kappa-delta)].

The reciprocal is rounded upward to a dyadic number before multiplication
and exact summation. The finite ledger covers every even spin 42 through
448 and the whole energy interval. The infinite-spin ledger covers every
even spin at least 450. Its remainder uses the minimum of the compact and
UV margins, because the real-norm cutoff leaves a small compact-energy
sliver in the remainder. These disjoint contributions are added once.

## Scope

This certifies new inequalities of the same underlying physical
relaxation when all component checks pass. It adds no physical constraint,
does not establish optimization convergence, and does not convert the
digitized published comparison curve into exact reference data.
