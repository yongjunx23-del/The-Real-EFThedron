# Review of the complete left-bound assembly

Reviewed by the tail agent on 2026-10-07, separately from the root implementation.
This is an internal mathematical and source-code review, not external peer review
or formal verification. The complete bound still requires the component records
to finish and their exact hashes to be checked by the assembler.

Files reviewed:

- `research/all_direction_certification_20261007/assemble_competitive_left.py`
- `research/all_direction_certification_20261007/real220_banded_coercivity.py`
- `research/all_direction_certification_20261007/real220_local_coercivity.py`

No mathematical sign or normalization defect was found in the inspected code.
The following are the substantive checks.

1. The fixed witness has objective `(-1,0)` and SHA-256
   `e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
   Therefore the sum of support upper bounds is an upper bound on `-g0`, and
   the final lower bound is correctly written as `g0 >= -C`.
2. The exact change of variables is
   `mu = 4/(1-beta^2)` and
   `tau = ell(ell+1)(1-beta^2)/(4 beta^2)`.
   The mapper reverses beta endpoints when taking the tau interval. The energy
   boundary `mu=49` is tested using `beta^2=45/49`, with rational comparisons.
   No rounded square-root endpoint is used to select a theorem.
3. Every real-norm rectangle takes the minimum coercivity among all intersecting
   closed spectral intervals. Including both intervals when an endpoint is
   shared is conservative. For a rectangle straddling `mu=49`, both compact and
   UV margins are included. For tau above256 the large-tau theorem is used;
   the finite-spin proof stops at `min(256,ell(ell+1)/45)`. The seam at ell108
   is correct because `106*107/45 < 256`.
4. If `V <= -kappa*p`, the disk support obeys
   `h(U,V) <= U^2/(2*kappa*p)`. Multiplying each real-norm cell by
   `1/(2*kappa)` is therefore valid. Taking the minimum over all spins before
   multiplying the aggregated mass is conservative. The separate per-spin
   weighted implementation uses the same inequality before summation.
5. The finite block42..448 and infinite block450..infinity are disjoint and
   exhaust all even spins above40. Together with the all-energy head0..40,
   the assembly covers every even spin. The infinite block uses compact
   coercivity on the compact real-norm cells and the weaker UV coercivity on
   the remaining cells, including the deliberately retained small compact
   sliver. This is valid because the asserted compact margin is at least the
   UV margin.
6. Producer and replay upper bounds are combined using the exact rational
   maximum. The banded producer rounds the positive reciprocal `1/(2*kappa)`
   upward to40 significant binary bits: if `e=floor(log2(x))`, its spacing is
   `2^(e-39)` and the stored value is `ceil(x/spacing)*spacing`. Both the ceiling
   and its inequality check use exact fractions. Its relative excess is less
   than `2^-39`. The assembly independently checks the stored reciprocal is
   no smaller than `1/(2*kappa)` and its product with the mass is exact.
   Sums and final decimal rounding use exact fractions. The
   outward six-decimal result is checked to be no smaller than the exact sum.
7. The weighted strip bound's selection between two individually valid upper
   bounds using float priority cannot make the result invalid: either selected
   bound is valid independently. Acceptance tests and coercivity maxima do not
   use float ordering.

The final inspected assembly explicitly asserts the head record spin list is
exactly `0,2,...,40`, and binds the additional in-memory rational source rebuild
at
`output/all_direction_certification_20261007/tail/left/refined/finite_source_exact_rebuild.json`:
all972 exact moments and9 source-remainder coefficient bounds reproduce the
frozen cache exactly. No source or compact cache bytes were modified.

For the later family with exact objective `(-1,-b)`, `b>=0`, the null witness
is unchanged and `V` becomes `V-b*p`. Thus every local coercivity increases
by exactly `b`, while `U` and the real-norm records remain unchanged. The head
must still be integrated or bounded for the new objective. Exact positive
combinations of witnesses are also admissible; rounding their coefficients
would require a new error analysis.
