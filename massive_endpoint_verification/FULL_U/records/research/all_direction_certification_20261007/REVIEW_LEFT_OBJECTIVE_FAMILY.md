# Review of the exact negative-g2 objective family

Reviewed `head/replay_left_family.py` on 7 October 2026.
The fixed base witness has SHA256
`e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
The additional witnesses retain its exact multipliers and use objective
(-1,-b), for b=1/4,1/2,1,2,4,8. No rounded new coefficient vector is used.

Since g2=<p Im f> with p>0, the change is exactly U_b=U_0,
V_b=V_0-b*p. The disk support H(U,V)=sqrt(U^2+V^2)+V is nondecreasing
in V (including the continuous U=0 case). Thus every original cell, gap,
and complete tail upper bound remains a valid upper bound after the shift.
This justifies both fallback bounds and reuse of source/threshold gap caps.

The beta-compacted direct fields are U/e^2 and V/e^2, e=1-beta^2.
Let z=beta^2, dbar=4z, Qbar=4-(4/3)e and
R=sum_k binom(ell,k)binom(ell+k,k)dbar^(ell-k)[(4/3)e]^k.
Then p/e^2=e*R/(dbar^ell Qbar^3), exactly the new `pscaled` expression.
Therefore the direct and endpoint support expressions subtract b*pscaled,
and the beta Taylor jet subtracts the jet of that same rational function.
The objective shift has no source Taylor/geometric residual, so the old
source error allowance remains valid. The original convex-support
trapezoidal curvature estimate applies to the shifted fields with the
shifted second derivative included.

For the normalized negative-V branch, the shifted imaginary field is
(V/p)-b and the real field and positive normalization factor are unchanged.
The code retains the original branch applicability (including the special
ell0 beta<1/2 restriction). The cell support upper bound is the minimum of
several individually valid upper bounds. Choosing one by floating ranking
does not affect validity; the chosen bound is then taken upward. Reusing a
previous smaller-b support bound is valid by the same monotonicity.

The cache retains field balls, not a shifted objective's endpoint support.
The model's direction stays immutable. Interval/jet/normalized caches are
reset each cell; exact-point field caches remain reusable. There is no
attempt to treat a different witness hash as the source-table witness.
The family record binds the base witness, exact shifts, original source
partition, gaps, and implementation dependencies explicitly.

The complete spin partition is checked as an exact rational union of all
ordinary cells and gap intervals, covering beta0..1 without omissions.
This review identifies no formula or coverage defect in the new family
driver. It does not assert completed family bounds before every spin
record and a valid base-tail certificate have been assembled. The driver
is an interval computation with internal review, not a formal proof or
an independent external audit.
