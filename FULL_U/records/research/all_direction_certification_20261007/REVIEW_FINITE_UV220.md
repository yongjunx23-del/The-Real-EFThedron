# Independent review of the finite-spin ultraviolet bridge

Parent-agent review, 2026-10-07. This is an internal mathematical and code
review, not an external audit or proof-assistant verification.

Witness: `g2_robust_stage1.npz`, SHA256
`933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92`.

Files reviewed:

- `tail/finite_source_corrections.py`
- `tail/finite_uv_bridge220.py`

No defect was found in the following construction.

The finite-energy source factor is a rational function R1=N(w)/D(w), with
the three exact denominator factors stated in the implementation. Expanding
the numerator verifies its coefficients through degree three. The recurrence
for the coefficients r_j is the exact polynomial-division recurrence.
For a truncation before order N, the residual numerator has coefficients
`r_N`, `-d2*r_(N-1)-d3*r_(N-2)`, and `-d3*r_(N-1)` at orders N through N+2.
Thus the implemented remainder is a polynomial identity; it does not assume
convergence of an infinite geometric series. The physical source ranges
give the positive denominator lower bound for `0<=w<=1/45`.

The exact moments include the folded factor two for Roy rows and the
taper factor for strip rows. Substituting the physical Roy source relation
gives the beta-integral factor `k! b!/(b+k+1)!` and energy denominator
`a+k+1`; for strips the energy denominator is `a+1`. These match the code.

For fixed spin, let `T=l(l+1)` and
`a_k=prod_(j<k)(1-j(j+1)/T)/(k!)^2`. Then
`P_l(1+2t/delta)=sum_k a_k*(t*tau)^k`, where `tau=T/delta`.
Combining this polynomial with `r_j*w^j`, where `w=tau/T`, gives the exact
source coefficient `(-1)^k*a_k*M[j,k]/T^j` multiplying `tau^(j+k)`.
The implementation therefore retains the finite-spin source dependence;
it does not replace the Legendre polynomial by a limiting Bessel function.

When l exceeds the stored moment degree 80, the omitted absolute Legendre
tail is bounded by

`a^81/(81!)^2 / (1-a/82^2)`, with `a=(79/20)*tau_max`.

The coefficient product is at most one and all later term ratios are at
most `a/82^2`, which is explicitly checked below one. On physical source
angles `|P_l|<=1`, so the residual rational R1 factor costs at most
`norm1*epsilon`. The retained R1 polynomial has modulus at most
`R1max+epsilon`, justifying its multiplier on the omitted Legendre tail.

The denominator Pstar is enclosed using the exact full finite-spin
positive polynomial and its two energy-interval endpoints. The crossing
polynomial keeps k below 64, and a separate exact Bernstein check proves
every omitted coefficient with k at least 64 nonpositive. This makes the
truncated crossing contribution a valid upper bound.

The K2 source correction uses its signed moments through order twelve,
with an explicit norm2 geometric remainder. Its zeroth moment agrees with
the stored Hs. The subtraction/addition of Hs in the field and the explicit
objective term +1 therefore have the correct signs and normalization.

The interval domain for each listed spin ends at
`min(256,T/45)`. This covers all energies at least 49 for spins through
106. At spin 108 and above, the separately checked unbounded-tau theorem
must start at spin 108 or lower; starting it at 110 would leave a gap for
spin 108. The final packet must verify this seam and include every even
spin from 42 through 448, as well as every interval in each partition.

This review approves the mathematical enclosure construction. It does not
assert that an unfinished batch passed, replace the 400/480-bit interval
records, or by itself establish a bound on g2.
