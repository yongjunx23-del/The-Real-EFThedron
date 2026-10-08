# Independent review of the left-witness uniform UV extension

Witness: `output/full_direction_20261007/left/angle180_final.npz`, SHA256
`e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.

Reviewed `tail/uniform_source_uv.py` on 2026-10-07. This is a review of the
enclosure formula, not by itself a certificate of negativity on a domain.

Let T=ell(ell+1), eta=1/T, tau=T/(mu-4). The retained Legendre coefficient is

    a_k(eta) = product_{j=0}^{k-1}(1-j(j+1)eta)/(k!)^2.

The recurrence in the code gives this exact rational polynomial, including
the zero linear factor at j=0. For ell>=450 and k<=80, every factor is
positive throughout the enlarged interval 0<=eta<=1/(450*451).

The source contraction is the exact polynomial

    sum_{j=0}^{N-1} sum_{k=0}^{80}
      (-1)^k a_k(eta) M[j,k] eta^j tau^(j+k).

The bivariate matrix entries reproduce this expression before interval
evaluation. The crossing contraction similarly expands the common positive
denominator q^D, q=1+(8/3)eta*tau. Its coefficient with indices j,r acquires
eta^(j+r) tau^(k+j+r) binomial(D-j,r)(8/3)^r, exactly as implemented.
The omitted crossing coefficients k>=64 are nonpositive by the existing
exact Bernstein argument and may be dropped in an upper bound.

For P_ell(1+(8/3)/(mu-4)), the positive retained terms increase with tau and
decrease with eta. Thus evaluation at (tau_lo,eta_hi) is a valid lower bound,
and at (tau_hi,eta_lo) plus the positive factorial-series remainder is an
upper bound. The geometric ratio used for the first omitted term n=81 is
(4/3)tau_hi/82^2<1 throughout tau<=256. For an actual integer spin, terms
beyond ell vanish; replacing the remaining positive terms by the infinite
factorial series is conservative.

The source R1 correction error is bounded using physical |P_ell|<=1 and
the exact polynomial-division remainder from the reviewed finite model.
The extra Legendre truncation allowance bounds the truncated R1 polynomial
by Rmax+epsilon. Its denominators stay positive on this domain. The signed
K2 correction and its geometric remainder are unchanged. The left objective
has b=0, so the g2 objective term +1 in the original finite evaluator is
correctly absent. Interval division by the positive P enclosure remains
valid for either sign of the numerator.

The eta/tau rectangle enlarges the discrete-spin physical set. A successful
cell cover proves the desired bound for every actual integer even spin in
it, without requiring the extension to be the physical function at every
noninteger eta. At ell>=450 and tau<=256, mu>=4+450*451/256>49; the separate
unbounded-tau proof must cover the remaining mu>=49 region.

Independent executable audit:
`audit_uniform_source_left.py` and
`output/competitive_20261007/left/uniform_source_independent_audit.json`.
The bivariate source and crossing polynomials specialize **exactly** to the
previously reviewed finite-spin polynomials at ell=450,452,600. Eighteen
512-bit point enclosures agree with the finite evaluator after subtracting
its g2 objective +1; the P enclosures agree with the full positive Legendre
polynomial. These checks support implementation correctness; they do not
replace an interval cover in tau and eta.

No mathematical defect was identified in this review. This is an internal
independent review, not external peer review or a theorem-prover proof.
