# Uniform finite-spin corrections for the left witness

The fixed coefficient file is `output/full_direction_20261007/left/angle180_final.npz`,
SHA256 `e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
Its objective is (-1,0). All caches in `tail/left/refined` are rebound to
that exact file; none of the upper-g2 coefficients or source contractions
are reused as left coefficients.

The finite UV source proof in `FINITE_UV_BRIDGE220.md` extends continuously
in eta=1/[ell(ell+1)]. Its retained Legendre coefficient is the exact polynomial

    a_k(eta) = product(j=0..k-1)(1-j(j+1)eta)/(k!)^2.

For integer ell>=L>=108 and k<=80, a_k is positive and decreases with eta
throughout 0<=eta<=1/[L(L+1)]. Put tau=ell(ell+1)/(mu-4), w=tau*eta.
Then tau<=256 implies mu>=4+L(L+1)/256>49. All source Legendre arguments
are physical, and the finite R1 identity/remainder from the earlier proof
applies unchanged.

`uniform_source_uv.py` forms the entire retained source as a bivariate
rational-coefficient polynomial:

    S(tau,eta) = sum(j=0..11) sum(k=0..80)
       (-1)^k M[j,k] tau^(j+k) eta^j a_k(eta).

For each eta power the tau polynomial is shifted to the exact cell midpoint
before Arb evaluation. Those interval coefficients are evaluated by Horner
in the eta interval. This preserves the large cancellations inside each
coefficient of eta; it replaces the earlier norm-times-Mehler error.

The crossing numerator is likewise expanded exactly. With q=1+(8/3)tau*eta
and F(k,z)=sum_j f_kj z^j, its retained common-q numerator is

    C(tau,eta)=sum(k=0..63) a_k(eta)(4*tau/3)^k
               sum_j f_kj (tau*eta)^j q^(D-j).

Its contribution is C/q^D. The exact Bernstein certificate of F(y+64,z)
on z in [0,3/8] proves that all discarded crossing terms have nonpositive
contribution for every integer ell>=L.

For the positive denominator Pstar, retain k<=80. On a rectangle
[tau_lo,tau_hi] x [eta_lo,eta_hi], monotonicity gives the lower bound from
(tau_lo,eta_hi) and the upper bound from (tau_hi,eta_lo) plus

    (4*tau_hi/3)^81/(81!)^2 / (1-(4*tau_hi/3)/82^2).

This is valid for every integer ell>=L; it does not assume an interpolated
noninteger Legendre polynomial exists. Positivity of the actual finite
Legendre series and a_k<=1/(k!)^2 suffice.

The source remainder is the same uniform R1 error using |P_ell|<=1, plus
the k>=81 physical-Legendre tail with tmax=79/20. The signed K2 expansion
retains12 terms and uses the geometric remainder with norm2. Forward
harmonic/rational/transfer-derivative kernels use the already reviewed
exact real contraction, retaining every finite w and eta factor. The
left objective adds zero (not the +1 of the upper-g2 proof).

The eta interval is divided into16 exact rational panels. For each panel,
`certify_uniform_uv.py` builds an adaptive, gapless rational tau partition
of [0,256]. Every accepted Arb value must be finite and have directed upper
strictly below -1/1000. Higher-precision replay uses the same closed cells.
The union of the16 panels covers every integer ell>=450, including the
unbounded spin limit eta=0. This is an actual finite-spin statement, not
only a joint asymptotic limit.

For tau>=256, `left_large_tau.py` reuses the reviewed positive-Legendre and
fifth-power-Jensen bound with all coefficients/norms replaced by the left
witness. The bound applies for every even ell>=108 and mu>=49. Its +1
objective offset conservatively overestimates the actual left b=0. The
400-bit directed value is below -3000; a second precision and explicit
scope assembly are required before a complete claim.

The exact specialization audit in
`output/competitive_20261007/left/uniform_source_independent_audit.json`
compares source and crossing polynomials with the separate finite-spin
implementation at ell450,452,600. All rational coefficients match;18
independent512-bit point enclosures and denominator intervals overlap.
This complements the algebraic proof and does not replace domain coverage.

Lower finite spins42..448 use `certify_left_finite_uv.py`. It records every
exact local directed upper bound. In addition to global negativity, the
producer asks for upper<=midpoint_upper/4; this extra refinement is only
for useful local margins. The final theorem depends on each stored upper
and its higher-precision replay. The real-kernel norm must be combined
using the weakest kappa among every overlapping spectral cell, separately
for each spin. A uniform minimum alone is too loose for this left witness.

None of these UV statements covers4<mu<49 or the finite-spin head integral.
Those are separate proof components and must be included in final assembly.
