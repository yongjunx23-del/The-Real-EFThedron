# Review of the b=6 uniform source proof

Reviewed `tail/b6_uniform_uv.py` on 7 October 2026. The objective is exactly
(1,6), with the original150-row witness and an exact107/5 increment to X4.
This changes neither U nor the Roy/strip source data. The original X4
coefficient is negative, and the ratio of its modified to original value is
strictly between zero and one. Rescaling its exact crossing polynomial by
that ratio therefore gives the intended modified crossing contribution.

The calculation uses eta=1/[ell(ell+1)] and tau=ell(ell+1)/(mu-4).
The independent integer-spin expansion has coefficients

    A_k(ell)=binom(ell,k) binom(ell+k,k)/[ell(ell+1)]^k.

For k<=40 and ell>=42 these are positive. They equal the product recurrence
used in the bivariate polynomial. Each retained coefficient decreases with
eta and increases with tau after multiplication by tau^k. Hence the endpoint
prefixes used by `p_bounds` are lower/upper bounds on the retained positive
Legendre series. The omitted part is bounded by the infinite positive series
with coefficient1/(k!)^2 and the displayed geometric ratio at k=41.

Source corrections are exact polynomials in eta and tau, followed by the
reviewed N12 rational-kernel remainder. For the largest possible w=tau*eta,
w<=8/(42*43); every source denominator remains positive. The source angular
tail uses the absolute positive series with parameter79*tau/20. Its ratio
after the first omitted term is less than one on tau0..8. The separate
second source remainder has ratio317*w/[60*(1+8*w/3)]<1. Neither bound drops
an unsigned residual.

The crossing contribution retains k0..40. For every k>=41 its coefficient
is proved nonpositive on z=w/(1+8*w/3) in[0,3/8] using exact Bernstein bounds.
Actual finite-spin Legendre coefficients are nonnegative up to k=ell and
zero beyond ell. Dropping those crossing terms is therefore an upper bound,
including for ell42 and44; it does not apply the polynomial product beyond
the physical finite-spin termination.

The inherited `ContractedRegion` object is initialized with450 to construct
its unused parent Legendre arrays. The called `numerator(tau,eta,w)` uses
explicit eta and w and contains no lower-spin-cutoff-dependent estimate.
Thus this call does not silently impose ell>=450. The b6 class supplies its
own L42 denominator, source remainder and crossing truncation bounds.

`audit_b6_uniform_source.py` independently specializes the bivariate source
and crossing polynomials at ell42,44,80,450 and compares them by exact
rational polynomial equality with the binomial-coefficient construction.
It also checks that the exact full finite Legendre polynomial is enclosed
at tau0,1,2,4,8. All checks passed at512 bits. Those point checks are diagnostics;
the paired320/400-bit complete eta/tau partitions prove continuum negativity.

The tau0..8 proof covers every even ell>=42 on that joint domain. It must be
combined with the separately certified tau>=8 theorem and the compact
energy certificate before an all-energy support bound is claimed. The
source and denominator remainders remain valid at eta=0 by continuity.
No formula or domain omission was found in this review. This is internal
review with an independent algebra check, not external review or formal
verification of all implementation dependencies.
