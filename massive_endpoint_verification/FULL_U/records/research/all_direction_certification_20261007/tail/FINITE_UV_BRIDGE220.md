# Finite-spin UV bridge for the frozen robust 220-row witness

This proof concerns every even spin 42 through 448, energy mu >= 49, and
tau = ell(ell+1)/(mu-4) <= 256. The separate `large_tau_highspin220.py`
certificate at L=108 covers tau >= 256. Spins below 108 cannot reach 256
when mu >= 49. This joins the finite bridge to the already reviewed
all-spin L=450 ultraviolet and compact-energy certificates.

The witness is `output/competitive_20261007/robust/g2_robust_stage1.npz`,
SHA256 `933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92`.
No coefficient is changed. The desired coercivity is V/p < -1/1000.

## Exact normalized field

Write T=ell(ell+1), w=1/(mu-4)=tau/T, c=4/3,
q=1+2cw, z=w/q, and Pstar=P_ell(1+2cw). The exact positive polynomial is

Pstar(tau) = sum(k=0..ell) a_ell,k (c tau)^k,
a_ell,k = product(j=0..k-1)[1-j(j+1)/T]/(k!)^2.

The forward harmonic, rational and transfer-derivative numerator and
outside terms are evaluated by the reviewed `ContractedRegion.numerator`.
All finite-w dependence is retained. The objective adds exactly 1.
The crossing contribution uses exact F(k,z) from the saved X4..X8
multipliers. Terms k=0..min(ell,63) are retained, with their positive
Legendre weights; terms k>=64 can be discarded in an upper bound.
The code proves the latter statement by expanding F(y+64,z) in y and
checking that every coefficient has nonpositive Bernstein coefficients
on 0<=z<=3/8. All polynomial expressions are rational and are shifted
to the exact cell midpoint before directed Arb evaluation.

## Exact source corrections

For each Roy or strip source point (s,t), normalized K1 differs from its
joint-limit value K1inf=2(s-c)(c-u) by the exact factor

R1(w) = [(1+(2+t/2)w)(1+2cw)^2] /
        [(1+(4-s)w)(1+(s+t)w)(1+(t+c)w)].

Here u=4-s-t, and all source angles 1+2tw are physical for mu>=49.
Set denominator=1+d1*w+d2*w^2+d3*w^3. Its numerator coefficients are
1, 22/3+t/2, 160/9+8t/3, 128/9+32t/9. Define the first N=12 formal
quotient coefficients r_j by the finite polynomial recurrence. The exact
identity, which requires no power-series convergence assumption, is

R1 - sum(j=0..N-1) r_j*w^j =
 w^N [r_N + (-d2*r_(N-1)-d3*r_(N-2))*w - d3*r_(N-1)*w^2]/denominator.

All three remainder numerator polynomials have exact Bernstein absolute
bounds on the source rectangles. The denominator is bounded below by
(1-44w)(1-157w/60)>0. The third factor uses t>=-79/20, the middle factor
is >=1, and w<=1/45. The resulting positive epsilon_R1(w) is increasing
there, so evaluating it at the upper cell endpoint is safe.

The exact source moments are
M[j,k] = integral W K1inf r_j(s,t) (-t)^k.
They are computed with rational beta integrals, not numerical quadrature.
Roy uses s=4+(79/10)x, t=-(79/20)x(1-z), x,z in [0,1], and includes
the folded-angular factor 2 in W. Strips use s=A+(B-A)x and
t=-(4/3)(1-z), and their taper and Legendre source weights. The saved
T label rounded to 1.33333 is never used as the physical endpoint.

The finite source approximation contracted before interval evaluation is
S_ell(tau) = sum(j=0..11) sum(k=0..min(ell,80))
 (-1)^k a_ell,k M[j,k] tau^(j+k)/T^j.

The rational R1 truncation error is bounded by norm1*epsilon_R1 using
|P_ell(1+2tw)|<=1 at physical angles. Here norm1 bounds integral |W*K1inf|
for the *combined* Roy/strip density, as in the reviewed source norm cache.
This avoids multiplying the rational error by a badly conditioned
monomial triangle bound for the Legendre polynomial.

If ell>80, the additional finite Legendre truncation is bounded by

sum(k>=81) (tmax*tau)^k/(k!)^2
 <= (tmax*tau)^81/(81!)^2 / [1-tmax*tau/82^2], tmax=79/20.

The ratio is below 1 for tau<=256. Multiply this by
norm1*(R1max+epsilon_R1), where
R1max=(1+2w)(1+2cw)^2/[(1-44w)(1-157w/60)].
For ell<=80 this truncation error is exactly zero.

The K2 source is also treated with signed moments, rather than a loose
O(w^2) absolute error. Its first 12 even-z terms are retained; the
remainder is bounded by norm2*r^24/(1-r^2), r=(317/60)z<1.
The constant Hs already present in the forward outside term is removed
and replaced by this exact signed K2 source. This is the `+Hs-s2`
combination in the implementation.

## Interval proof and replay

`finite_uv_bridge220.py` evaluates a *one-sided upper envelope* of V/p.
The added source remainders are positive error bounds; its point output
must not be described as a two-sided enclosure of the exact V/p.
For each spin its exact rational tau partition covers
[0,min(256,T/45)] without gaps. A cell is accepted only if the Arb value
is finite and its directed upper endpoint is strictly below -1/1000.
The endpoint tau=0 compactifies mu=infinity by a rational limit.

The producer hashes the fixed witness, source moments and builder,
source norm/component packet, crossing builder and contracted kernels.
`--replay FILE --precision 480` checks every saved cell independently at
higher precision, and binds the reference certificate's SHA256. The
finite-source cache can also be rebuilt exactly; all moments and all
Bernstein remainder bounds must agree as rational numbers.

This bridge by itself says nothing about mu<49. Compact spins42..108 are
a separate finite proof, while compact spins>=110 are covered by the
reviewed sharp all-spin theorem. The proof also does not certify the
finite-spin head integral, which is assembled by the head worker.
