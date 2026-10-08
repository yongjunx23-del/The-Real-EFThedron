# Scoped finite-spin UV sign certificate for blend90

For the exact witness `output/lower_g0_upper_g2_20261006/blend_0.9.npz`,
SHA256 `d37101e716392336c8f359dbb59b3d14d04394cb133c9832607b322c6507e6d3`,
the new interval calculation proves

    V_l(mu)/p_l(mu) < 0

for **every even l>=10000 and every mu>=4+l(l+1)/256**. The endpoint at
mu=infinity is included by continuous compactification. Equivalently,
tau=l(l+1)/(mu-4) belongs to[0,256]. This is a uniform finite-spin result,
not merely the sign of the limiting function.

It is only a component of a possible coupling certificate: spins below10000,
the complementary tau>256 region, the real-kernel support integral, and
all head integration errors remain to be handled for this witness.

## Parameters and interval coverage

Set T=l(l+1), eta=1/T, w=1/(mu-4)=tau*eta, c=4/3, and
Q=mu-c=(1+8w/3)/w. On each rational tau cell[a,b], the implementation uses

    tau in[a,b], eta in[0,1/(10000*10001)], w=tau*eta.

Treating these intervals independently only enlarges the domain. Source
endpoints and rational harmonic parameters are exact decimals; T for strip
transfer is exactly4/3; witness multipliers are their exact binary64 values.

## Legendre denominator and crossing tail

The exact positive expansion is

    P_l(1+2c tau/T)
      =sum_k (c tau)^k/(k!)^2 prod_{j=0}^{k-1}(1-j(j+1)/T).

At fixed tau the first64 factors increase with T. Hence the first64 terms,
evaluated at T=10000*10001 and tau=a, are a rigorous lower bound for P_l.
The upper bound is I_0(2sqrt(c b)), obtained by replacing every product by1.
Both are enclosed with Arb; their interval hull encloses the denominator
at every permitted spin and energy.

The exact crossing numerator is the same positive series weighted by
the polynomial F(k,1/Q) from `matrix_crossing_series.py`. The exact
shifted-power/Bernstein certificate proves F(k,z)<=0 for every k>=10 and
0<=z<=3/8. Therefore dropping all k>=64 gives an upper bound for the
crossing numerator. The first64 terms are evaluated directly with interval
eta and w. No spin interpolation or cutoff approximation is used here.

## Exact finite-energy forward kernels and derivatives

All factors of mu are canceled before evaluation. With a formal transfer
variable t, write A=2+t/2, q=1+(2+t/2)w, and

    omega=(A w+i sqrt(1+(4+t)w))^2/q^2,
    Q^3 F=8 A^2(1+8w/3)^3/q^3.

The first four Taylor coefficients of P_l(1+2t/(mu-4)) are

    1, tau, tau^2(1-2eta)/4,
    tau^3(1-2eta)(1-6eta)/36.

These are exact, not asymptotic. Ordinary truncated Taylor algebra with
Arb complex intervals supplies every derivative-harmonic contribution
through order3. The subtraction kernels are factored exactly as

    Q^3 Ka=-(3t-4)^2(8w+3)^2/
             [27(tw+4w+2)(3tw+4w+3)],
    Q^3 Kb=-2(3t-4)^2(8w+3)^2/
             [27(tw-4w-1)(3tw+4w+3)].

The Kb terms multiply P_l and are combined with the outside contribution.
All ordinary and rational harmonic rows are included. An independent
assembly comparison against the existing finite-energy code is saved in
`tail_blend90_finite_forward_replay.json`; its largest discrepancy was
5.6e-7 at a far-UV point, where the ordinary floating implementation loses
precision. That comparison is a diagnostic, not part of the interval proof.

## Source kernel errors

The limiting source transform

    S(tau)=integral W K1_infinity J0(2sqrt(-t tau))

is enclosed with the exact source moments and centered Taylor remainder
from `tail_certify_joint_limit.py`. To compare it with the finite source,
use

    K1_finite P_l(y)-K1_infinity J0
      =(K1_finite-K1_infinity)P_l(y)
         +K1_infinity(P_l(y)-J0).

Every source argument y is physical in this UV region, so |P_l(y)|<=1.
The explicit uniform Mehler error is derived in
`tail_uniform_error_route.md`. Its worst b=(-t)tau uses -t<=79/20 and
the minimum T=10000*10001.

For mu>=M>48, the inherited rational kernel estimate is

    |K1_finite/K1_infinity-1|
       <= (8/M+2153/M^2+14088/M^3)
              /[(1-48/M)(1-(20/3)/M)].

This is the expanded-numerator bound from
`research/tail_certification/certify_uv.py`; it uses only the same source
domains, not the old witness's multipliers. Here M=4+T_min/b, with b the
upper endpoint of the tau cell.

For the other source kernel,

    K2_finite=K2_infinity/[1-((t-c)/Q)^2].

Thus its absolute integrated error is bounded by norm2*r^2/(1-r^2), where
r=(79/20+c)w_max/(1+(8/3)w_max). The exact rational source norms supplied
by `matrix_source_norms.py` are stored in `matrix_blend90_source_norms.json`:

    norm1 <=1602.8737962161426...
    norm2 <=422.3040578451336...

The code uses the exact rational fields, not these rounded decimals.

## Numerical proof and limits of independence

`tail_finite_uv_region.py` adaptively covers[0,256] by rational dyadic
cells, accepting a cell only when its complete interval upper bound is
strictly negative. The256-bit run passed with1421 leaves and2841 cell
evaluations. A separate320-bit run reevaluates the whole domain with the
same implementation and also passes with1421 leaves and2841 evaluations.
Exact rational endpoint checks confirm each accepted partition is gapless.
Both runs use the exact source-moment and Bernstein
norm certificates; this is not a formally verified or independently
implemented proof checker.

A separate agent reviewed the inequalities and independently checked the
two factored Ka,Kb identities symbolically. It also reconstructed the eps1
coefficient estimates: the mu^-1 coefficient is2-3t/2<=8, the mu^-2
absolute coefficient is at most2124.445<2153 using s*u<=0, and the cubic
coefficient is bounded by2112*(20/3)+64/9<14088. It found no issue in the
Mehler estimate, source-domain restrictions, omitted crossing signs,
Legendre denominator bounds, or the K2 correction. This internal review
does not replace an external mathematical audit.

No negative lower bound on g0 follows from this file alone. In particular,
the unbounded-spin statement here has an explicit energy restriction,
and the restricted sign does not by itself bound the remaining dual cost.
