# Independent formula review of the real-norm branch

Reviewed `research/all_direction_certification_20261007/real220.py` and
`real220_finite.py` against the transfer-derivative definition and finite-spin
head implementation. This is a formula review, not a completed coupling
certificate or an independent reimplementation of every interval operation.

For the fixed 220-row witnesses, the formulas examined are consistent.

Let z=beta^2, e=1-z, delta=4z/e, T_l=l(l+1), and
B_k(T_l)=product_{j=0}^{k-1}(T_l-j(j+1))/(k!)^2. The normalization is
Uref=16e^4/[beta(1+z)^5]. For physical forward transfer t,

```
F_t Im(omega_t)/Uref
 = beta^2 (1+t/4)^3 sqrt(1+t/mu)
   (1+t/[2(mu-2)])^(-5).
```

The sine part of `omega^j*(1-omega)^(2p-2)`, divided by Im(omega), is exactly
the binomial sum of Chebyshev U polynomials constructed as Dp. The regular-G
subtraction is purely real, so it contributes no real-field U term. Multiplying
this by P_l(1+2t/delta)=sum B_k(T_l)t^k/delta^k gives the stored rational
coefficient fields:

```
U_l = Uref * sum_{k=0}^3 w_k(beta) B_k(T_l)/delta^k.
```

For the real norm, the measure is
dnu_l=4(2l+1)/(pi e^2) d beta. Since
p_l=P_l(1+8/[3delta])/(mu-4/3)^3,

```
dnu_l U_l^2/p_l
 = 65536/(27pi) * e^3(2+z)^3/(1+z)^10
   * (2l+1)/(z P_l) * [sum B_k zeta_k]^2 d beta,
zeta_k = w_k/delta^k.
```

This matches the finite-spin code. Its second-order centered Taylor enclosure
of each rational zeta_k, followed by the exact spin-polynomial contraction,
preserves cancellation in the value and first derivative. Combining second
derivative balls with positive B_k is conservative.

For K<=l, the positive Legendre expansion gives

```
1/(z P_l) <= z^(K-1) /
  sum_{k=0}^K B_k (e/3)^k z^(K-k).
```

Thus the K=20 expression is valid for every spin42..400. The endpoint z=0 is
handled by the factored power and z=1 has denominator1. Any interval failure
must still be resolved by subdivision; no sampled sign fills a failed cell.

For all even spins greater than L, the exact reference coefficient test proves
P_l >= (1+a*T_l/delta)^N. Also B_i B_j <=T_l^(i+j)/(i!j!)^2. With m=i+j,
integrate x over[ell-2,ell]. On this interval T_x<=T_ell and
T_ell/T_x <=R=(L+2)(L+3)/[L(L+1)]. The weight integral is
T_ell-T_(ell-2)=4ell-2. Therefore the discrete sum is bounded using

```
(2ell+1)/(4ell-2) <= (2L+5)/(4L+6),
sum_even (2ell+1) T_ell^m/(1+aT_ell/delta)^N
 <= [(2L+5)/(4L+6)] R^m
    integral_{T_L}^infinity T^m/(1+aT/delta)^N dT.
```

This independently confirms the improved factor in `real220.py`, which is
slightly smaller than the older matrix-Gram helper's factor. Substituting
y=1/(1+aT/delta) gives exactly the finite beta-polynomial integral in that
file. The extra factor delta converts the common prefactor65536e^3/z into
262144e^2, as used there. N=16 is safely above the convergence requirement
N>m+1 for m<=6.

The strip term uses the same existing physical-strip Bernstein/L2 lemma:

```
|U_strip| <= pi/(B-A) * taper * sqrt(L2)
 * min(1, sqrt(2delta/[pi(ell+1/2)sqrt(T(delta-T))])).
```

The finite-spin code retains the actual ell. The infinite-spin code uses the
largest factor at ell=L+2. Adding the resulting strip bound to the absolute
k=0 coefficient is a valid triangle bound. The two source windows are disjoint;
an interval overlapping both conservatively includes both bounds. Roy real
source terms vanish for the considered spins since J<=6; their g0 subtraction
mass remains in alpha.

Suggested defensive checks for future use: require even first and last spins,
first>=K, first>6, last>=first, and verify the known canonical T=4/3 strip
inventory. The current fixed220,42..400 and L=400 applications meet these
requirements. This review does not replace the separate V-tail sign or
positive-V support remainder proof, nor the finite head integration.
