# Independent review of `threshold_certificate.py`

Reviewed on 6 October 2026. This is a derivation/code review, not an independent complete reimplementation or a proof-assistant check. The scope is the limits of `V_l/p_l` as spectral energy approaches 4 from above for the listed finite even spins.

The formulas and polynomial quadrature error accounting in the reviewed version are consistent with the existing exact kernel definitions.

For `c=4/3`, the threshold value `Q=8/3` supplies the factor `Q^3=512/27`. For fixed even spin at least two,

\[
\frac{P_\ell(1+2t/(\mu-4))}{P_\ell(1+2c/(\mu-4))}
\longrightarrow (t/c)^\ell.
\]

The two subtraction terms combine before integration. Writing `u=4-s-t`, their threshold kernel is

\[
Q^3\left[
(t/c)^\ell\left(\frac1{4-s}+\frac1{4-u}-\frac1{2c}\right)
+\frac{1-(t/c)^\ell}{c+t}+\frac1{4-t}-\frac1c
\right].
\]

For even spin the divided difference is a finite alternating polynomial, exactly the `dd` recurrence in the code. This removes the apparent subtraction pole at `t=-c` algebraically. In the Roy source domain, `t` contains `s-4`; consequently `(t/c)^ell/(4-s)` is also a polynomial for every spin being certified. The principal-value endpoint limit can be justified by subtracting the source polynomial: its potential logarithmic term is multiplied by a positive power of `mu-4` and vanishes. The threshold formula does not introduce a new pole integral.

The forward and rational rows vanish after division by `p_l` for these spins. The threshold-regular harmonic derivative terms vanish as well. The surviving `DtG_regular` term is the negative derivative of the subtraction expression times `Q^3`; the reviewed code contains this factor. The crossing limit is exactly `F(ell,3/8)` from the independently derived positive Legendre-series identity.

The geometric approximant `(1-(1-den/center)^N)/den` is mathematically a polynomial of degree `N-1` in `den`, despite its rational evaluation form. Source denominator intervals do not contain zero in its uses. Its omitted tail is bounded by `r^N/[|center|(1-r)]`. The code integrates this truncated polynomial itself, so one remainder allowance is sufficient; a factor of two would apply to a different procedure that quadrature-evaluated the original rational function and bounded both its integral and quadrature remainders.

For Roy sources, `t` and `u` are bilinear in the two mapped source coordinates. The retained numerator therefore has degree at most

`max(nmax, jmax) + ell + N - 1`

in either coordinate. For strips, the energy taper adds degree two, and

`max(nmax+2, jmax) + ell + N - 1`

is a conservative coordinatewise bound. The chosen Gauss order satisfies `2*order-1 >= degree`, so tensor Gauss exactness applies to the truncated polynomial. The imported Gauss constructor encloses all roots and weights with Arb and explicitly proves completeness through disjoint sign brackets.

The reference absolute source mass `2*sum |lambda|` for Roy sources and `sum |lambda|` for strips are conservative. Legendre factors lie in `[-1,1]` on their source arguments, and the strip taper lies between zero and one. The degree-zero source masses used in the subtraction constant match the exact normalized source integrals, including the strip `n=2` term `-2/15`.

Two limitations should remain explicit. First, strict negativity of a limiting value does not provide an explicit finite-width threshold neighborhood. Second, finite spins `2..40` are not an all-spin certificate. The tool's current scope statement correctly preserves both limitations.

The strip parameter is interpreted as the intended exact `T=4/3`, as in the earlier analytic certificates, even when its name contains the usual decimal floating representation. A final witness ledger should state this parameter convention explicitly.
