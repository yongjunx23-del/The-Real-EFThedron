# Independent review of the all-spin Gram remainder

This reviews `tail_m011_order8.json` for the rational12-row witness
`null_support_h8_d4_dropcross_m0.11_highspin.json`. It proves a reference-norm
remainder and its contribution to a support bound; it is not by itself a
complete bound on g2.

Write \(\delta=s-4\), \(L_\ell=\ell(\ell+1)\), and

\[
U_\ell=U_0\,[a_0(x)+a_1(x)L_\ell/\delta],
\quad x=2/(s-2),\quad
U_0=128\sqrt{s(s-4)}/(s-2)^5.
\]

The witness fixes a0 and a1 exactly. Crossing has no real kernel and therefore
does not enter this reference norm.

For every even ell at least20, positive Legendre Taylor coefficients imply

\[
P_\ell\left(1+\frac8{3\delta}\right)
\ge\left(1+\frac{aL_\ell}{\delta}\right)^8,
\qquad a=\frac{339}{4096}.
\]

The exact verification compares each binomial coefficient
\(\binom8k a^k\) to

\[
\frac{(4/3)^k}{(k!)^2}
\prod_{j=0}^{k-1}\left(1-\frac{j(j+1)}{420}\right).
\]

These coefficients are lower bounds because \(L_\ell\ge420\), all factors
are positive, and each factor increases with L. No large-spin expansion is
truncated without a sign: all omitted Legendre Taylor terms are positive.

For each even ell at least20, use the interval
\([L_{\ell-2},L_\ell]\). These intervals tile \([342,\infty)\). On them,

\[
\frac{2\ell+1}{L_\ell-L_{\ell-2}}\le\frac{41}{78},
\qquad L_\ell\le\frac{70}{57}\,t.
\]

The respective differences reduce to positive rational multiples of
\((\ell-20)\) and \((\ell-20)(13\ell-7)\), checked symbolically in the
review script. Dropping only a negative cross term gives

\[
\begin{aligned}
&\sum_{\ell\ge20,\,\mathrm{even}}
\frac{(2\ell+1)[a_0+a_1L_\ell/\delta]^2}{P_\ell}
\\[-2pt]
&\quad\le\frac{41}{78}\int_{342}^{\infty}
\frac{a_0^2+2(a_0a_1)_+(70/57)t/\delta
+a_1^2(70/57)^2t^2/\delta^2}
{(1+at/\delta)^8}\,dt.
\end{aligned}
\]

Thus all omitted spins are covered by one elementary integral. With
\(z=\delta/(\delta+342a)\), its three moments are proportional to

\[
z^7/7,\qquad z^6/6-z^7/7,\qquad z^5/5-z^6/3+z^7/7.
\]

These are exact antiderivatives, not quadrature approximations.

After changing energy to \(\beta=\sqrt{1-4/s}\), all factors are rational
functions of beta and the known constant pi. The remaining integral is
enclosed by directed Arb range evaluation over914 gapless dyadic beta
intervals covering [0,1]. Each outward endpoint is converted directly to its
exact rational `fmpq`, multiplied by the exact interval width, and summed
exactly. There is no numerical integration error omitted between nodes.

The original192-bit certificate and the independent320-bit fresh upper sum
both prove

\[
\boxed{\sum_{\ell\ge20,\,\mathrm{even}}
\int\frac{U_\ell^2}{p_\ell}\,d\nu_\ell<\frac9{100}.}
\]

The independent upper sum is approximately0.08999524790601926. Some
individual Arb radius enclosures are not nested as mantissa precision is
increased; the two complete directed upper sums independently lie below
9/100. The review does not require nesting or rely on rounded decimal
endpoint strings.

The separate exact absorptive-sign proof establishes
\(V_{\mathrm{raw},\ell}\ge(11/10)p_\ell\) for all these spins and energies.
For the objective g2 minus the null combination, the absorptive coefficient
is therefore at most \(-p_\ell/10\). Its disk support is bounded by

\[
\sqrt{U_\ell^2+V_\ell^2}+V_\ell
\le\frac{U_\ell^2}{2|V_\ell|}
\le5\frac{U_\ell^2}{p_\ell}.
\]

Equivalently, choosing \(h=5U_\ell^2/p_\ell\) makes the support certificate
matrix \(\bigl(\begin{smallmatrix}h-2V&-U\\-U&h\end{smallmatrix}\bigr)\)
positive semidefinite. Hence the entire omitted-spin contribution is at most

\[
\boxed{\frac9{20}=0.45.}
\]

The retained spins0 through18 need their own complete energy certificates
before this can be added to a physical coupling inequality. The proof is
conditional on the manuscript's dispersion/null identities and amplitude
assumptions. It is an exact-arithmetic and directed-interval computer-assisted
certificate, not a formal proof-assistant derivation of those physical inputs.
