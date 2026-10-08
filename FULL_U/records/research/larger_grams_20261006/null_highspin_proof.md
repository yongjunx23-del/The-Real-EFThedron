# Larger derivative/harmonic nulls: exact high-spin sign certificate

This note proves only the absorptive sign component. It does not certify the finite-spin support integral or the all-spin real reference norm.

The baseline candidates are `null_support_h8_d4_dropcross_m0.11_highspin.json` and `null_support_h8_d4_dropcross_m0.26_highspin.json` in the output folder. They prove, respectively,

\[
V_\ell(s)\ge\tfrac{11}{10}p_\ell(s),\qquad
V_\ell(s)\ge\tfrac54p_\ell(s)
\]

for every `s>4` and every even `ell>=20`. Here `0=<U*r+V*b>` is the exact null combination with the rational `raw_null_weights` written in the JSON. These weights are the negatives of the discovery support multipliers, rounded to eight decimal places. The rounded weights, not the unrounded discovery file, are certified.

Both use 12 identities: `Dt1G_regular`, `H2,...,H8`, `Dt1H_j2,...,Dt1H_j4`, and `X4`. The same verifier now certifies three larger candidates with `H16/Dt12` or `H32/Dt16` in their filenames. Their definitive summary is below.

| Identity count | Harmonics / derivative harmonics | Certified high-spin raw V/p floor | Exact proof cells | Independently audited support through spin 40 |
|---:|---|---:|---:|---:|
|12|H2–H8 / Dt1H2–Dt1H4|11/10|17|11.30362421844939|
|12|H2–H8 / Dt1H2–Dt1H4|5/4|20|12.17315577674115|
|28|H2–H16 / Dt1H2–Dt1H12|11/10|39|8.98466374629704|
|28|H2–H16 / Dt1H2–Dt1H12|5/4|40|9.57658450245199|
|48|H2–H32 / Dt1H2–Dt1H16|11/10|38|8.94582446626099|

The last column is **not a coupling bound**. It uses ordinary quadrature at 180 nodes per energy panel. The exact continuum sign is certified, while the finite-spin integration error and the omitted all-spin support still require independent certificates. The 28-row larger-margin candidate is useful because its tail majorant has coefficient 2 rather than 5.

## Reduction

Put `x=2/(s-2)`, `delta=s-4`, `Q=s-4/3`, `Z=Q/delta`, `L=ell*(ell+1)`, and `w=L*(Z-1)/2=4L/(3delta)`. The candidate absorptive kernel has the form

\[
\frac{V_\ell}{p_\ell}=C(x)+\frac{C_0(x)+C_1(x)w}{P_\ell(Z)}
+\kappa\frac{2Z^2P_\ell''(Z)-3ZP_\ell'(Z)}{Q^2P_\ell(Z)},
\qquad\kappa\ge0.
\]

All three rational functions are reconstructed symbolically from the physical first-transfer and harmonic kernels and written in the certificate. In particular

\[
C=\lambda_{Dt1G}Q^3[(s-8/3)^{-2}-s^{-2}].
\]

The script discards the last, nonnegative crossing term. Its sign follows directly from the positive Taylor expansion at `Z=1`. Write `z=(Z-1)/2` and `B_j=prod_{k<j}(L-k(k+1))/(j!)^2`. The first three coefficients of `2Z²P''-3ZP'` as a polynomial in `z` are

\[
\frac{L(L-8)}4,\quad
\frac{L(L^2-5L-30)}{12},\quad
\frac{L(L-2)(L^2+2L-96)}{96}.
\]

They are positive for `L>=420`. Every later coefficient is

\[
\tfrac12(k+2)(k+1)B_{k+2}
+\tfrac12(k+1)(4k-3)B_{k+1}
+k(2k-5)B_k\ge0,\quad k\ge3.
\]

There is no large-spin truncation in this step.

## Positive polynomial reference

For all `ell>=20`, the first thirteen Taylor terms obey

\[
P_\ell(Z)\ge\sum_{j=0}^{12}a_jw^j,
\qquad
a_j=\frac1{(j!)^2}\prod_{k<j}\left(1-\frac{k(k+1)}{420}\right)>0.
\]

On each exact rational interval `[lo,hi]` in `x`, rational Bernstein enclosures give `C0>=c0`, `C1>=c1`, and `C>=cc>cap`. Also `w>=(280*lo)/(1-lo)=w_min`. Therefore

\[
(V/p-cap)P\ge(cc-cap)\sum_{j=0}^{12}a_jw^j+c_0+c_1w.
\]

The latter degree-12 polynomial is certified nonnegative on `[w_min,infinity)` by the exact substitution `w=w_min+t/(1-t)`, multiplication by `(1-t)^12`, and nonnegative Bernstein coefficients on a rational partition of `[0,1]`.

All retained comparisons use exact `Fraction`/SymPy rational arithmetic. Floating point samples only reject prospective cells; they cannot accept one. The accepted `x` cells exactly cover `[0,1]` without gaps, as do all inner compactified `w` covers. The certificates require 17 and 20 `x` cells, respectively.

## Physical assumptions and present limitation

The proof is conditional on the manuscript's first-transfer derivative identities, requiring the local uniform convergence that justifies differentiating the dispersion relation. `X4` uses the manuscript's crossing normalization. The certification is exact arithmetic and analytic inequalities, not a formal proof-assistant development.

For an upper `g2` bound the support witness is `U_support=-U`, `V_support=p-V`. Hence `V_support<=-(cap-1)*p` on every omitted spin and the full unitarity matrix gives

\[
V_{support}+\sqrt{U_{support}^2+V_{support}^2}
\le \frac{U^2}{2(cap-1)p}.
\]

A separate all-spin reference integral and certified finite-spin support cost are necessary before reporting a coupling bound.

## Discovery results, not certified bounds

The full-support discovery objective through spin 40 is approximately 11.30250 for the original `cap=11/10` candidate and 12.17277 for the original `cap=5/4` candidate. The discovery code omits nodes with very large `p`, so these values are only candidate-selection diagnostics. The table's independent audit includes all quadrature nodes and covers threshold behavior with separate probes, but still does not certify integration error. A less restricted 28-row witness reaches approximately 7.424, but its high-spin continuum sign has not been certified here.

The independent `null_highspin_independent_review.py` checks every stored outer rational-function enclosure, every Taylor reference coefficient, the complete compactified polynomial, every inner exact Bernstein coefficient, and gapless coverage at both levels. It also compares the reconstructed symbolic kernels to the original physical kernel code at 18 points. The relative discrepancy through energy 1000 is below `3.5e-12`; at energy one million the ordinary floating derivative kernels lose precision through subtraction, producing differences up to `5.2e-8`. The exact certificates do not use those floating evaluations.

An earlier informal assertion that `C(infinity)<1` by itself invalidates a witness was incorrect: the `X4` contribution survives a joint large-spin/large-energy limit through `P''~ell^4`. It only obstructs the particular proof that discards this positive crossing contribution.
