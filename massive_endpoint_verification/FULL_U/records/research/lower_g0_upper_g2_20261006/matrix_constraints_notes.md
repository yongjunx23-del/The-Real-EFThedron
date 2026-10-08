# Hankel, Gram, and all-spin crossing work, 6 October 2026

This note concerns the identical massive scalar problem in four dimensions, in the normalization of the supplied manuscript. It introduces no additional analyticity assumptions. The results below are component lemmas, not a completed new lower bound for `g0`.

## What matrices add to the existing solver

At each spectral energy and even spin, write the normalized partial wave as `(r,b)`. The disk condition is exactly

\[
W=\begin{pmatrix}b&r\\r&2-b\end{pmatrix}\succeq0.
\]

For any finite positive reference measure, its block moments `Y_k=integral x^k W dnu` obey

\[
[Y_{i+j}]_{i,j=0}^d\succeq0,\qquad
[Y_{i+j+1}-Y_{i+j+2}]_{i,j=0}^{d-1}\succeq0.
\]

These are the Hankel and interval-localizing conditions. General weights give the Gram matrix `G_ab=integral psi_a^T W psi_b dnu`, also positive semidefinite. This is the construction in manuscript A1.1 and A5.3.

For a chosen direction and null multipliers, let the dual integrand be `U r+V b`. A scalar reference function `h` gives the equivalent matrix certificate

\[
C=\begin{pmatrix}h-2V&-U\\-U&h\end{pmatrix}\succeq0,
\qquad \tfrac12\operatorname{tr}(CW)=h-Ur-Vb\ge0.
\]

Its exact pointwise minimum is `h=V+sqrt(U^2+V^2)`, precisely the disk support already used in the current solver. Therefore adding integrated Hankel or Gram constraints formed from the same spectrum and the same nulls cannot beat the exact continuum disk-support optimum. A finite matrix relaxation is usually weaker. Matrices remain useful for:

1. Including additional valid harmonic, derivative, crossing, or resummed null equations efficiently.
2. Proving positivity on continuum cells after polynomial/rational transformations.
3. Bounding the infinite-spin contribution without discarding it.

Adding a PSD condition that was already implied by the disk is not a new physical assumption and should not be advertised as new independent information.

## Exact higher-crossing series certificate

Set

\[
c=4/3,\quad\delta=\mu-4,\quad Q=\mu-c,\quad
T=\ell(\ell+1),\quad x=c/\delta,\quad z=1/Q,
\]
\[
B_k(\ell)=\frac{\prod_{i=0}^{k-1}[T-i(i+1)]}{(k!)^2},\qquad
P_\ell(1+2x)=\sum_{k=0}^{\ell}B_k(\ell)x^k.
\]

All the weights `B_k x^k` are nonnegative. In the existing crossing catalogue, `beta_Dj` multiplies `2^j P_l^(j)/[j! delta^j Q^(D+1-j)]`. Differentiating the positive Legendre series gives the exact identity

\[
\frac{\sum_D\lambda_D X_D}{p_\ell}
=\frac{\sum_{k=0}^{\ell}B_kx^k F(k,z)}{P_\ell(1+2x)},
\]
\[
F(k,z)=\sum_D\lambda_D\sum_j
\beta_{Dj}\binom{k}{j}c^{-j}z^{D-2-j}.
\]

Repeated degrees use their separate catalogue rows. For the existing seven rows through degree eight, `F` has degree six in `k` and degree five in `z`. Its coefficients are exact rationals when the saved binary64 multipliers are interpreted as exact binary rationals.

For `output/angle_constraints/ming0_scaled_threshold_repaired.npz`, SHA256

`c7efae0020c9883d6949978ff924d11c3b828f470f10e24ff813f993e47cc06a`,

the new exact calculation proves

\[
F(k,z)\le-\eta,\qquad
\eta=\frac{380785068244075875273}{2882303761517117440}
\approx132.111359436886,
\]

for **every real `k>=10` and every `0<=z<=3/8`**. This covers all physical spectral energies, including the threshold and ultraviolet limiting values. It is not a sampled result or only an asymptotic expansion.

The proof substitutes `k=10+y`. Every coefficient of every nonnegative power of `y` has nonpositive Bernstein coefficients over the entire interval `z in [0,3/8]`. The constant term's largest Bernstein coefficient is `-eta`. All arithmetic in this sign test is exact rational arithmetic.

Consequently an infinite Legendre-series remainder can be removed with a proved sign:

\[
\frac{\sum_D\lambda_DX_D}{p_\ell}
\le-\eta+
\frac{\sum_{k=0}^{\min(9,\ell)}B_kx^k[F(k,z)+\eta]}{P_\ell(1+2x)}.
\]

This replaces a potentially difficult high-spin derivative ratio by a finite polynomial divided by a positive Legendre polynomial. In the joint limit `tau=T/delta`, the same formula has `z=0`, `B_k x^k -> (c tau)^k/(k!)^2`, and denominator `I_0(2 sqrt(c tau))`. It is useful for closing the large-`tau` part of an interval proof.

Files:

- `matrix_crossing_series.py`: derives and verifies the exact coefficient certificate.
- `../../output/lower_g0_upper_g2_20261006/matrix_crossing_series_certificate.json`: the polynomial, all rational Bernstein coefficients, witness hash, and diagnostic identity replay.

The independent diagnostic replay uses 80-digit positive Legendre sums and compares them to the existing SciPy Jacobi formula at 56 spin-energy points. The largest relative-scale discrepancy was `2.4201148511747848e-14`. This diagnostic is distinct from the exact rational positivity proof. The proof covers only the **combined crossing contribution**, not the Roy, derivative, forward, or strip parts of `V/p`.

## Real-kernel Gram tail including transfer derivatives

Above the largest local Roy spin, derivatives through order three make the forward part of the real kernel a polynomial of degree three in `T`:

\[
U_\ell(\mu)=\sum_{q=0}^3 a_q(\mu)T^q.
\]

This follows directly from the exact transfer Taylor expansion

`P_l(1+2t/delta)=sum_q B_q(l)(t/delta)^q`.

The energy coefficients can be obtained with the same exact truncated Taylor algebra as `forward_transfer_derivatives.py`; no finite differences are needed. Physical strip rows add their separate supported real kernel and require the existing strip majorants or a combined estimate.

Two exact Gram constructions are implemented in `matrix_real_tail_gram.py`.

### Compact-energy construction

Use `P_l >= B_n(l) (c/delta)^n`, with `n>=8` and all spins `l>L>=n`. The reference matrix

\[
R_{qr}=\sum_{\ell>L,\;\mathrm{even}}(2\ell+1)
\frac{B_q(\ell)B_r(\ell)}{B_n(\ell)}
\]

is positive semidefinite and has **exact rational entries**. To see this, decompose the rational function in `T` into simple fractions at `i(i+1)`, and use

\[
\frac{2\ell+1}{T-i(i+1)}
=\frac1{\ell-i}+\frac1{\ell+i+1}.
\]

The even-spin sums reduce to integer and half-integer digamma differences. Their common Euler constant and `log 2` terms cancel because the residues sum to zero; only rational harmonic numbers remain. The script stores the complete `4x4` matrix for `L=40,n=8`.

This construction is useful on bounded energy intervals. A single high Taylor term can produce a nonintegrable ultraviolet upper envelope even when the exact tail is integrable, so it must not be used blindly at infinity.

### Ultraviolet-integrable construction

Choose an exact positive rational `a` such that

\[
P_\ell(1+2c/\delta)\ge(1+aT/\delta)^n,
\qquad \ell>L,
\]

with `n=2d+2`. The code verifies every coefficient of this inequality exactly using `T>=(L+2)(L+3)`. For `L=40,d=3,n=8`, its rational choice is `a=347/4096`.

Define

\[
R=\frac{(L+2)(L+3)}{L(L+1)},\quad
D=\frac{2L+5}{2(2L+1)},\quad
y=\frac1{1+aL(L+1)/\delta}.
\]

An elementary comparison of each spin term with the interval `[l-2,l]` gives

\[
\sum_{\ell>L,\;\mathrm{even}}(2\ell+1)\frac{U_\ell^2}{P_\ell}
\le |a|^T G(\delta)|a|,
\]
\[
G_{qr}=DR^{q+r}\left(\frac\delta a\right)^{q+r+1}
\int_0^y h^{n-q-r-2}(1-h)^{q+r}\,dh.
\]

Here `|a|` denotes the vector of absolute **kernel coefficients**; the scalar rational `a=347/4096` in the denominator is separate. The beta integral is an elementary polynomial with exact rational coefficients. The matrix is a Gram moment matrix, and every entry is nonnegative.

Unlike the single-Taylor-term estimate, its high-energy scaling is the correct `delta^(q+r+1)`. The desired real-tail norm is obtained by integrating this upper envelope with `Q^3 dmu/(2 pi rho)`. This still requires validated energy integration and validated energy-coefficient bounds. Those integrations have **not** been completed in this file.

The script's example comparisons with finite spin sums check the implementation but do not certify any witness integral. These envelopes are conservative. They provide an implementable path to an all-spin derivative-tail certificate; they are not themselves a new numerical `g0` bound.

## Remaining obligations for a lower-g0 certificate

The existing approximately `7.78278` dual cost is not made rigorous by these component lemmas. One still needs:

1. Directed enclosures of the complete finite-spin support integrals.
2. Bounds for all remaining spins over all energy ranges, including the joint large-spin and high-energy regime.
3. The full combined absorptive kernel, not crossing alone, controlled with finite-spin errors.
4. A validated real-tail integral, including the derivative and strip terms.
5. All component costs summed with outward rounding for exactly the same fixed witness.

The new crossing certificate is already exact and complete in its stated domain. The real-tail formulas are analytic component bounds with executable exact rational reference matrices; they have not yet been assembled into a completed coupling certificate.

## Whole-domain Roy source guards added to the search

The enriched search developed a negative Roy source-weight patch between its sampled guard points. A sampled source margin was therefore insufficient. The new `matrix_roy_guards.py` implements constraints over the complete source domain.

Write `y=s-4`, `y_min=8/3`, `y_max=79/10`, and use coordinates

\[
y=y_{\min}+(y_{\max}-y_{\min})a,\qquad
z=q\frac{y-y_{\min}}{y},\qquad 0\le a,q\le1.
\]

For the source basis through `J=6,n=14`, define the polynomial

\[
P(a,q)=\left(\frac y{y_{\max}}\right)^6
\sum_{J,n}\lambda_{Jn}P_n(2y/y_{\max}-1)P_J(z).
\]

It has bidegree `(20,6)`. On a fixed partition of 16 energy cells and 4 angular cells, convert `P` to Bernstein form exactly. Let `R=(y/y_max)^6`, whose Bernstein coefficients are strictly positive. For every cell and coefficient index impose

`Bernstein(P) >= margin * Bernstein(R)`.

These are **9,408 linear inequalities** in the source multipliers. They prove `W>=margin` throughout the dangerous source domain, including between all numerical sample points. The search matrix divides by the corresponding positive coefficient of `R`, making each row `A lambda >= margin`. Search arithmetic is floating point; the postcheck reconstructs the polynomial and every subdivision with exact rational arithmetic using the final binary64 multipliers. Thus the actual certificate does not rely on floating feasibility tolerance.

The first projected witness, `lower_bernstein_projected.npz`, has a training head cost about `7.73954078623` and passes the exact whole-domain source check with

`W >= 0.00010001010003`.

Its SHA256 is `b87b1be5694f05a3c0be6807935933a82d2d020331b968e30e70ee7fffd78cf5`. It also passes the exact crossing remainder certificate with `eta=117.12391616976231`. The head cost is numerical and excludes the remaining obligations listed above; the two source/crossing component certificates do not turn that number into a complete certified coupling bound.

Relevant files are `matrix_roy_guard16x4.npz` for numerical search, `lower_bernstein_projected_roy.json` for the exact source certificate, and `matrix_crossing_bernstein_projected_certificate.json` for the exact crossing certificate, all under `output/lower_g0_upper_g2_20261006`.
