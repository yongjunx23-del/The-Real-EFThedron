# Full unitarity Gram certificate for the finite head

The kernel is the exact rational combination

\[
(U,V)=(0,p)-\frac14\partial_tG_t\big|_0-\frac12H_2-\frac52X_4.
\]

Its integral against `(r,b)` is `g2`, because the other three rows are null equations. This note certifies its support for spins `0,2,...,18` over all physical energies. A separate all-spin tail certificate is required to get a bound on the complete amplitude.

## Positive matrix and support inequality

For the unitarity matrix

\[
W=\begin{pmatrix}b&r\\r&2-b\end{pmatrix}\succeq0,
\]

choose any function `h` satisfying

\[
\mathcal C=\begin{pmatrix}h-2V&-U\\-U&h\end{pmatrix}\succeq0.
\]

Then

\[
\frac12\operatorname{tr}(\mathcal CW)=h-Ur-Vb\ge0.
\]

The certificate uses this complete matrix condition. It does not replace `2-b` by `2`. In particular it retains the correlation discarded by the earlier small reference-Gram estimate.

## Compactification and rational coefficients

Use the velocity `beta=sqrt(1-4/s)` and `e=1-beta^2`. On a rational cell `[a,b]` choose

\[
h(\beta)=e^2H_{[a,b]},
\]

with `H` an explicitly recorded nonnegative rational constant. The kernel expressions are rational functions of `beta`. The program divides the matrix by `e^2`, performs the congruence `diag(beta^ell,1)`, and multiplies by a positive common polynomial denominator `d(beta)`. The resulting polynomial matrix is

\[
\widetilde{\mathcal C}(\beta)=
\begin{pmatrix}
H A(\beta)+B(\beta)&C(\beta)\\
C(\beta)&H D(\beta)
\end{pmatrix},
\]

where

\[
A=\beta^{2\ell}d,\quad
B=-2\beta^{2\ell}d\,V/e^2,\quad
C=-\beta^\ell d\,U/e^2,\quad D=d.
\]

All entries are exact rational polynomials. The denominator is proved positive by checking a positive endpoint value and the absence of roots on `[0,1]`. The congruence is invertible throughout the open physical interval; endpoint singularities have therefore been removed without cutting away any energy region.

## Algebraic certificate on every cell

Expand this polynomial matrix in the Bernstein basis of its recorded cell:

\[
\widetilde{\mathcal C}(\beta)
=\sum_{j=0}^{n}M_j\binom nj t^j(1-t)^{n-j},
\quad t=(\beta-a)/(b-a).
\]

Every coefficient matrix `M_j` is checked positive semidefinite in exact rational arithmetic. For a `2x2` coefficient this means both diagonal entries are nonnegative and its determinant is nonnegative. Since the scalar Bernstein weights are nonnegative throughout the cell, this proves the entire matrix inequality on a continuum of energies.

The numerical calculations only choose trial values of `H` and where to subdivide. Acceptance depends entirely on exact rational determinant and diagonal checks. There is no optimizer tolerance, sampled positivity assumption, or numerical quadrature error in an accepted cell.

The matrices may be arranged in one large block diagonal positive matrix. Their interpretation is a collection of local polynomial certificates, not one monomial Hankel matrix. This distinction matters: increasing their number refines the continuum matrix certificate; it does not add new physical null equations.

## Reference integral and resulting number

The exact measure is

\[
d\nu_\ell=\frac{4(2\ell+1)}{\pi(1-\beta^2)^2}\,d\beta.
\]

Hence a cell contributes exactly

\[
\int_a^b h\,d\nu_\ell
=\frac{4(2\ell+1)}\pi H_{[a,b]}(b-a).
\]

The classical exact inequality `333/106<pi` makes `4/pi<=424/333`, yielding a rational outward bound. Summing all 3775 gapless cells gives

\[
\boxed{\sum_{\ell=0,2,\ldots,18}\int h_\ell\,d\nu_\ell
\le21.297281501.}
\]

The unrounded exact rational sum is recorded in the JSON certificate. This is a finite-head result, not a standalone bound on `g2`. Combining it with the separately proved tail `51/8` gives the explicit algebraic inequality

\[
g_2\le27.672281501<27.673.
\]

## Verification artifacts

- `certificate_full_w_head.py` constructs the kernels and accepted matrix coefficients.
- `certificate_full_w_head.json` records their rational polynomial coefficients and every rational cell and reference value.
- `certificate_replay_head.py` uses a different homogeneous-Horner/projective transformation to rebuild every Bernstein coefficient matrix, check positivity, check all partitions, and reconstruct the exact objective.
- `certificate_full_w_head_replay.json` records the completed replay and source hash.

The mathematical assumptions remain those of the original manuscript: its amplitude class, normalized disk unitarity, valid dispersion/null relations and convergence. These arithmetic checks do not independently prove those physical assumptions.
