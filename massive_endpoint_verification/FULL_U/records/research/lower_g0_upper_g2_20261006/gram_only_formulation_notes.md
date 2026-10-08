# A formulation using only positive Gram matrices

This note distinguishes a useful algebraic reformulation from additional physical assumptions. It records no new numerical bound.

## 1. What the matrices encode

Write the phase-space-normalized partial wave as `f = r + i b`, so the two-to-two unitarity disk is

\[
W=\begin{pmatrix}b&r\\r&2-b\end{pmatrix}\succeq0.
\]

This is equivalent to `r^2+(b-1)^2 <= 1`. For a positive, finite reference measure on `0 <= x <= 1`, define

\[
I_k=\int x^k b\,d\nu,\quad R_k=\int x^k r\,d\nu,
\quad \nu_k=\int x^k\,d\nu,
\quad Y_k=\begin{pmatrix}I_k&R_k\\R_k&2\nu_k-I_k\end{pmatrix}.
\]

Examples of necessary algebraic inequalities are

\[
0\le I_0\le2\nu_0,\qquad
R_0^2\le I_0(2\nu_0-I_0),\qquad
I_1^2\le I_0I_2,
\]

\[
(2\nu_1-I_1)^2\le(2\nu_0-I_0)(2\nu_2-I_2).
\]

The degree-one block matrix keeps further correlations:

\[
\begin{pmatrix}
I_0&R_0&I_1&R_1\\
R_0&2\nu_0-I_0&R_1&2\nu_1-I_1\\
I_1&R_1&I_2&R_2\\
R_1&2\nu_1-I_1&R_2&2\nu_2-I_2
\end{pmatrix}\succeq0,
\qquad Y_1-Y_2\succeq0.
\]

The second condition inserts the nonnegative factor `x(1-x)`. Higher orders impose `[Y_{i+j}] >= 0` and `[Y_{i+j+1}-Y_{i+j+2}] >= 0`. General weights give `G_ab = sum_l integral psi_a^T W_l psi_b dnu_l >= 0`; monomial weights give these Hankel matrices.

These are consequences of the pointwise disk. A finite matrix hierarchy can approach or efficiently represent that disk-constrained moment problem; it cannot supply additional unitarity information beyond the exact continuum disk with the same linear identities.

## 2. The crucial extra linear equations

Use **one shared set of moments** in the coupling maps, the crossing equations and the harmonic/resummed-harmonic identities. The finite problem then has the form

\[
G_j(g_0,g_2,m)\succeq0\quad(j=1,\ldots,J),
\qquad A m+B(g_0,g_2)=c,
\]

where `m` contains unresolved spectral moments. This is a semidefinite program when the matrix entries and equations are affine. Null equations are essential: unitarity matrices by themselves do not impose crossing or connect otherwise independent subtraction data to the spectrum.

For massive four-dimensional scattering, the manuscript uses `x=2/(s-2)` and its reference moments in A5.8. Those known moments contain beta-function constants. Algebraic in the unknown variables does not imply rational coefficients; enclosing known constants outwards, or treating their exact expressions correctly, is a separate numerical-certification requirement.

Do not silently remove high spins or energy regions. A finite all-spin relaxation needs blocks/majorants for all omitted sectors. The unweighted sum of reference measures over infinitely many spins need not be finite. Every actual Gram weight must have a convergent reference integral, or be part of an explicitly controlled tail construction.

## 3. Obtaining the requested list in g0 and g2

There are two legitimate forms of a deliverable.

1. A list of PSD matrices and linear equations with auxiliary moments, followed by existential projection onto `(g0,g2)`.
2. A list of consequences `a_i g0+b_i g2 <= C_i`, each accompanied by a nonnegative matrix combination and null multipliers proving it. Their intersection is an explicit outer polygon in coupling space. Increasing the direction list resolves its shape; improving the matrix/kernel hierarchy tightens it.

The second is usually much more practical than eliminating every moment symbolically. Merely discarding auxiliary moments is invalid. A finite projected semidefinite region is semialgebraic over its known real coefficients, but a simple short conjunction of determinant inequalities only in `(g0,g2)` is not generally guaranteed. Quantifier elimination can produce large formulas with cases. PSD certificates avoid that blow-up. Computable exact certificates additionally need exact expressions or validated enclosures for the reference constants.

All matrices should share the same auxiliary moments before projection. Intersecting individually projected matrix regions can be weaker than projecting their simultaneous intersection: each separate projection may use a different, mutually inconsistent choice of moments. Keeping the shared variables during discovery preserves this information; the final certified directional inequalities can then eliminate them.

Keep the full lower-right block `2 nu - I` whenever feasible. Replacing it by `2 nu` drops absorptive correlations: even the smallest determinant changes from `R0^2 <= I0(2nu0-I0)` to the weaker `R0^2 <= 2nu0 I0`. Some simple Gram-reference enclosures make such relaxations for convenience. Their weakness is not an inherent limitation of Gram matrices.

Checking only the diagonal entries or selected 2-by-2 minors of a large matrix loses information. Exact PSD is equivalent to nonnegativity of **all principal minors**, whose enumeration is usually less convenient than retaining the matrix.

## 4. Direct relation to the current support-dual method

For a coupling direction plus valid null multipliers, let the integrand be `U r+V b`. Introduce a reference function `h` and require

\[
C=\begin{pmatrix}h-2V&-U\\-U&h\end{pmatrix}\succeq0.
\]

Equivalently,

\[
h\ge0,\qquad h-2V\ge0,\qquad h(h-2V)\ge U^2.
\]

The exact identity

\[
\tfrac12\operatorname{tr}(CW)=h-U r-V b
\]

therefore gives `a g0+b g2 <= sum_l integral h_l dnu_l`. The least possible pointwise `h` is `V+sqrt(U^2+V^2)`, precisely the support function used by the existing solver.

Thus an all-matrix approach is viable. It is not restricted to the very coarse small Gram ellipses. It can seek polynomial/rational `h`, clear positive denominators and encode sector positivity by finite sum-of-squares Gram matrices. A matrix-polynomial interval certificate can take the form

\[
C(x)=Z_d(x)^TQZ_d(x)+x(1-x)Z_{d-1}(x)^TRZ_{d-1}(x),
\qquad Q,R\succeq0,
\]

with `Z_d=(1,x,...,x^d)^T tensor I_2`. Polynomial coefficient matching is algebraic. On several variables one uses appropriate domain factors/sectors. Logs and other nonalgebraic kernel pieces require certified envelopes or additional exact treatment before declaring the certificate polynomial.

This reformulation can certify a strong discovered candidate. Matching an existing published finite dual basis is plausible; exceeding the exact optimum of the identical continuum problem is impossible. More effective crossing/dispersion identities and better finite approximation can beat a published *finite computation*.

## 5. What the supplied PDF already demonstrates

The forward family in A5.11-A5.15 gives

\[
-\frac{128}{27}g_2-\sqrt{2K_\infty g_2}
\le g_0\le
-\frac{64}{27}g_2+\sqrt{2K_\infty g_2},
\quad g_2\ge0,
\]

where `K_infinity = 233.7516883972...`. This is a useful illustrative Gram-derived inequality, but it leaves `g2` unbounded and `g0` unbounded below even at infinite matrix order. This family alone cannot solve the requested targets.

For any `t>0`, the inequality `sqrt(2 K g2) <= t g2+K/(2t)` turns the same family into explicit linear cuts

\[
g_0+(64/27-t)g_2\le K/(2t),
\qquad
-g_0-(128/27+t)g_2\le K/(2t).
\]

Intersecting all positive `t` recovers the two curved boundaries. For a certified rational list, choose positive rational `t` and a certified rational upper bound on `K`.

A5.21-A5.25 describes a substantially richer construction: 231 PSD matrices, maximum size 48, and 69 directional consequences. The quoted enclosure includes `g0 >= -26.40757956` and `g2 <= 4.36947880`. These are much weaker than the current fixed-t dual benchmarks. The missing coefficient arrays prevent treating that entire stated construction as already independently reproduced here.

The obstacle is not that Gram matrices are incapable of strong bounds. The particular forward-family/reference estimates lose information; independent auxiliary sectors further relax physical correlations. A stronger construction should retain correlations, include effective nonforward/resummed equations, and control all spins.

## 6. What the previous phrase about remainder bounds actually meant

For the existing lower-g0 witness, derivative harmonics make the real kernel a low-degree polynomial in `l(l+1)`. A 4-by-4 **known reference Gram matrix** bounds the complete sum of squared real kernels over infinitely many spins by a quadratic form in its four coefficient functions. This replaces an uncontrolled omitted spin sum by an explicit upper envelope. Together with a strictly negative imaginary kernel it upper-bounds the omitted disk-support cost.

That is mathematically useful: it addresses an infinite set of constraints with finitely many entries. However, the required energy integration and the full complementary-domain analysis have not been completed for the latest lower-g0 candidate. These component lemmas have therefore **not yet delivered a certified lower-g0 number or an independent numerical improvement**. That limitation must accompany any claim of usefulness.
