# Larger shared Gram matrices: a numerical experiment

This report establishes a computational mechanism, not a new continuum bound.
Every number below is an optimum of a finite quadrature and finite spin model.
They must not be compared to a published rigorous bound as if they were new
physical exclusions. The current certified physical bounds are recorded
separately by the main research task.

## The essential distinction

Increasing the old forward matrix beyond order128 can only approach that
family's limiting curve. That family has already discarded absorptive
correlations and remains unbounded in the requested lower-g0/upper-g2
directions. The experiment below changes what the matrix retains.

The physical disk matrix is

\[
W_\ell(s)=\begin{pmatrix}b_\ell&r_\ell\\r_\ell&2-b_\ell\end{pmatrix}\succeq0.
\]

For real basis functions \(\phi_i\), retain the complete shared block

\[
G=\begin{pmatrix}B&R\\R&2T-D\end{pmatrix}\succeq0,
\]

where all of \(B,R,D\) are moments of the **same** absorptive and real
variables. Replacing the lower block by \(2T\) drops the disk's finite radius
information. Adding more reference-only rows need not recover it.

## Exactly what was solved

Let quadrature node index \(n\) include an even spin and energy, and let
\(w_n>0\) be the common spectral quadrature weight. Define

\[
q_n=w_np_n,\quad a_n=\min(q_n,1),\quad
\gamma_n=a_n/q_n,\quad
t_n=q_nb_n/a_n,\quad z_n=\sqrt{q_n/a_n}\,r_n.
\]

These are positive diagonal congruences of the original disk:

\[
Y_n=\begin{pmatrix}t_n&z_n\\z_n&2-\gamma_nt_n\end{pmatrix}\succeq0.
\]

The complete nodewise condition is equivalently
\(z_n^2+\gamma_nt_n^2\le2t_n\). The moment hierarchy instead starts from
the weaker necessary inequalities

\[
t_n\ge0,\qquad \gamma_nt_n\le2,\qquad z_n^2\le2t_n
\]

and adds, for each spin and each chosen energy sector,

\[
\boxed{
\begin{pmatrix}
\Phi^T\operatorname{diag}(t)\Phi&
\Phi^T\operatorname{diag}(z)\Phi\\
\Phi^T\operatorname{diag}(z)\Phi&
2\Phi^T\Phi-\Phi^T\operatorname{diag}(\gamma t)\Phi
\end{pmatrix}\succeq0.
}
\]

Here \(\Phi\) is a QR-orthogonalized Chebyshev evaluation basis through
degree \(d\). Each block has dimension \(2(d+1)\). Orthogonalization changes
conditioning, not the spanned test functions. The coupling and null equations
use the **same** arrays \(t,z\), so no moment is independently reoptimized in
different projected families.

The linear coupling maps are

\[
g_2=\sum_na_nt_n,
\qquad
g_0=\sum_n\left[
u_{0n}\sqrt{w_n/p_n}\sqrt{a_n}\,z_n
+(v_{0n}/p_n)a_nt_n\right].
\]

Every null row has the identical linear map with \((u_0,v_0)\) replaced by
its own kernel pair. A singular value decomposition equilibrates those
equations; singular directions below relative threshold1e-10 are discarded.
Discarding them relaxes the discrete problem, but does not control its
continuum approximation error.

If \(\Phi\) has full row rank on its sector, congruence by its inverse makes
the block condition equivalent to all nodewise \(Y_n\succeq0\) there. Thus the
finite hierarchy has an exact, checkable endpoint: the full disk on that same
grid. It does not merely add arbitrarily many cuts from one unchanged weak
ellipse.

## Controlled matrix ablation

Keep exactly the same36 energy nodes, spins0,2,4,6, and10 nulls:
H2 throughH8, Dt1G_regular, Dt1H_j2, andX4. Velocity quadrature covers the
full energy range by open Gauss nodes on9 intervals. These are still only
quadrature evaluations, not rigorous integral enclosures.

| Largest block | Retain absorptive term | Replace lower block by reference only |
|---:|---:|---:|
|2|11.813116|11.814137|
|10|10.715278|—|
|18|8.690848|—|
|34|8.169355|11.814095|
|50|7.966201|—|
|64|7.949458|—|
|Full nodewise disk|7.949461|—|

Entries are the maximum g2 in this finite model. Several runs have solver
status `optimal_inaccurate`; differences in the last few digits are not
significant. The finite result at block64 agrees with the nodewise endpoint
to about3e-6. The reference-only hierarchy fails to improve because its
additional blocks follow from the already imposed nodewise relaxed
2x2 matrices.

The substantial change is therefore the information retained by the block,
not matrix dimension by itself.

## Localized matrices can be more efficient

On the exact same36-node model, split the Chebyshev test basis into the9
energy panels, retaining common variables and common null equations. There
are36 additional Gram blocks:9 panels times4 spins.

| Dimension of each block | Finite-model maximum g2 |
|---:|---:|
|2|10.673207|
|4|9.341846|
|6|8.089132|
|8|7.949462|

These36 blocks of size8 recover the finite nodewise disk value, while the
global construction needed four blocks of approximately size64. Every run
also retains the144 nodewise necessary2x2 reference cones and scalar bounds
described above; the table does not hide those auxiliary constraints.

This suggests using many interval/localizing matrices and effective weights,
alongside larger blocks when needed. One monolithic huge Hankel matrix is
not always the most efficient choice.

## Richer identities: promising but not yet stable enough for a bound

The richer model uses68 nulls:

-19 forward harmonics H2 throughH20;
-7 crossing kernels through degree8;
-3 regular transfer derivatives of G and21 derivative harmonics;
-18 Roy equations, with B=11.9, J=0,2,4, and polynomial orders0 through5.

All local blocks continue to share the same spectral variables.

| Nodes per energy panel | Maximum spin | Block dimension | Finite-model maximum g2 |
|---:|---:|---:|---:|
|24|16|6|1.120296|
|24|16|10|1.024600|
|32|20|6|1.337818|
|32|20|10|1.144810|
|48|20|10|1.144914|
|48|20|18|1.095844|
|48|24|10|1.236596|

The two block10 results at maximum spin20 agree well as quadrature order
increases from32 to48. However, increasing maximum spin20 to24 has a large
effect. Consequently, neither1.0246 nor1.0958 is an all-spin bound or evidence
of beating a published dual. These experiments directly show why an apparent
improvement from a small truncated computation can be misleading.

The rich runs mostly return `optimal_inaccurate`, with scaled null residuals
around1e-7. The pointwise rich solve encountered a solver failure. A stable
dual residual check and all omitted-sector bounds remain necessary.

## What would turn this into explicit certified physical inequalities

The useful output of discovery is a dual positive matrix combination and null
multipliers. They must be converted to polynomial or rational matrix
inequalities valid on continuous energy sectors and every omitted spin,
with positive denominators proved and all reference integrals enclosed.
This is precisely where continuum matrix positivity, interval localization,
Legendre remainder bounds, and exact or interval reference moments enter.

The missing continuum work cannot be replaced by increasing the quadrature
order, reporting a favorable solver status, or drawing the finite result
against a published curve. None of the finite numbers in this note are
exported as physical halfplanes in g0,g2.

## Reproduction and outputs

- `formulation_scaled_grams.py`: all well-scaled experiments above.
- `formulation_report.py`: CSV aggregation and matrix ablation plot.
- `../../output/larger_grams_20261006/formulation_scaled_grams.json`: global hierarchy.
- `../../output/larger_grams_20261006/formulation_scaled_grams_relaxed.json`: ablation.
- `../../output/larger_grams_20261006/formulation_local_basic.json`: localized hierarchy.
- `../../output/larger_grams_20261006/formulation_local_rich_*.json`: richer studies.
- `../../output/larger_grams_20261006/formulation_matrix_ablation.png`: comparison plot.
- `../../output/larger_grams_20261006/formulation_study.csv`: all logged results and statuses.

The earlier `formulation_shared_moments.py` was an intentionally finite
compact energy-window prototype. It spuriously drove amplitudes nearly to
zero when global nulls were imposed without outside-window contributions;
its numbers are not used above. This prototype is kept as a reproducible
failure case, not as evidence for a strong bound.
