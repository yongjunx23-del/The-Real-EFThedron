# Head integration and replication workflow supplement

This supplement documents the heads of the current seven inequalities, not a new optimization. The certified witness is the 150-row file `output/angles_resummed/angle_000p000.npz`. The 266-row numerical candidate must not be substituted into these fixed-witness certificates.

## 1. Exact claim and acceptance rule

Use a rational objective normal `(1,b)`, not rounded sine/cosine values. Set `d=0` if `b<=0` and `d=4*b` if `b>0`. The modified fields are exactly

\[
U_b=U_0,\qquad V_b=V_0+b p+d X_4.
\]

The attempted slopes are `-1,-1/2,-1/4,-1/10,0,1/10,1/4,1/2`. The last fails the sufficient finite-spin and infinite-spin estimates and has no complete certificate. The correct passing list has seven entries. For each, every component must refer to the same witness hash and the same exact `(b,d)`.

Every retained bulk cell, threshold interval, logarithmic source-endpoint interval and ultraviolet endpoint must have a proved finite upper contribution. The final claim is an outward-rounded rational upper bound on their sum. An optimizer status, small residual, sampled sign, matching high-precision decimal, or a plot is not this acceptance test.

## 2. Partition and numerical settings for the heads

Use `beta=sqrt(1-4/mu)`, `e=1-beta^2`, and `x=e/4=1/mu`. For even spin ell,

\[
d\nu_\ell=\frac{4(2\ell+1)}{\pi e^2}d\beta,
\qquad
D_\ell=\frac{4(2\ell+1)}\pi\int_0^1 H(U_\ell/e^2,V_\ell/e^2)d\beta.
\]

The exact rational beta partitions are saved in the old `*_tight_component.json` files and reused by `head_family.py`. Read the endpoints from those files; reconstructing them from rounded displayed energies is incorrect. There are five separately bounded endpoint gaps per spin, covering beta=0 and the four distinct source endpoints mu=11.9,12,24,48. The final beta=1 endpoint is handled by factored ultraviolet expressions.

| Spin | Old bulk cells | Refined bulk cells | Threshold cap | Source-endpoint gap construction |
|---:|---:|---:|---|---|
| 0 | 5,487 | 21,948 | beta in `[0,2^-26]` | exact dyadic gaps, original `--gap-bits 26` |
| 2 | 4,654 | 18,616 | beta in `[0,2^-14]` | exact dyadic gaps, original `--gap-bits 30` |
| 4 | 4,361 | 17,444 | beta in `[0,2^-9]` | exact dyadic gaps, original `--gap-bits 30` |
| 6 | 3,258 | 13,032 | beta in `[0,2^-9]` | exact dyadic gaps, original `--gap-bits 30` |

Refinement divides each bulk interval into four equal rational subintervals; gaps are unchanged. The eight-direction family is evaluated on the original partition at 256 and320 bits. The refined b=0 calculation is evaluated on the fourfold partition at both precisions. This is 64 original directional components plus8 refined components, or72 component records.

Scalar and low-spin source expansions retain N=64 terms. Scalar UV factoring starts at mu=96. For spins2,4,6 the UV factoring starts at mu=128 and uses Nuv=32 terms. The endpoint formulas include mu=infinity through x=0; no arbitrary energy cutoff is substituted for infinity.

## 3. Exact directional additions and threshold transfer

With c=4/3, Q=mu-c, Z=Q/(mu-4), the additions to the scaled absorptive field are

\[
\frac{p_\ell}{e^2}=\frac{xP_\ell(Z)}{16(1-cx)^3},\qquad
\frac{X_{4,\ell}}{e^2}=\frac{x^3[2Z^2P_\ell''(Z)-3ZP_\ell'(Z)]}{16(1-cx)^5}.
\]

The same exact rational expressions enter the derivative-series evaluator. They must not be added only at quadrature nodes while leaving the derivative remainder unchanged.

For the low-spin threshold cap, let the old interval proof give `Vold/p <= v_old <0`. Recompute

\[
v_{new}=v_{old}+b+d\,\sup_{\mathrm{cap}}(X_4/p).
\]

For d>=0 this is a valid upper enclosure. Require `v_new<0`; otherwise reject the cap or subdivide/rederive it. The unchanged real numerator lets the old support cap be multiplied by `(-v_old)/(-v_new)`. This is implemented in `DirectionalLow.threshold_cap`; it is not an assumption that the original threshold check automatically survives a change of direction.

At spin2, `(X4/p)(4+)=-9/32`; the selected d=4b gives an added threshold constant `-b/8` for b>0. This is only a local design motivation. It does not replace checks on the rest of the domain or on other spins.

## 4. Bulk integration error, not just grid convergence

For scaled fields f=(Uhat,Vhat), interval series of order2 enclose their second derivatives throughout a cell. If its width is h, use convexity and

\[
H(U+\delta U,V+\delta V)\le H(U,V)+|\delta U|+2|\delta V|.
\]

A uniform source-series error eta in Vhat contributes explicitly:

\[
\int_a^b H(f)\,d\beta\le
\frac h2[H(f(a))+H(f(b))]
+\frac{h^3}{8}\left(\sup|\widehat U''|+2\sup|\widehat V''|\right)
+4h\eta.
\]

Multiply by `4(2ell+1)/pi` with directed interval arithmetic. The series coefficient of order2 is half the second derivative, so the code's curvature expression is `2*abs(U_series[2])+4*abs(V_series[2])`. Source truncation errors are uniform errors, not differentiable constant offsets. Low-spin code may use the minimum of this bound and an independently valid interval Riemann/coercivity bound.

Denominators containing zero or nonfinite interval results must lead to subdivision or an analytic endpoint bound. They must never be silently dropped. Midpoint values are used only to prioritize adaptive subdivision and are not lower-bound certificates.

## 5. Principal values, removable singularities and logarithmic gaps

The physical `mu=s` pole is treated by exact polynomial division plus real principal-value logarithms. Do not delete a nearby quadrature node. The crossed subtraction denominator `mu-8/3+t` cancels in the combined even-spin kernel; use the divided-difference identity or exact polynomial division before interval evaluation.

For an endpoint beta-root r within a gap of width h<1,

\[
\int_{\mathrm{gap}}|\log|\beta-r||d\beta\le h[1-\log(h/2)].
\]

Use the exact identity

\[
\log\left|\frac{\mu-A}{\mu-B}\right|
=\log(A/B)+\log|\beta-r_A|-\log|\beta-r_B|
+\log|\beta+r_A|-\log|\beta+r_B|,
\quad r_A=\sqrt{1-4/A}.
\]

Bound root positions by Arb intervals and prove they lie in the intended exact rational gap. Separate each logarithmic coefficient, enclose the regular terms, and use `H<=abs(U)+2abs(V)`. At source-window boundaries the real local term is enclosed together with zero. For analytic derivative-series branches, decide the active source branch using exact beta endpoints, not an inflated interval image of mu.

## 6. Combining heads and interpolating final inequalities

For each spin and direction use the larger of its256-bit and320-bit certified upper totals. For b<=0,d=0, positivity p>0 gives `H(U,V+b*p)<=H(U,V)`; therefore a refined b=0 head bound is also valid at any negative slope. The assembler takes the minimum of that refined cap and the direction's own cap. This operation is valid separately for each spin.

After adding the finite8..40 block and the infinite tail, preserve the direct witness total. Search all pairs of passing slopes b1<b<b2 and set

\[
w_1=(b_2-b)/(b_2-b_1),\qquad C_{mix}=w_1C_1+(1-w_1)C_2.
\]

Take the smaller valid cap. The current b=-1/4 andb=-1/10 values are improved using b=-1/2 andb=0 with exact weights1/2 and1/5 on the left, respectively. The b=1/10 value uses b=0 andb=1/4 with exact weights3/5 and2/5. The ledger records these weights. This is positive addition of valid inequalities; it is not interpolation of optimizing coefficients followed by an assumed cost.

The angle is `atan(b)` and the unit-normal support is `C/sqrt(1+b*b)`, for normals with positive g0 component. These conversions are for plots. They do not certify the opposite-facing halfplane or a complete360-degree curve.

## 7. What the new replication checks establish

`python research/verify_replication_records.py` is a read-only standard-library check of the supplied frozen archive. It checks the expected ledger and witness hashes, all102 manifest entries (84 distinct paths after slash normalization), all72 head partitions/sums/threshold signs, exact assembly and pairwise interpolation. It printed all seven expected inequalities on6October2026. It does not freshly evaluate the kernel enclosures or prove the analytic lemmas.

`python research/review_head_replay_records.py --out PATH` audits all36 producer/replay pairs, including refined spin0. It verifies the256/320bit labels, matching exact partitions, dependencies, sums, signs and diagnostic agreement below1e-12. The original campaign's inline general-head audit did not have a retained producer script; this new file supplies a reproducible replacement record audit. It uses the same interval kernel calculations after they have been run separately, and is not a second kernel implementation.

The new helper was run against the frozen inputs and wrote only `output/replication_docs_20261006/head_replay_records_audit.json`. It did not change the original certification ledger or any frozen dependency.

The original refined-head review script checked spins2,4,6 only. Refined spin0 nevertheless has both complete kernel runs and is consumed by the exact assembler; the new helper explicitly includes it in the record audit.

## 8. Honest verification boundary

The primary finite8..40 interval producers evaluated every saved cell. Their separate higher-precision kernel audit sampled130 spectral cells and97 real cells; its partition/merge and affine-sign checks were exhaustive. Head and tail replays reused the same kernels. Hashes identify files and detect changes; they do not prove the formulas. This is a conditional computer-assisted certificate, not a formal theorem-prover result, proof of optimality, or construction of a primal-feasible amplitude.

For an independent replication, derive the kernels and all analytic tail inequalities again, implement them independently, evaluate the exact witness below, and bound the entire integral. Identical adaptive partitions or optimizer trajectories are unnecessary. A verified total no larger than the published six-decimal cap reproduces that inequality. A larger conservative total establishes only the larger inequality.
