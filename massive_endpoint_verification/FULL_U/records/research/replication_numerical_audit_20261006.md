# Replication audit: the 266-row numerical endpoint candidate

Audit date: 2026-10-06. Workspace: `C:/Users/laure/Documents/ChatGPT/Unitarity codex`. All relative paths below are relative to that workspace. This document describes the **uncertified 266-row numerical candidate** with measured partial cost **2.6988815404066395**. It does **not** describe the separate, certified 150-row witness. In particular, no certificate for the 150-row witness can be attached to the 266-row coefficients by changing a filename.

The authoritative final numerical witness is `output/directional_certification/numerical/bridge/bp0_stopped.npz`, SHA-256 `fab65bffd85740b1d6043124bd1f2a0324547e1c033f98274340ab0b92beb2e4`. Its arrays are `names`, `multipliers`, and `direction`; the direction is `(1,0)`. All 266 multipliers are nonzero. The main replication handoff includes the exact ordered coefficients as binary64 hexadecimal strings; section11 below also gives a lossless extraction/reconstruction recipe for the supplied NPZ. This removes any dependence on an unpublished spectral-kernel/data release: the spectral kernels are rebuilt by the supplied code.

## 1. What the number means

| Integrated domain | Numerical contribution | Record |
|---|---:|---|
| Every even spin 0..40, quadrature mapped over the entire energy ray 4<mu<infinity | 2.698881478706302 | `output/directional_certification/numerical/fine/bp0_stopped_head.json` |
| Every even spin 42..100, 4.02<mu<49 | 3.4472427435560786e-8 | `output/directional_certification/numerical/fine/bp0_stopped_source.json` |
| Every even spin 42..100, entire mu>49 ray | 2.722791006942304e-8 | `output/directional_certification/numerical/fine/bp0_stopped_uv.json` |
| Sum of those numerical contributions | **2.6988815404066395** | `output/directional_certification/numerical/endpoint_summary.json` |

The previous 266-row candidate's analogous measured sum was 2.698943676116444. The measured improvement is 6.213570980451166e-5. The newest candidate's fine spin<=40 evaluation exceeds its coarser training objective by 4.049988322485376e-7.

**No rigorous error bar is supplied for any of these quadratures. The domains 4<mu<4.02 at even spins42..100 and all energies at even spins>100 remain unmeasured.** Because the disk support is nonnegative, omitting such domains does not by itself give an upper bound on the complete dual cost. Thus the sum is a partial numerical dual cost, not a certified upper bound on physical g0. The small change under one grid refinement is a convergence diagnostic, not a certificate or a proved uncertainty interval. There is also no assertion of optimizer convergence or global optimality for the final checkpoint.

## 2. Physical assumptions, normalization, and dual objective

The mathematical identities assume the supplied manuscript's class of identical massive scalar amplitudes in four spacetime dimensions: m=1; s+t+u=4; crossing symmetry and real analyticity; no additional subthreshold poles in the assumed amplitude class; and the stated fixed-t analytic domain, twice-subtracted dispersion relations, growth/convergence assumptions, and legitimacy of the smearing/sum/integral exchanges. These assumptions must be retained when interpreting either a numerical witness or a later certificate. This computation does not independently prove the underlying dispersion theorem for a broader amplitude class.

Let f_l(mu)=r_l(mu)+i b_l(mu), rho(mu)=sqrt(1-4/mu), and

```
dnu_l = (2*l+1) dmu / (2*pi*rho(mu)),     l = 0,2,4,...,
Q = mu-4/3,     Z = Q/(mu-4),     p_l = P_l(Z)/Q^3 > 0.
```

The pointwise unitarity condition actually used is the **full disk inequality**

```
r_l^2 + (b_l-1)^2 <= 1.
```

It includes 0<=b_l<=2, but it is stronger than merely imposing that box. Elastic equality is not imposed. No restriction on the amount of inelasticity beyond the disk is added. No published numerical upper bound on g2 is used as a prior.

For a real pair (U,V), the pointwise disk support is

```
H(U,V) = V + sqrt(U^2+V^2).
```

When V<0, the implementation evaluates the identical expression `U*(U/(hypot(U,V)-V))` to avoid subtractive cancellation. For a direction `(a,b)`, the formal complete dual bound is

```
a*g0 + b*g2 <= sum_even_l integral_4^infinity H(U_l,V_l) dnu_l,
U = a*u0 + sum_i lambda_i*u_i,
V = a*v0_raw + b*p + sum_i lambda_i*v_i_raw.
```

Here every `(u_i,v_i)` is an exact null after the subtraction-coefficient conversion described below. The saved final direction is a=1,b=0. A direction `(1,b)` is not a unit vector: comparison with a unit-angle support function requires dividing the cost by sqrt(1+b^2).

The code generally assembles `V/p` first, sums all terms there, and only then multiplies by p. This avoids catastrophic cancellation between individually enormous threshold kernels. `DiskSupport` also computes analytic gradients and Hessians of the discretized support; it does not finite-difference them. Its numerical radius denominator is floored at1e-300 to avoid division by zero; this is an implementation regularization, not an interval proof. Sources: `research/directional_certification/numerical_dual_newton.py:16`, `:32`, `:36`; `research/angle_audit.py:68`, `:139`, `:146`.

## 3. Complete physical row inventory: exactly 266 rows

### 3.1 Coupling kernels, which are the objective rather than extra null rows

Define

```
F(mu) = 32/(mu-2)^3,
omega(mu) = (2/(mu-2) + i*sqrt(mu*(mu-4))/(mu-2))^2,
C_l(mu) = 8/((3*mu-8)*(3*mu-4)) * (1/(mu-2)-4*P_l(Z)/mu).
u0 = F*Im(omega),
v0_raw = F*Re(omega) + C_l,
g2 kernel = (0,p_l).
```

These are the normalization implemented in `research/existing_solver/massive_forward_kernels.py:31` and `:52`, identified there with manuscript equations A5.1 and A5.3. The common measure already includes `(2*l+1)/(2*pi*rho)`; do not add an extra spin factor inside the raw kernels.

### 3.2 Harmonic nulls: 39 rows

The saved names are **H2,H3,...,H40**, with complex kernel

```
q_Hj(mu) = F(mu)*omega(mu)^j,     (u,v_raw)=(Im q, Re q).
```

The actual normalized inventory uses these individual harmonics (`research/angle_audit.py:80`). Do not confuse the saved `Hj` names with the alternative UV-subtracted bank `Hj-(-1)^j H2` returned by `harmonic_bank` in `massive_forward_kernels.py:96`; those span the same space only when H2 is retained and coefficients are transformed correctly.

### 3.3 Rational resummations: 33 cosine plus18 sine rows

Each row is an exact rational resummation of harmonic nulls, not a truncation of its infinite series:

```
q = F*omega^2*(m(omega)-m(-1)),
m_cos(z) = (1-r*c*z)/(1-2*r*c*z+r^2*z^2),
m_sin(z) = r*z/(1-2*r*c*z+r^2*z^2).
```

All listed radii are less than1. The evaluator combines the difference into one rational expression containing an explicit `(omega+1)` factor, preserving UV cancellation (`research/existing_solver/massive_forward_kernels.py:116`).

| Kind | Radius r | Complete set of cosine parameters c | Count |
|---|---:|---|---:|
| cosine | .8 | -1,-.75,-.5,0,.5,.75,1 | 7 |
| cosine | .95 | -1,-.95,-.9,-.75,-.5,-.25,0,.25,.5,.75,.9,.95,1 | 13 |
| cosine | .98 | -1,-.95,-.9,-.75,-.5,-.25,0,.25,.5,.75,.9,.95,1 | 13 |
| sine | .95 | -.9,-.75,-.5,-.25,0,.25,.5,.75,.9 | 9 |
| sine | .98 | -.9,-.75,-.5,-.25,0,.25,.5,.75,.9 | 9 |

Literal row labels, including differences such as `-1.0` versus `-1`, are preserved in the supplied NPZ and the main handoff's exact coefficient appendix. They matter for byte-identical inventories and name-based warm starts. They do not indicate different numerical parameters when parsed.

### 3.4 Crossing null: one X4 row

The only standalone algebraic crossing row is

```
X4 = < [2*Z^2*P_l''(Z)-3*Z*P_l'(Z)]/Q^5 * b_l > = 0.
```

Its real-part kernel is zero. The active builder uses exact rational crossing-coefficient elimination and a normalized Jacobi-polynomial formula (`research/crossing_kernels.py:24`, `:69`); its first-row normalization agrees with `research/existing_solver/massive_forward_kernels.py:74`. The final multiplier is **-33.088397116313395**. X6, X8, higher-degree crossing rows, and their possible multiplicities are not part of this witness.

### 3.5 Roy nulls: 95 rows

Every combination

```
B = 11.9,
J in {0,2,4,6,8},
n in {0,1,...,18}
```

is present, with saved name `Roy_B11.9_J{J}_n{n}`. These are physical even-J projections, not the spectral-spin cutoff. The source variable s ranges over (4,B),

```
x = 2*(s-4)/(B-4)-1,
w_n(s) = P_n(x)/(B-4),
z in [0,1],
t = -(s-4)*(1-z)/2,
u_ext = 4-s-t.
```

Even angular parity folds the original angular integral to `2*integral_0^1 P_J(z) dz`; this factor2 is retained in the numerical weights. The raw Roy row obeys `<u_R*r+v_R*b> + c_R*g0 = 0`, where c_R=-2 only for J=n=0 and is zero otherwise. Its local raw real-part kernel is `2*pi*w_n(mu)/(2*J+1)` on 4<mu<B at the single spectral spin l=J, and zero at other spectral spins.

For clarity, define the divided difference `D_l(x,y)=(P_l(x)-P_l(y))/(x-y)` with its polynomial continuation at equality, and `y=1+2*t/(mu-4)`. Before smearing, the absorptive dispersion expression is

```
K_l(mu;s,t) = P_l(y)/(mu-s)
             + P_l(y)*(1/(mu-u_ext)-1/Q)
             + P_l(Z)*(1/(mu-t)-2/Q)
             - 2/(mu-4)*D_l(y,-Z).
```

The raw absorptive kernel is the **negative** smeared principal value of this expression. The `mu=s` pole is treated analytically: subtract `P_l(z)` before numerical integration, use the polynomial divided difference for the remainder, and restore the Legendre-Q principal value only when l=J. The added normalized-P_l integral is `4*Q_n(2*(mu-4)/(B-4)-1)/((2*J+1)*(B-4))`, then divided by P_l(Z) in the normalized builder. For mu>=B+1 the evaluator uses an algebraically factored nonsingular expression that exposes UV falloff.

These definitions, including the sign and local coefficient, are implemented at `research/roy_kernels.py:71`, `:102`, `:120`, `:133`, `:148`, `:155`. The real PV branch Q_n and cancellation-safe recurrence are at `:38` and `:17`. There are **no odd-J Roy rows** in the final witness.

### 3.6 Tapered fixed-transfer-strip nulls: 80 rows

Every combination

```
(A,B) in {(12,24),(24,48)},
T = 4/3,
j in {0,1,2,3,4},
n in {0,1,...,7},
taper power = 1
```

is present, with labels `Taper1_Strip_A{A}_B{B}_T1.33333_j{j}_n{n}`. **The rounded label T1.33333 is parsed back to the canonical T=4/3**, not used literally (`research/angle_audit.py:27`). The source variables range over A<s<B and -T<t<0. Set x=2*(s-A)/(B-A)-1 and eta=2*(t+T)/T-1. The source test is

```
w(s,t) = (1-x^2)*P_n(x)/(B-A) * P_j(eta)/T.
```

The raw strip identity is `<u_S*r+v_S*b> + c_S*g0 = 0`. Before tapering, c_S=-1 only for j=n=0. Its local raw real part is `pi*w_n(mu)*H_lj(mu)` on A<mu<B, with `H_lj(mu)=integral_-T^0 P_j(eta)*P_l(1+2*t/(mu-4)) dt/T`. Thus strip real kernels extend to **every even spectral spin**, an essential tail issue.

The absorptive part is minus the same K_l above smeared over the strip. The source energy integral's two poles are evaluated with real Legendre-Q formulas near the source window (`research/roy_kernels.py:171`, `:207`, `:218`, `:223`, `:242`). The taper is an exact finite Legendre combination: repeatedly multiply coefficients by `[2/3,0,-2/3]`, the expansion of 1-x^2; power1 requires untapered energy degrees through n=9. Sources: `research/angle_audit.py:35`, `:104`, `:128`; `research/transfer_kernels.py:9` and `:39`.

For the tapered rows, c_S is -2/3 at j=0,n=0, +2/15 at j=0,n=2, and zero otherwise. All local and absorptive terms are tapered by the same finite combination.

### 3.7 Subtraction coefficients are converted into true zero rows

The Roy/strip builders do not, by themselves, return a zero row: they return `<raw row>+c*g0=0`. The numerical inventory therefore adds `c*(u0,v0_raw)` to that raw row. In normalized output this reads `u_i=bank.u_i*sqrt(p)+c_i*u0` and `(v_i/p)=bank.v_i+c_i*(v0_raw/p)`. Omitting these terms changes the problem. Source: `research/angle_audit.py:122` through `:132`.

The endpoint's effective forward-g0 coefficient in a tail implementation is consequently

```
c0 = 1 - 2*lambda_Roy(0,0)
       - (2/3)*sum_windows lambda_Strip(j=0,n=0)
       + (2/15)*sum_windows lambda_Strip(j=0,n=2).
```

### 3.8 Explicit scope exclusions

No fractional finite-difference row, fractional derivative row, forward t-derivative row `Dt...`, odd-Roy row, higher crossing row than X4, extra fixed-energy angular sampling identity, external g2 prior, primal ansatz, or imposed elastic saturation enters this 266-row result. Some such code exists in the repository and some was investigated earlier; existence in the repository is not evidence that it contributed to this saved witness. The rows above and the saved ordered266-name inventory are exhaustive.

## 4. Additional optimizer restrictions: all of them

These restrictions constrain the **choice of dual multipliers**. They do not add assumptions about allowed physical amplitudes. They make the search numerically safer but their finite sampling is not a continuum proof.

1. **Physical high-p cells:** at every training node with p_l(mu)>10^4, require `V_l/p_l <= -eps`, with eps=2e-6. There is no claim that V is negative everywhere; spin zero has sampled positive V/p and its full disk support is integrated.
2. **Zero-weight threshold probes:** append spectral energy mu=4+1e-10 at every training head spin l=0,2,...,24 with quadrature weight zero. The same high-p selection decides which of these impose sign constraints. Zero weight does not suppress their constraint. The l=0 point has p<10^4 and is not forced negative by this rule. The separate bridge never receives this dangerously small delta probe, because raw high-spin Legendre values would overflow.
3. **325 Roy-source guards:** use25 equally spaced s values in [20/3,11.9]; at each s use13 equally spaced z values in `[0,1-(8/3)/(s-4)]`. Let `W(s,z)=sum_Jn lambda_Jn*P_n(2*(s-4)/7.9-1)*P_J(z)`. The artificial row vector is minus the coefficient vector of W and has base zero, enforcing W>=a small positive margin after normalization. At s=20/3 the13 angular points coincide; these duplicates are retained.
4. **One X4 sign guard:** an artificial row with coefficient1 for lambda_X4 and0 otherwise, base zero, effectively imposes lambda_X4<=-2e-6 (stronger than the advertised necessary nonpositive sign).
5. **88 additional spin/energy guards:** all even spins26..40 at mu in `{4+1e-6,4+.01,4.1,5,8,12.1,20,50,100,1000,1e5}`. Their base is the actual g0 normalized kernel plus b for direction `(1,b)`, and row coefficients are the actual normalized null kernels. These guards remain present in the bridge-trained run even though the same spins are also included in its integral objective.

For **each** artificial column, compute `scale=max(max_i abs(vector_i),1)`; divide the row vector and base by that scale; give it artificial p=10^5, real-part/base U=0, and cost weight0. Since 10^5>10^4, the standard high-p mask includes all414 artificial columns. Their strict normalized target means the unscaled absorptive inequality is `base+lambda.vector <= -eps*scale`. For the Roy guard it is `W>=eps*scale`; scale is1 for these bounded Legendre factors. The additional-spin scale can exceed1. This precise scaling matters if reproducing the search.

Sources: `research/directional_certification/numerical_direction_search.py:103`, `:123`, `:130`, `:131`, `:136`; high-p mask and target in `research/directional_certification/numerical_dual_newton.py:91`. The separate151x151 source screen (`numerical_direction_search.py:22`) is a diagnostic report, not325 additional constraints and not an exact polynomial sign certificate. The final stopped checkpoint has no newly emitted full optimizer screen record; its saved accepted-step feasibility is recorded in its history.

## 5. Conditioning, solver, tolerances, and stopping

### 5.1 Exact change of finite row basis used for conditioning

For each rational row, subtract its known expansion coefficients along all retained H2..H40. For cosine, the coefficient at H_(k+2) is r^k*T_k(c); for sine it is r^k*U_(k-1)(c) for k>=1 and0 for k=0. At k=0 subtract m(-1) as well. This produces an invertible triangular row transform T. Optimize with row matrices `T@Urows,T@Vrows` and initial coefficients `solve(T.T,lambda_physical)`; save final physical coefficients `T.T@lambda_transformed`. It introduces no nulls and removes none. Source: `numerical_direction_search.py:35`, `:141`, `:145`, `:152`.

### 5.2 Newton/SOCP subproblem

The algorithm minimizes the finite-grid convex disk support in lambda. At each accepted iterate:

- form the exact analytic numerical gradient and Hessian of that finite sum;
- diagonal-equilibrate the Hessian using `sqrt(max(diag(H),1e-40))`;
- eigendecompose the symmetrized equilibrated Hessian and floor eigenvalues at `max(max_eigenvalue*1e-10,1e-12)`;
- solve a whitened quadratic step problem with the sign inequalities and a Euclidean trust region, using CVXPY/CLARABEL;
- scale each QP constraint by `max(norm(row),1e-30)`;
- initialize the trust radius at10, increase it by1.4 up to100 after alpha>.9, halve it down to .05 after alpha<.1;
- impose a .995 safety factor on the step length that would hit a sign constraint;
- backtrack by .5 for at most70 trials until finite objective, feasibility violation<=max(eps*.01,1e-13)=2e-8, and Armijo decrease with coefficient1e-4 all hold;
- stop if the descent is nonnegative, feasible alpha<=1e-15, QP failure, or line-search failure; these statuses are not silently treated as convergence.

The QP solver uses max_iter=250, absolute and relative gap tolerances1e-10, and feasibility tolerance1e-11. The **bridge-stage copy** also sets CLARABEL `max_threads=2` (`research/directional_certification/numerical_dual_newton.py:68`). The preceding16-step refinement used `research/improvement_dual_newton.py`, whose CLARABEL call leaves max_threads at the solver default. The current generic driver was extended during the campaign; recreating the first stage with the current generic driver can therefore use a different thread setting. The original fixed-g0 driver remains available for that stage.

If the initial point violates sign constraints, the code first solves a row-scaled minimum-step projection with extra target margin max(eps*.01,1e-12). For coefficient i, its projection scale is the maximum of the Euclidean norm of `u_i*sqrt(weight/p)`, the norm of `v_i/p`, and1e-20. Projection feasibility is accepted only if the residual is<=max(eps*.1,1e-10)=2e-7. Projection failure or excessive residual returns a failure status. The final bridge start was inherited from a feasible nearby witness. Sources for these details: `numerical_dual_newton.py:54`, `:77`, `:91`, `:95`, `:120` through `:173`.

The generic direction driver also has a warm-start-only adjustment `lambda_X4 += (32/9)*max(b_new-b_old,0)` because X4/p tends to -9/32 at the spin-two threshold. This is an initialization choice, not a constraint on the optimized coefficient. It is identically zero for both endpoint stages documented here, whose old and new b are0 (`numerical_direction_search.py:66` through `:72`). No nonzero-angle optimization was completed in this newest numerical campaign.

The convergence test requires **both** `old_value-value < ftol*max(1,abs(value))` and `-descent < gtol^2`. A small accepted decrease alone is not the convergence criterion. The newest stage did not reach or claim this test: it was stopped after six accepted steps and the saved checkpoint was frozen. All six accepted alpha values are1; the last recorded threshold violation is -1.191935439237568e-11. Its original planned max_iter was8; max_iter=6 is the convenient replay equivalent for the saved accepted count.

## 6. Complete spectral grids and internal source quadratures

### 6.1 Spectral variable and positive weights

Use beta=sqrt(1-4/mu), mapping mu in (4,infinity) to beta in (0,1). Then

```
dnu_l = 4*(2*l+1)/pi * d beta/(1-beta^2)^2.
```

Each beta panel uses Gauss-Legendre nodes x,w on (-1,1). The first panel [0,.12] maps linearly. Every subsequent [a,b] panel uses `angle=pi*(x+1)/4`, `beta=a+(b-a)*sin(angle)^2`, with derivative `(b-a)*pi/4*sin(2*angle)`. The grids contain no endpoint nodes and no node at infinity. Mapping an entire domain is different from rigorously bounding the missed behavior between finite nodes.

### 6.2 Rational-peak panel boundaries

The initial set is `{0,.12,.35,.95,.995,1}` plus beta values for source energies11.9,12,24,48. For every rational row with r>=.94 and abs(c)<1, add

```
center = tan(acos(c)/4),
width = (1+center^2)*(-log(r))/4,
center-2*width, center, center+2*width,
```

keeping additional values strictly between .12 and .995 and deduplicating. The final inventory gives **56 panels**. Source: `research/improvement_peak_audit.py:16`, `:34`. The full ordered numerical edge list is in the final head JSON and repeated in the appendix below.

### 6.3 Training stages

| Stage | Spectral spins | First head panel | Later panels | Separate bridge first panel | Source energy/angular orders | Physical cells incl. zero-cost threshold probes |
|---|---|---:|---:|---:|---|---:|
| Original enriched266-row training | 0..24 even |80|16|not integrated|64/56|12493|
| First refinement `b0_refined` |0..24 even|80|24|not integrated|80/80|18213|
| Final bridge stage |0..24 and26..40 even|80|16|20|64/64|19693|

The bridge stage has13*(80+55*16+1)=12493 head cells including13 zero-weight threshold probes, plus8*(20+55*16)=7200 bridge cells. Each training stage also appends414 artificial guard cells. The final planned run had20107 total cells. Physical integration weights for artificial/threshold probes are exactly zero.

### 6.4 Finest head replay

The final independently assembled head replay uses n=32 on the55 later panels, first head order160, source energy/angular orders96/96, and all even spins0..40. Thus there are1920 spectral nodes for each of13 head spins and1780 nodes for each of8 bridge spins, for39200 integrated cells. The head's smallest sampled energy is4.0000000001814335; the bridge's is4.000000679912664; the largest is86617343.18372007. Those extrema are quadrature nodes, not physical cutoff declarations.

**Metadata caveat:** the nested `head.first_panel_order` in the head JSON says32 because `angle_audit.panel_support` records its input n after the peak-grid function has been substituted. The actual head first-panel order is160: `numerical_campaign_audit.py:33` sets `peaks.HEAD_FIRST_ORDER=160`, `improvement_peak_audit.py:38` uses it, and top-level `head_first_panel_order=160` plus1920 nodes verify this. Do not replay with a32-node first panel on the basis of that stale nested field. The bridge nested first-panel value20 is correct.

### 6.5 Source-window added-spin replay

Spins42,44,...,100,8 panels with physical energy edges `[4.02,4.1,6,8,11.9,12,24,48,49]`; each panel maps to beta and uses the sine-squared map above with48 Gauss nodes, for384 energy nodes per spin and11520 integrated cells. Source/angular orders are80/80. Raw threshold Legendre values prohibit extending this unscaled routine arbitrarily close to4 at high spin. Source: `research/improvement_source_tail_audit.py:12`.

### 6.6 Entire-UV added-spin replay

Use x=1/mu and7 panels with exact intended endpoints `[0,1e-6,1e-4,1e-3,1/200,1/100,1/60,1/49]`. Each panel has128 linearly mapped Gauss nodes; spins42..100 even give896 nodes per spin and26880 integrated cells. Weight is `dx/[x^2*2*pi*sqrt(1-4/mu)]` times `(2*l+1)`. Source/angular orders are80/80. No finite UV cutoff is applied.

The dedicated evaluator `research/existing_solver/scaled_spin_audit.py:19` combines source weights first and uses the nonsingular, normalized Legendre ratio expression. It requires mu above every source window plus1 and spectral spin above every Roy J; both hold here. Sources: `research/improvement_uv_tail_audit.py:9`; `scaled_spin_audit.py:62`. Do not add the earlier restricted tau=l^2/mu tail diagnostic to this result: it is contained in the full UV integral.

### 6.7 Internal source quadratures and pole treatment

`build_roy_bank` uses tensor-product Gauss-Legendre source energy and folded-angular quadrature with analytic PV pole extraction. `build_fixed_transfer_bank` uses tensor-product source energy and transfer quadrature, with analytic Legendre-Q energy integrals in/near source windows. Taper combinations are made algebraically before calling/summing untapered rows. Internal source orders64,80,96 above are separate from spectral quadrature orders16,24,32,48,128. Raising one does not raise the other. Neither the tensor quadratures nor the special-function evaluations are interval-enclosed in this numerical campaign.

## 7. Every final numerical check and its scope

1. **Low-spin threshold probes:** delta={1e-4,1e-6,1e-8,1e-10,1e-12}, spins0..20 even; pass/fail uses the final delta for nonzero spins with zero sign tolerance. Final maximum is -1.999071669089858e-6. These probes have no quadrature weight.
2. **Extended threshold probes:** delta={1e-2,1e-4,1e-6}, spins2..40 even; pass/fail uses delta1e-6 with zero sign tolerance. Final maximum is -1.197925189000415e-5.
3. **Integrated high-p guard replay:** on every integrated head/bridge cell with p>1e4 and nonzero spin, no V/p>0 was found. All integrated support arrays are finite. These checks do not sample all energies.
4. **Source-window added spins:** all30 added spins are finite; the sampled positive-V contribution is0 and max sampled V/p=-410.802652 approximately.
5. **UV added spins:** sampled positive-V contribution0, max sampled V/p=-15.5995592 approximately.
6. **Training versus finer replay:** +4.049988322485376e-7 at the same spin cutoff40. The training and fine source/angular/spectral quadratures all differ; this is not a one-parameter error estimate.

Sources: `research/angle_audit.py:168` through `:229`, all three final JSON audit records. The head status is `finite_head_audited_tail_unproved`, with `all_spin_certified=false` and `continuum_certified=false`.

The earlier exact positive-Bernstein Roy certificate `output/improvement_20261006/right_final_roy_sign_certificate.json` and Arb real-tail record `output/certification_20261006/right_enriched_A100_arb.json` apply to the **earlier** SHA `032ed315...` witness, not the present SHA `fab65bff...` witness. They cannot be reused as final-witness certificates without recomputation. The extended transfer-row validation `output/improvement_20261006/extended_transfer_validation.json` tests the generic row implementation through j=4 and untapered n=9; it is a numerical algebra validation, not an interval error proof.

## 8. Historical chain, inputs, and optimizer stages

| Stage and saved witness | Row count | Recorded training objective | Status / changes |
|---|---:|---:|---|
| `output/angles_resummed/angle_000p000.npz` |150|see original records|starting witness; independently certified in a separate campaign|
| `output/improvement_20261006/right_peak_conditioned.npz` |187|2.70562043787099|adds37 rational rows, peak-aware conditioning;3 accepted steps; recorded numerical convergence|
| `output/improvement_20261006/right_roy18j8_stopped.npz` |222|2.705563363436272|adds Roy J=8 and n through18;17 accepted steps; stopped checkpoint|
| `output/improvement_20261006/right_roy18j8_strip7j4.npz` |266|2.698815096042|adds strip j through4,n through7;20 accepted steps; iteration limit|
| `output/directional_certification/numerical/b0_refined.npz` |266|2.698812924859148|same rows; n24/source80 refinement;16 accepted steps; iteration limit|
| `output/directional_certification/numerical/bridge/bp0_stopped.npz` |266|2.6988810737074695 transformed objective|same rows; trains spins26..40 too;6 accepted steps; stopped checkpoint|

Training objectives across different spin cutoffs must not be compared as if they were the same objective. Adding the bridge raises the starting objective to2.6989400302853745 before reoptimization. After the first bridge step it is2.698881593155981; after six it is2.6988810737074695.

Recorded row additions are not complete inventories. For example, `right_roy18j8_strip7j4.json` records roy_spin=6 and roy_order=-1 because that stage added no Roy rows; the warm witness already contained J=8,n=18. The saved names are authoritative.

The rejected `right_enriched.npz` ordinary-grid result2.6702534 is **not** part of this chain: rational-peak replay raised it to about2.769334. It was quadrature aliasing. See `output/improvement_20261006/rejected_naive_rational_enrichment.json` and `research/improvement_numerical_report.md`.

### Stage-specific settings

- 187-row stage: peak n16, lmax24, hmax40, rational enrichment enabled, source64/angular56, eps2e-6, eigen_floor1e-10, max_iter40. Its record does not retain explicit ftol/gtol. The current original driver defaults are3e-11/1e-6; these are a reconstructive setting, not newly discovered evidence of an unlogged historical CLI value.
- 222-row stage: same n16/source64/angular56, Roy order18/spin8; saved record explicitly retains a stopped accepted iterate16 (17 steps). Not every original CLI option is retained. Use the saved input/checkpoint for exact fixed-witness reproduction.
- First266-row stage: n16/source64/angular56, strip order7/j4, eps2e-6, eigen_floor1e-10, ftol1e-7, gtol5e-4, max_iter20.
- `b0_refined`: n24/source80/angular80, eps2e-6, eigen_floor1e-10, ftol1e-9, gtol1e-4, max_iter16; transformed objective2.698812924859046 versus physical-basis2.698812924859148, discrepancy1.021405182655144e-13.
- Final bridge: n16/source64/angular64, eps2e-6, eigen_floor1e-10, ftol1e-8, gtol1e-4, planned max_iter8, interrupted after6 saved accepted iterates. The compact stopped record omits some options, but the launched campaign and source `numerical_direction_campaign.py:28` through `:37` fix these settings. The final frozen NPZ is copied from `bp0p000000_checkpoint.npz`; there is no claim that a full normal-termination optimizer result was written.

## 9. Reproduction commands

### 9.1 Environment

The executed numerical stack was Python3.12.14, NumPy2.3.5, SciPy1.18.1, CVXPY1.9.3, CLARABEL0.11.1, Windows64-bit. Python executable:

```powershell
$py='C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
Set-Location 'C:/Users/laure/Documents/ChatGPT/Unitarity codex'
```

The search modules set OPENBLAS_NUM_THREADS and OMP_NUM_THREADS to2 only if those variables are not already set. Optimizer final-bit trajectories may vary with solver/library versions and threading. The authoritative object is the fixed binary64 witness, whose complete numerical replay is independent of the optimization trajectory.

### 9.2 Preferred reproduction: fixed witness to all three measured components

Use a fresh output directory; do not overwrite archived records.

```powershell
& $py research/directional_certification/numerical_campaign_audit.py output/directional_certification/numerical/bridge/bp0_stopped.npz --outdir output/replication_numerical/fine --n 32 --threshold-order 160 --kernel-order 96 --angular-order 96 --tails
```

Equivalent separate calls, useful for diagnosing a mismatch:

```powershell
& $py research/improvement_peak_audit.py output/directional_certification/numerical/bridge/bp0_stopped.npz --n 32 --threshold-order 160 --kernel-order 96 --angular-order 96 --out output/replication_numerical/head.json
& $py research/improvement_source_tail_audit.py output/directional_certification/numerical/bridge/bp0_stopped.npz --n 48 --kernel-order 80 --angular-order 80 --out output/replication_numerical/source.json
& $py research/improvement_uv_tail_audit.py output/directional_certification/numerical/bridge/bp0_stopped.npz --n 128 --kernel-order 80 --angular-order 80 --out output/replication_numerical/uv.json
```

Sum `head.support + source.full_support + uv.full_support` exactly once. Do not add the previous restricted-tau UV audit, the same-spin training objective, or any prior witness's tail.

### 9.3 Reproduce the two newest optimization stages

The first line intentionally uses the original fixed-g0 driver/optimizer, preserving its default CLARABEL threading behavior. The next line uses the bridge-capable copy and an explicit6-step limit in place of a manual interruption.

```powershell
& $py research/improvement_right_search.py --warm output/improvement_20261006/right_roy18j8_strip7j4.npz --peak-panels --n 24 --kernel-order 80 --angular-order 80 --strip-rational-harmonics --eigen-floor 1e-10 --max-iter 16 --ftol 1e-9 --gtol 1e-4 --out output/replication_numerical/b0_refined
& $py research/directional_certification/numerical_direction_search.py --warm output/directional_certification/numerical/b0_refined.npz --b 0 --peak-panels --train-bridge --n 16 --kernel-order 64 --angular-order 64 --strip-rational-harmonics --eigen-floor 1e-10 --ftol 1e-8 --gtol 1e-4 --max-iter 6 --out output/replication_numerical/bridge_from_frozen_intermediate
```

To test an entirely regenerated trajectory, change the second command's warm path to `output/replication_numerical/b0_refined.npz`. A numerically nearby result is expected; exact matching of the archived NPZ's ZIP bytes is not required to reproduce the optimization problem. For exact archived coefficients, replay the frozen witness or reconstruct the main handoff's coefficient appendix instead.

### 9.4 Regenerate the larger row bank from the original150-row seed

These commands reconstruct the reported mathematical stages using retained defaults and explicit conservative stopping counts. They do not assert that missing historical CLI settings or solver final bits have been recovered. For strict input continuity, each command below starts from the **archived** preceding checkpoint. To run a new connected trajectory, substitute the preceding newly generated output.

```powershell
& $py research/improvement_right_search.py --warm output/angles_resummed/angle_000p000.npz --enrich --peak-panels --n 16 --kernel-order 64 --angular-order 56 --strip-rational-harmonics --eigen-floor 1e-10 --eps 2e-6 --ftol 3e-11 --gtol 1e-6 --max-iter 40 --out output/replication_numerical/rational187
& $py research/improvement_right_search.py --warm output/improvement_20261006/right_peak_conditioned.npz --roy-order 18 --roy-spin 8 --peak-panels --n 16 --kernel-order 64 --angular-order 56 --strip-rational-harmonics --eigen-floor 1e-10 --eps 2e-6 --ftol 3e-11 --gtol 1e-6 --max-iter 17 --out output/replication_numerical/roy222
& $py research/improvement_right_search.py --warm output/improvement_20261006/right_roy18j8_stopped.npz --strip-order 7 --strip-j 4 --peak-panels --n 16 --kernel-order 64 --angular-order 56 --strip-rational-harmonics --eigen-floor 1e-10 --eps 2e-6 --ftol 1e-7 --gtol 5e-4 --max-iter 20 --out output/replication_numerical/strip266
```

This document does not claim an independently logged, bit-reproducible zero-initialization path to the original150-row seed. The exact final266 saved coefficients suffice to reconstruct and audit the reported numerical candidate without that search history. A truly new search from zero is a new optimization experiment, even if it uses exactly the same266 rows and restrictions.

### 9.5 Kernel validation before changing constraints

Read and, in an isolated copy or with new output paths, rerun the existing independent algebra checks: `roy_kernels.py:258` (divided differences, real Legendre-Q PV, scalar smearing), `transfer_kernels.py:56` (independent scalar PV and order refinement), `crossing_kernels.py:116` (exact coefficient elimination and crossing-symmetric pole test), `improvement_transfer_audit.py` (extension through strip j4/untapered n9), and `numerical_dual_newton.py:207` (finite-difference gradient/Hessian checks and independent small conic optimum). Some self-tests write adjacent historical audit JSON paths; inspect their destinations first instead of overwriting a frozen record. These checks validate identities/implementation numerically, not the all-spin outer bound.

## 10. What is required to certify this266-row candidate

The minimum mathematical target is a **complete, outward-rounded upper bound** on its disk-support sum over every even spin and every mu>4, under the stated amplitude assumptions. An optimizer convergence certificate or a finite-grid duality gap alone cannot supply this.

1. **Freeze the candidate and interpretations.** Bind all component records to the exact witness SHA and code hashes. Interpret each saved multiplier as its exact binary64 rational. For a new exact proof use the row-label parameters consistently as exact decimals and T=4/3; do not silently substitute the rounded literal1.33333. Compare the resulting exact kernel convention with the numerical reconstruction and state it explicitly.
2. **Extend the interval kernel parsers and source moments.** The old150-row `head_certification_spin0.py` parser is restricted and does not simply accept all new rational-sine rows. Extend U/V and derivative formulas for sine rows, Roy J through8/n through18, and strip j through4/n through7. Regenerate exact/Arb source and UV moments for the **new coefficients** and independently verify their enclosures. Do not reuse an old-witness moment cache merely because the source windows match.
3. **Certify finite low-spin integrals, including endpoints and all source transitions.** A natural split is direct interval integration for l=0,2,4,6,8, because the new Roy family has a local real-part contribution at l=8. If retaining a separate8..40 block instead, explicitly include that nonzero Roy real part there. The old finite-V proof assumed no Roy local real part for l>=8 because its largest J was6. Cover threshold gaps, source-window boundaries/PV logarithms, rational peaks, and the entire UV endpoint. Stable analytic PV cancellation and geometric source remainders must be validated, not only sampled. Use interval upper Riemann or the established convex-interpolation/Taylor bounds, with directed rounding and independent partition replay.
4. **Certify every remaining finite spin.** For example, handle10..40 via interval integration or a proved local coercivity inequality V_l<=-kappa_l(mu)*p_l together with a certified real norm. When kappa varies by spectral cell and spin, merge partitions using exact overlap tests and retain a safe verified minimum. The existing `certified_finite_v.py:150` hardcodes the old witness hash and `:174` asserts cosine-only rational rows. It is not a generic verifier for this file.
5. **Prove the all-spin tail over the whole energy ray.** Establish an exact dangerous-source Roy polynomial sign, a compact-energy absorptive envelope, a joint large-spin/large-energy UV envelope, and a strictly positive coercivity floor on the chosen remaining domain. The new sine rational rows and enlarged J,j degrees change the analytic error budgets; the old fixed-tail proof's cosine-only bounds do not automatically cover them. No finite set of tested spins or tau values proves this statement. A possible efficient split starts beyond l=40; choosing l=100 only moves more finite work into the head and still requires the infinite remainder.
6. **Bound the real-kernel tail.** Once V<=-kappa*p is proved, use `H(U,V)<=U^2/(2*kappa*p)` and rigorously integrate/sum the real norm. Strip real kernels persist at arbitrarily high spin. `research/certify_real_tail.py` already parses cosine/sine rational rows and the larger source families, but its output by itself is only a real-norm component. It does not prove the needed V coercivity or the head. For other tail techniques, every use of a global g2 cap must be explicitly supplied and certified; none is part of the present numerical result.
7. **Close every omitted domain and assemble exactly.** In particular include4<mu<4.02 for spins42..100, all even spins>100, all cell boundaries of measure zero handled consistently, and every threshold/UV gap. Add only disjoint domain contributions, use outward rounding for each bound and final decimal, and perform an independent replay with hashes/precision/partition/NaN checks. A final record must state every assumption, component, exact coefficient convention, and remaining limitation. If any component is absent, call it a partial certificate rather than a bound.

The separate150-row certification framework is a useful implementation template, not a certificate for the new witness. Its current absorptive kernel/UV evaluators are deliberately specialized; copying the new filename into old CLI options cannot make their mathematical assumptions true. The relevant conceptual requirement, also underlying a rigorous S-matrix dual bound, is admissibility and integral control over the **complete continuum**. Sampling, source-sign checks, real-tail bounds, and solver feasibility are useful components but individually insufficient.

## 11. Numerical record and code hashes

The following appendices were extracted directly from the current saved files during this audit. They supplement the source path/line citations above. The source paths identify exact implementations; line numbers are meaningful for the listed SHA versions.


### 11.1 Input, source, and record SHA-256 values

| Path | SHA-256 |
|---|---|
| `output/angles_resummed/angle_000p000.npz` | `bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c` |
| `output/improvement_20261006/right_peak_conditioned.npz` | `865439b3c42ac3689b7c52379cdb4bad8f84dd15d87e61c6372bef3ecdff9665` |
| `output/improvement_20261006/right_roy18j8_stopped.npz` | `0cc140a713be3b7884959fa186981cf75dfcf9ef2c9e303463f79eb758c1b604` |
| `output/improvement_20261006/right_roy18j8_strip7j4.npz` | `032ed315425850abd5ac75557795d959141672ba26e2712f2b46bb42a61c03c8` |
| `output/directional_certification/numerical/b0_refined.npz` | `94e17d38dd9d716fb16db65f856e1a7ed13d8180e90d9321cea8b0144b8a0baa` |
| `output/directional_certification/numerical/bridge/bp0_stopped.npz` | `fab65bffd85740b1d6043124bd1f2a0324547e1c033f98274340ab0b92beb2e4` |
| `research/directional_certification/numerical_direction_search.py` | `dac16dce1bcf0bc632bbf1c4f883c63f813cf6d6c5109d93938836cff90ed6f4` |
| `research/directional_certification/numerical_direction_campaign.py` | `dee736a0e5d2f8d14cac59b241b001749fff82e20e782e1e59bdda5876695cd5` |
| `research/directional_certification/numerical_dual_newton.py` | `9b7c4b6a795505ddb3af6800f39bee7e5b5f9ebca363e13b427318651f79b740` |
| `research/directional_certification/numerical_campaign_audit.py` | `b33f70921489e1dd75835de5b054ffe2cf0a28cc8206bf3c70e47c12d269a0cd` |
| `research/angle_audit.py` | `5da3fa09e791fc84f0a50a91c5db0d99d22014999a89c01433d44a306dbcc7ac` |
| `research/roy_kernels.py` | `900dad1c5885a7f86a83a179a70d865ba8cf71df3fe1c8f71537094307fc650e` |
| `research/transfer_kernels.py` | `20a622d8291e42114dbaf944358d777a010a6288df50a1a302680343a1ad358c` |
| `research/existing_solver/massive_forward_kernels.py` | `5df66ff7995d501333c05a0b5a9037fdf91c34eb9d69f959e46af3e19152000a` |
| `research/crossing_kernels.py` | `5c2921c37d62faa12f7c267328789e420479375abc7c847da952c362f303066d` |
| `research/improvement_right_search.py` | `6c9ea1966a82da3e53e6a2ee306b308722fb46d74927142522c64f8238159cf4` |
| `research/improvement_dual_newton.py` | `07d692f4cd5cf25a2ad4b036c2f1eecb5c121267f64d03fe0b34dc214731bf36` |
| `research/improvement_peak_audit.py` | `fd6c7194b3f9069487ce0997680aed142f0e768559a27ac65563e39775037aef` |
| `research/improvement_source_tail_audit.py` | `4039b158e84ef412ff048747c3045a8af82330f2aa77850fff7203cc39e811ca` |
| `research/improvement_uv_tail_audit.py` | `d52c84197ed6fd0d470a1c5886cd4a95153361006184a43be8aa1a54c3d83843` |
| `research/existing_solver/scaled_spin_audit.py` | `a80043fc4644c9485268eeedf4975983a130effd6eabd587585bcc5d31200a35` |
| `output/directional_certification/numerical/fine/bp0_stopped_head.json` | `da7a097120b77dccf49e7f1ed29bce4fc6e7170cf91d232cebb41f96db9c732b` |
| `output/directional_certification/numerical/fine/bp0_stopped_source.json` | `a512418b06d665df67faa22ae95cd696d85ad289fb43a1d1a1047167e3f3ef4d` |
| `output/directional_certification/numerical/fine/bp0_stopped_uv.json` | `0fb0bfc3f41755ca584864a5e2ae63318a1e3dd04b4f834acdf40e886f011c3a` |
| `output/directional_certification/numerical/bridge/bp0p000000.history.json` | `ea38267e675e3f1e24399275ed38862bd1131e01fad1af6dcdb3e32122dbadd7` |

### 11.2 Ordered beta-panel edges (57 edges,56 panels)

```json
[
  0.0,
  0.12,
  0.12346789094483634,
  0.1392125158835187,
  0.15617332625549285,
  0.17223737327630112,
  0.18267581368159952,
  0.19311425408689792,
  0.2091783011077062,
  0.24046119881457073,
  0.2571225942098194,
  0.26794919243112275,
  0.2787757906524261,
  0.29543718604767477,
  0.31335235405830575,
  0.33071587362300253,
  0.3419987131196464,
  0.35,
  0.3532815526162903,
  0.3706450721809871,
  0.3841666461796637,
  0.4023790904231463,
  0.41421356237309503,
  0.4260480343230438,
  0.4442604785665264,
  0.45850208319613267,
  0.477784592858123,
  0.49031439482586836,
  0.5028441967936137,
  0.5221267064556041,
  0.5431547395979255,
  0.5638817976446129,
  0.5773502691896258,
  0.5908187407346388,
  0.6115457987813262,
  0.6531852075066038,
  0.6761548127590452,
  0.6910804946404703,
  0.7060061765218953,
  0.7289757817743368,
  0.754642865495115,
  0.7800518343952566,
  0.7965626182145923,
  0.8083208415456143,
  0.813073402033928,
  0.814779446364774,
  0.816496580927726,
  0.8351667221849756,
  0.8384823709340696,
  0.8526112132785795,
  0.8700557043721834,
  0.8969015850115446,
  0.9128709291752769,
  0.95,
  0.9574271077563381,
  0.995,
  1.0
]
```

### 11.3 Lossless coefficient extraction and reconstruction

This is an alternative to copying all266 decimal coefficients by hand. It preserves row order and every binary64 multiplier bit; decimal formatting should not be used as a truncation step. A re-created NPZ may have a different ZIP-container hash even if its arrays agree bit for bit. For certificate records tied to a byte-level SHA, retain the archived original NPZ or explicitly bind a newly verified container.

```python
import json
from pathlib import Path
import numpy as np

src = Path('output/directional_certification/numerical/bridge/bp0_stopped.npz')
d = np.load(src)
rows = [{'name':str(n), 'coefficient_hex':float(c).hex()}
        for n,c in zip(d['names'],d['multipliers'])]
archive = {'direction_hex':[float(x).hex() for x in d['direction']],
           'rows':rows}
# Write to a NEW path:
Path('output/replication_numerical').mkdir(parents=True,exist_ok=True)
Path('output/replication_numerical/exact_coefficients.json').write_text(
    json.dumps(archive,indent=2), encoding='utf8')
# Reconstruct from the text archive:
names = np.array([r['name'] for r in archive['rows']])
coefficients = np.array([float.fromhex(r['coefficient_hex'])
                         for r in archive['rows']],dtype=np.float64)
direction = np.array([float.fromhex(x) for x in archive['direction_hex']])
assert np.array_equal(names,d['names'])
assert coefficients.tobytes()==d['multipliers'].tobytes()
assert direction.tobytes()==d['direction'].tobytes()
np.savez_compressed('output/replication_numerical/reconstructed_witness.npz',
                    names=names,multipliers=coefficients,direction=direction)
```

### 11.4 Scope of the documentary audit

This audit read the row builders, optimizer, saved configurations, accepted-step histories, final NPZ arrays, and all three final numerical replay records. It did not rerun an expensive optimization, change any frozen code/record, or create a new continuum certificate. The source line numbers and hashes above identify the implementations actually inspected. Unlogged early-stage stopping tolerances and bit-identical optimizer history are identified as such rather than silently filled in.
