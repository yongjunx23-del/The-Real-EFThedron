# Lower-g0 numerical enrichment, 6 October 2026

This note records the new lower-g0 search. None of the head objectives in this note is, by itself, a complete exclusion bound. The root task separately assembles the fixed-witness audits and certified components.

## Objective and inventory

The normal is exactly `(-1,0)`: minimizing the complete dual support would prove `g0 >= -C`. Smaller positive `C` is stronger. The pointwise support of the exact partial-wave disk remains `V + sqrt(U^2+V^2)`; negative-V evaluation uses the cancellation-free identity `U^2/(sqrt(U^2+V^2)-V)`.

The input was the authoritative 131-row witness `output/angle_constraints/ming0_scaled_threshold_repaired.npz`. Its SHA256 is `c7efae0020c9883d6949978ff924d11c3b828f470f10e24ff813f993e47cc06a`.

The new inventory has **220 physical null rows**:

- Forward harmonics H2 through H40: 39 rows.
- Crossing X4, X5, X6, X7_1, X7_2, X8_1, X8_2: seven rows.
- Roy smearing at endpoint 11.9, even source spins J=0,2,4,6, Legendre orders n=0,...,14: 60 rows.
- Two transfer strips [12,24] and [24,48], width exactly 4/3, power-one taper, angular indices j=0,...,3, source orders n=0,...,7: 64 rows.
- Threshold-regular derivatives Dt1, Dt2, Dt3, H_j for j=2,...,12 and one G_regular per order: 36 rows.
- Resummed rational cosine and sine generators at r=.8 and c=-.9,-.75,-.5,0,.5,.75,.9: 14 rows.

No fractional harmonic generators were added. No new physical analyticity assumption was added. These rows are the established manuscript-derived sum rules, with larger basis sizes. An invertible harmonic-tail change of coordinates subtracts the H2...H40 components of rational generators during optimization; saved multipliers are always transformed back to the original physical rows.

## Numerical objective and search restrictions

`lower_search.py` assembles kernels with `angle_audit.normalized_inventory` at kernel and angular quadrature orders64. The objective contains every even spin0,...,40. It uses the existing split velocity grid, with 48 Gauss nodes per panel. The first threshold panel for spins26,...,40 has order20 to avoid overflowing raw Legendre polynomials. This gives numerical integrals over the full spectral energy interval for those spins; it is not a proof of quadrature error or of their threshold endpoint neighborhoods.

The new head training support started at7.781783183182079 for the old witness on this grid. The objective does not contain physical spin42+ costs. Zero-weight artificial cells restrict the searched multipliers; they are not additional physical sum rules.

Restrictions comprise the original threshold-dominated head points, Roy source positivity samples, the necessary highest-degree crossing sign, dense moderate-energy grids through spin200, and joint-scale grids `mu=4+ell(ell+1)/tau` for every even42,...,100 and ell=120,160,200,300,500,1000. The joint grid covers tau=.001,...,150, with dense steps near the previously troublesome bands. These tail grids require `V/p <= -1e-4` in the later runs. This sampled margin does not bound energies or spins between samples.

Initialization uses a scaled L1 linear-program projection with HiGHS. The Newton objective Hessian is diagonally equilibrated and whitened, with eigenvalue floor1e-7. Convex QP subproblems use CLARABEL; a full-grid line search retains feasible decreases. Solver failure or stagnation is reported without an optimality claim.

## Important rejected candidates

The richer basis rapidly lowered the head objective. `lower_enriched_step4.npz` has head support7.735893303638799 in an independent n96 replay (kernel and angular orders80). It is **not an acceptable all-spin witness yet**: exact Roy Bernstein checking located an inter-node source-sign failure. Minimizing the source weight at z=0 near s=7.0605887 gives W approximately -3.19e-6, even though a251x251 grid had returned positive values.

The subsequent filename `lower_source_safe_frozen.npz` is also misleading and must **not** be interpreted as certified source-safe. Merely imposing W>=1e-4 on a41x21 source grid still allowed W approximately -0.00748228859 at s=7.12479428,z=0. Its name records the attempted repair, not a successful certification. `lower_source_safe.npz` and its checkpoint likewise require rejection or repair before any all-spin conclusion.

This is why a dense grid alone was not accepted as a source-sign certificate.

## Continuum source and crossing cone restrictions

The matrix agent generated9408 linear Bernstein inequalities on16 source-energy cells times4 angular cells. The cleared Roy polynomial is

`(s-4)^6 sum_{J,n} lambda_{J,n} P_n(2(s-4)/(79/10)-1) P_J(q*(1-(8/3)/(s-4)))`,

with `20/3 <= s <=119/10`, `0<=q<=1`. Its degrees are20 and6. Each inequality is divided by the corresponding positive Bernstein coefficient of `((s-4)/(79/10))^6`, so a shared search margin directly implies the same lower bound on W. The floating matrix is only used to search; saved multipliers are checked again with exact rational arithmetic.

The driver also accepts42 exact crossing-series inequalities from `matrix_crossing_series.py`. They prove that all positive-Legendre-series coefficients beyond k>=10 have a negative sign, uniformly over all physical spectral energies. This is a rigorous component sign statement, not a sign proof for the complete V kernel.

`lower_bernstein_projected.npz` is the first fixed enriched candidate that passed the new exact continuum source check:

- SHA256: `b87b1be5694f05a3c0be6807935933a82d2d020331b968e30e70ee7fffd78cf5`.
- Training head support:7.739540786230227.
- Exact source certificate: `lower_bernstein_projected_roy.json`, W>=0.00010001010002979653.
- Exact crossing certificate: `matrix_crossing_bernstein_projected_certificate.json`, F(k,1/Q)<=-117.12391616976231 for all k>=10 and all energies.

The main task also formed a90% convex combination of the enriched frozen candidate with the original authoritative witness. This combination is `blend_0.9.npz`, and its source sign was independently certified. The root task's fixed-witness audits found a smaller measured support than the projected Bernstein point, so the blend is the comparison candidate unless a subsequent complete audit establishes a better one.

The Bernstein Newton refinement initially appeared stalled, then returned four accepted steps before it was stopped at the parent task's request to finalize the audited result. The retained future candidate is `lower_bernstein_partial_frozen.npz`, SHA256 `eb06cd5857961a9bf104710e9af475eb92e61222b5d5019b45a69a343d21e6d2`. Its n48 training support is7.733304374326462, and exact source replay proves W>=0.00010001010293763446 in `lower_bernstein_partial_roy.json`. It is not used to supersede the audited blend without a complete fixed-witness replay, and there is no claim of optimizer convergence.

A subsequent fixed-witness numerical replay found head n96=7.733664192832357, all mu>49 at even42,...,100 contributing0.0005733352287886019, and 4.02<mu<49 at those spins contributing0.0000130640820696505. The measured total is **7.734250592143215**. Its UV integrand includes small positive-V islands: maximum sampled V/p=0.0006544747024578475, included positive-V support0.00004975264439261966. Thus this slight numerical improvement over the root's7.7348690521 blend does not improve the simple negative-V tail-certification route. The main reported comparison remains the more thoroughly audited blend. See `lower_bernstein_partial_audit_summary.json` for the exact omitted domains and files.

## Reproduction

Use the workspace Python with NumPy, SciPy, SymPy, CVXPY/CLARABEL and the existing research modules. The first unconstrained-source campaign can be reproduced with:

```text
python research/lower_g0_upper_g2_20261006/lower_search.py --out output/lower_g0_upper_g2_20261006/reproduction_enriched --source-margin 0
```

The source cone projection, using the saved matrix and authoritative row order, is:

```text
python research/lower_g0_upper_g2_20261006/lower_search.py --warm output/lower_g0_upper_g2_20261006/lower_source_safe.npz --out output/lower_g0_upper_g2_20261006/reproduction_projected --crossing-matrix --source-bernstein output/lower_g0_upper_g2_20261006/matrix_roy_guard16x4.npz --project-only
python research/lower_g0_upper_g2_20261006/matrix_roy_guards.py output/lower_g0_upper_g2_20261006/reproduction_projected.npz --verify --out output/lower_g0_upper_g2_20261006/reproduction_roy.json --margin 0
```

The cache `lower_enriched_cache.npz` stores the n48 objective and base grid guard matrices for exactly the220-row inventory. To change objective quadrature, use a new cache path as well as `--n`; never relabel this cache as a different quadrature order.

A complete theorem still needs a bound on the entire omitted spin/energy domain and directed control of finite-head integrals. Source, crossing, threshold-limit and joint-UV-limit certificates are distinct ingredients and must not be conflated with that completed theorem.
