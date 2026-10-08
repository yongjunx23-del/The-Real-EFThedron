# Left-sector 220-row search and audited candidates, 7 October 2026

The ten directly optimized normals are 120,135,150,165,180,195,210,225,240,255 degrees. `output/full_direction_20261007/left/direct_supports.json` and `.csv` are the joined results. They are fixed-witness numerical support integrals, **not support values of an intersected polygon and not complete continuum exclusion bounds**. All old research and certificates were read-only inputs.

For angle theta, the actual saved normal is `(cos(theta),sin(theta))` as binary64; at180 degrees it is exactly `(-1,0)`. The normalized amplitude convention is mass squared one and `M/(32 pi)=g0+g2*(sbar^2+tbar^2+ubar^2)+...`. A support C in the180-degree direction would imply `g0>=-C` only after the entire dual integral had been bounded.

## Results and their limits

| Angle | Measured direct support | Published dual, graphical | Published primal, graphical |
|---:|---:|---:|---:|
|120|4.768806046|4.944045398|4.333092628|
|135|6.202830540|6.432014866|5.646673353|
|150|7.219825532|7.481653151|6.575442621|
|165|7.737657204|8.021429137|7.056106341|
|180|7.733518866|8.014557984|7.055908074|
|195|7.196686130|7.461507949|6.574861333|
|210|6.175728434|6.399968478|5.645748657|
|225|4.728230348|4.902281731|4.331887540|
|240|2.962156528|3.070512586|2.722815446|
|255|0.993667958|1.039322823|0.928910673|

Benchmark values come from the existing vector extraction and `research/angle_benchmarks.py`, based on [The Phases of the Scalar S-Matrix Island](https://arxiv.org/html/2605.06613v1). They are graphical supports with extraction uncertainty, not exact numerical optimizer data. The primal and axiomatic dual analyticity classes differ as described in the existing benchmark notes. No universal comparison against every published dual method is asserted.

The180-degree measured value improves the previous main lower-g0 candidate7.734869052068362 by0.001350185636605. The new fixed witness is `angle180_final.npz`, SHA256 `e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.

The measured domains are:

- every even spin0..40, all spectral energies in endpoint-split numerical quadrature;
- every even spin42..100 on4.02<mu<49;
- every even spin42..100 on49<mu<infinity, using inverse-energy quadrature.

Remaining obligations are all even spins102 and above at every energy,4<mu<4.02 for spins42..100, and rigorous quadrature-error enclosures. The contributions omitted from a measured support are nonnegative. Consequently its displayed value alone is not an upper bound on the complete dual cost. No small unmeasured remainder is assumed.

## Exact row inventory

The authoritative order is the `names` array in each NPZ, matching `output/lower_g0_upper_g2_20261006/lower_enriched_cache.npz`.

| Family | Rows | Count |
|---|---|---:|
|Forward harmonics|H2..H40|39|
|Crossing|X4,X5,X6,X7_1,X7_2,X8_1,X8_2|7|
|Roy|B=11.9,J=0,2,4,6,n=0..14|60|
|Tapered transfer strips|[12,24],[24,48],T=4/3,power1,j=0..3,n=0..7|64|
|Regular transfer derivatives|Dt1,Dt2,Dt3, each Hj2..12 and G_regular|36|
|Rational resummations|cosine and sine,r=.8,c=-.9,-.75,-.5,0,.5,.75,.9|14|

There are no fractional-power generators. `T1.33333` is the legacy label for exactly4/3. Fixed-witness exact checks treat the multipliers as exact binary64 rationals and decimal row parameters as exact decimals. The invertible harmonic-tail coordinate transform is used only during optimization; saved multipliers refer to the physical rows above.

## Training objective and restrictions

`search_left.py` reuses the old n48/k64/a64 physical kernel cache for all even spins0..40. Its threshold bridge first panel has20 nodes at spins26..40. Cache load time is explicitly separated from new guard construction: the original cache-building cost was incurred in the previous campaign and is not relabelled as a new kernel build.

The objective is `sum W*(V+sqrt(U^2+V^2))`. For negative V it is evaluated as `U^2/(sqrt(U^2+V^2)-V)` to avoid cancellation. `U=a*u0+lambda@u`, `V/p=a*v0/p+b+lambda@(v/p)`. The generic driver reconstructs the original guard normalization from nonzero crossing rows and explicitly checks agreement with the old(-1,0) bases. Maximum observed reconstruction error was approximately1.1e-16.

The artificial guard cells have zero objective weight and p=100000. They restrict the search for a dual witness; they impose no extra physical constraint on an amplitude.

- Physical columns with p>10000 require V/p<=-1e-8.
-861 source samples cover s in[20/3,11.9], mapped angular coordinate in[0,1], with search source margin1e-4.
- One necessary highest-degree crossing coefficient guard is included.
- Moderate-energy grids include spins26..72 and80,100,140,200.
- Joint-scale grids use every even spin42..100 and120,160,200,300,500,1000 at tau=ell(ell+1)/(mu-4), with logarithmic and dense linear nodes from.001 to150. The new tail search margin is.002.
-9408 normalized Bernstein inequalities constrain the dangerous Roy source polynomial on16x4 boxes. The final saved multipliers are checked exactly again; floating search feasibility is not accepted as the source certificate.
-42 exact crossing-series inequalities come from `matrix_crossing_series.optimizer_guard_matrix`, cutoff k10, degree8.

The added artificial guard count is18744; with the physical high-p columns the optimizer checks21697 inequalities. No blanket X4<=0 restriction is imposed: higher crossing rows can compensate a positive X4 multiplier.

## Solver and runtime interpretation

The Newton Hessian is diagonally equilibrated and whitened with eigenvalue floor1e-7. CLARABEL solves trust-region QPs. The local wrapper starts with800 nearby inequalities, checks every omitted row, adds violated rows and repeats. The final full-grid line search still enforces every original inequality. A1e-9 normalized inward QP cushion avoids rounding-induced zero-length steps. A cutting-plane L1/HiGHS feasibility projection similarly checks the complete constraint set before returning.

The optional persistent active-set mode was used for the final255-degree continuation. It changes search efficiency, not the final feasibility check. The maximum step counts were finite budgets, generally40; no global optimization or convergence claim is made.

Every final-run JSON records cache load, guard construction, projection, optimizer, and total seconds. The selected continuation stages are separately preserved. In particular:

- The first full-constraint180-degree accepted step was recorded at49.584952seconds, including0.140222seconds cache load,0.570367 guard construction and2.899771 projection. Its subsequent continuations are in `angle180_cut2.json`, `angle180_refined.json`, `angle180_final.json`.
- The final180-degree continuation took54.417815seconds; the selected210-degree continuation took59.042783seconds.
- The complete120- and240-degree final runs took277.587758 and484.510seconds respectively.
- The255-degree first-stage checkpoint records503.923037seconds; the final28-step continuation took40.628512seconds.
- An initial full-row150-degree LP was stopped while searching. A first constraint-generation run had a JSON circular-reference logging error after saving an accepted NPZ; a first persistent-set restart had a local-variable naming error before optimization. Logs are retained. These development failures are not hidden inside a purported exact campaign-total timing; their full elapsed costs were not all instrumented.

The driver evolved during these experiments. Saved fixed witnesses and their fresh audits are authoritative. The current driver reproduces the mathematical search problem but is not claimed to reproduce each historic solver trajectory bit for bit.

## Dense audit and component certificates

`audit_left.py` uses the independent `angle_audit` assembler at spectral n96 and source/angular order80. It adds source-region quadrature n48 and entire-UV quadrature n96 at spins42..100. `batch_audit_left.py` retains the identical dense kernel bank in memory across witnesses; later head/source timings near.02/.01seconds are cache reuse, not fresh construction. The UV contraction remains witness-specific.

All ten final points pass dense head threshold/high-p screens, exact whole-domain Roy Bernstein checks, and sampled added-spin source/UV sign checks. No positive-V added-spin contribution was found. The largest head change from training among these ten is about0.0003703, unlike the separate rejected upper-g2 trial; training costs are not used in the result table.

The180-degree head was additionally replayed at n160/source96:7.7330645781685865 versus7.733064719720818 at n96/source80, a difference of-1.4155e-7. This numerical agreement is not a rigorous quadrature error bound.

Each witness has an exact all-energy crossing-series sign component and Arb256/N120 threshold-limit component for even spins2..40. For a general normal(a,b), the reviewed old threshold formula changes only by starting its effective subtraction coefficient at a and adding b to V/p. These are limit statements, not finite-width threshold neighborhoods.

For the new180-degree witness, `angle180_joint_certificate_256.json` additionally proves the limiting normalized imaginary kernel<=-1/50 for all tau>=0 in the simultaneous large-spin/large-energy limit. It has56680 compact cells and took370.053693seconds. This limit theorem is not a uniform finite-spin bound or an integral certificate. Exact source moments are cached under the new left folder. The new competitive right-sector fifteen-direction certificate belongs to the separate150-row witness; it must not be applied to these220-row candidates.

## Reproduction

Run from the workspace root with the supplied NumPy/SciPy/SymPy/CVXPY/CLARABEL and python-flint dependencies. Preserve the archived outputs and use new output names for fresh runs.

```powershell
$Python = 'C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $Python research/full_direction_20261007/left/search_left.py --angle 180 --warm output/lower_g0_upper_g2_20261006/lower_bernstein_partial_frozen.npz --constraint-generation --persist-constraints --max-iter 40 --out output/full_direction_20261007/left/reproduction180
& $Python research/full_direction_20261007/left/audit_left.py output/full_direction_20261007/left/reproduction180.npz --out output/full_direction_20261007/left/reproduction180_audit
& $Python research/full_direction_20261007/left/threshold_component.py output/full_direction_20261007/left/reproduction180.npz --out output/full_direction_20261007/left/reproduction180_threshold.json
& $Python research/lower_g0_upper_g2_20261006/matrix_crossing_series.py output/full_direction_20261007/left/reproduction180.npz --out output/full_direction_20261007/left/reproduction180_crossing.json
```

For a fixed archived witness, skip optimization and audit its NPZ directly. The exact selected names are in `direct_supports.csv`. To reuse kernels across a list of witnesses, pass all their paths to `batch_audit_left.py --out NEW_DIRECTORY`. To rebuild the joined result after all canonical audited files exist, run `summarize_left.py --require-all`; this writes only the new campaign's joined JSON/CSV. Every support is tied to its saved NPZ SHA256.

For the180-degree joint-limit proof, use `tail_candidate_joint_certificate.py` with the desired witness and an **absolute prefix inside the new output folder**, `--prec 256 --margin 1/50`. Its legacy default parent folder should not receive new outputs. These exact components remain conditional on the manuscript dispersion/null assumptions and reviewed formulas; they do not complete the omitted-domain support estimate.
