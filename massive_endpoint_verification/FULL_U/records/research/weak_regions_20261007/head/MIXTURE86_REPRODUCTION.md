# Reproduce the exact86.492742-degree mixture certificate

This construction uses a new exact rational positive combination of two already completed witnesses. It is distinct from the more strongly optimized but uncertified `angle86_refined.npz` coefficients: that file supplies only the desired normal direction.

Let the exact target normal be `(a,b)`, the stored binary64 direction in `angle86_refined.npz`. The two base objectives are `(a85,b85)` and `(0,1)`. Define exact rational positive weights

`x = a/a85`, `y = b - x*b85`.

Then the new objective is exactly `(a,b)`, and every new multiplier is exactly `x*lambda85 + y*lambda_g2`. No binary64 rounding is applied to the mixed multipliers. The final `mixture_definition.json` contains all220 rational multipliers and their exact source definitions.

For the disk support `H(U,V) = sqrt(U^2+V^2) + V`, positive homogeneity and subadditivity give

`H(x*U85+y*Ug2, x*V85+y*Vg2) <= x*H(U85,V85) + y*H(Ug2,Vg2)`.

Thus the positive sum of completed base costs certifies the mixed witness on every unrecomputed domain. Spins0 and2 are directly reintegrated because their support values improve under mixing. All other spins and both infinite-domain tails retain the proved weighted base costs. New tail-coercivity estimates are unnecessary for this specific mixture construction.

Base completed certificates:

- `output/boundary_length_20261007/complete_angle85_banded_certificate.json`.
- `output/competitive_20261007/complete_g2_certificate.json`.

The completed mixture support is assembled as

`C_mix = x*C85 + y*Cg2 - x*(H85_0+H85_2) - y*(Hg2_0+Hg2_2) + Hmix_0 + Hmix_2`.

Each subtracted term is exactly a retained additive component of its base completed certificate. The assembler reconstructs and checks all21 head components before subtraction. The85 certificate uses per-spin maxima of original/replay enclosures; the g2 certificate uses the refinedspin0 replay plus the other original replays. Using producer values indiscriminately would give an incorrect exact decomposition.

The two bases have identical exact threshold/source gap intervals for each recomputed spin. Their gap costs combine by subadditivity. Ordinary cells use the union of both saved partitions. Direct interval H, weighted convex-endpoint Taylor bounds with explicit source errors, and negative-V normalized estimates use the frozen mathematical evaluators. The minimum of these and the conservative whole-base-cell fallback is selected. Both256-bit and384-bit calculations evaluate the full exact partition; assembly retains the larger whole-spin enclosure. Every source moment file and its UV counterpart in both moment directories is hash-checked against the bound source manifests.

Use the pinned runtime and run from the workspace root:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
```

```text
PYTHON research/weak_regions_20261007/head/mixture_head.py --spin 0 --bits 256 --out output/weak_regions_20261007/mixture86/head_l0_256.json
PYTHON research/weak_regions_20261007/head/mixture_head.py --spin 0 --bits 384 --out output/weak_regions_20261007/mixture86/head_l0_384.json
PYTHON research/weak_regions_20261007/head/mixture_head.py --spin 2 --bits 256 --out output/weak_regions_20261007/mixture86/head_l2_256.json
PYTHON research/weak_regions_20261007/head/mixture_head.py --spin 2 --bits 384 --out output/weak_regions_20261007/mixture86/head_l2_384.json
PYTHON research/weak_regions_20261007/head/repair_mixture_activity.py
PYTHON research/weak_regions_20261007/assemble_mixture86.py --out output/weak_regions_20261007/complete_mixture86_certificate.json
```

Run only into new empty output locations in a disposable reconstruction; existing outputs are intentionally protected against overwrite. Both complete base proof trees are needed, not just their top-level JSON summaries. The output certificate preserves an expanded dependency hash manifest. `REVIEW_MIXTURE_HEAD.md` states the mathematical and source-binding review. The normalized support is rounded upward using an exact rational lower bound on the target normal's norm, while the inequality retains its original exact rational normal.

This is internal directed-interval certification with exact accounting under the same physical assumptions as the base witnesses. It is not proof-assistant verification or independent external review. The separate numerical86 witness has different multipliers and remains uncertified unless its own complete proof is supplied.

A source-activity audit rejected nine raw spin2 normalized candidates outside the Roy window at each precision. The mandatory repair step recomputes these cells using actual fields and convex-endpoint bounds. Assembly accepts only the repaired spin2 records and checks that all other cells are unchanged. Raw head_l2 files are diagnostic, not complete certificates. See the review document for why the already completed base certificates are unaffected.
