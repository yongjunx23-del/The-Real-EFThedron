# Three lower-boundary mixture certificates: method and reproduction

This campaign targets the lower branch at 15%, 45%, and 60% of its displayed, normalized boundary length. These are boundary stations, not angles. The corresponding physical outward normal angles are approximately 258.5185525804567, 260.4648793723453, and 262.3334856668359 degrees. The target map records the exact binary64 coordinates that define the requested normals. Axes are scaled by their observed published-dual ranges before boundary length is measured.

This file describes the calculation. Only the separate complete certificate packets, once successfully assembled, establish bounds. Numerical probe files alone do not establish bounds.

## Exact witnesses

| Station | First completed witness | Second completed witness | Directly reintegrated spins |
|---|---:|---:|---|
| lower15 | 255 degrees | 259.294118509649 degrees | 0, 2, 4 |
| lower45 | 259.294118509649 degrees | 261.3814 degrees, as exactly stored | 0, 2, 4 |
| lower60 | 261.3814 degrees, as exactly stored | 263.27473709797 degrees | 0, 2 |

The exact objective vectors, angles, source paths, source hashes, and all 220 multipliers are stored in `output/six_targets_20261007/lower/<tag>/definition.json`. Angles in prose are descriptive; the rational objective vector in that definition is authoritative.

Let the two existing objective vectors be `d1=(a1,b1)` and `d2=(a2,b2)`, with exact rational multipliers `lambda1` and `lambda2`. For target `d=(a,b)`, define

```text
det = a1*b2 - b1*a2
w1 = (a*b2 - b*a2)/det
w2 = (a1*b - b1*a)/det
lambda = w1*lambda1 + w2*lambda2
```

The preparation and final assembly require both weights to be strictly positive and verify `w1*d1+w2*d2=d` exactly. They check all 220 coefficient identities and the complete ordered row-name lists. Binary64 base coefficients and target-normal components are interpreted as exact dyadic rationals. The mixture coefficients are exact general rationals. A rounded replacement NPZ is never used to define a certified witness.

## Constraints

No new physical constraint is introduced, and no original row is dropped. Each witness uses the identical 220-row constraint bank of its two previously certified parents. The exact ordered row names and every coefficient are included in each definition. Their interpretation and physical assumptions are inherited from the complete prior constraint catalogue and certificates; this addendum does not replace those required dependencies.

Positive combination preserves all equality identities and the required signs of multipliers on inequality constraints. It also preserves the positivity bounds used for omitted spins. The new work tightens the evaluation of the same disk-support dual objective. It is not a new independent optimization of all 220 multipliers at each target.

## Why the inherited infinite domains remain bounded

For the disk-support integrand

```text
H(U,V) = sqrt(U^2+V^2) + V,
```

positive homogeneity and subadditivity give

```text
H(w1*U1+w2*U2, w1*V1+w2*V2)
    <= w1*H(U1,V1) + w2*H(U2,V2).
```

The fields depend linearly on the objective and multipliers. Thus the positive weighted sum of the two completed certificate costs bounds the exact mixture for every energy and every even spin. This inherits all compact-energy, ultraviolet, large-spin, source-endpoint, and infinite-energy proofs from the completed parents. No new truncation of those domains is made.

To improve that valid initial result, the assembler removes only the exactly identified retained costs of the explicitly reintegrated spins. It replaces each with the smaller of (i) its paired direct interval enclosure and (ii) its original weighted convex bound. The remaining contribution is the exact weighted total minus exactly those retained spin costs. The complete parent proof packets and all their hashes are retained as dependencies.

## Direct low-spin interval enclosures

The spin0/2 evaluator is a transparent adaptation of the previously reviewed exact-mixture evaluator. It reads both original source-moment packages and checks their witness identities, manifests, saved cell records, and source hashes. The source moments are reused without changing their witness hash: fields are evaluated separately for each base and only then combined with the exact rational weights.

Both selected bases have identical threshold/source gap intervals for every directly reintegrated spin in this campaign. The code requires equality. It uses the union of both original cell endpoint lists outside those gaps, giving an exact partition of beta from zero to one. The complete gap contribution is bounded by the positive weighted sum of the saved base gap bounds.

Every new cell is bounded by the minimum of applicable valid options:

1. The positive weighted sum of the two saved whole-cell bounds containing the new cell. Charging the entire enclosing old cell, even for a smaller subcell, is conservative because the integrand is nonnegative.
2. The support function evaluated on interval enclosures of the combined fields, multiplied by the cell width and the original spin factor.
3. The same convex endpoint/Taylor remainder formula as the prior reviewed mixture calculation. Combined derivatives and source errors are formed by exact positive weighting before applying the enclosure.
4. A negative-V normalized bound, **only when the upper beta endpoint is strictly below 1/2**. The code does not extend the original Centered normalized branch past its source-activity domain.

Spins0/2 use the original `FiniteCentered` fields. Spin4 uses the original source-aware `FiniteFast` fields. The latter's optional Taylor construction can encounter a denominator interval containing zero on a wide cell. The spin4 v2 evaluator catches only `ValueError` and `ZeroDivisionError` from that optional construction and retains the other valid candidates, including the unconditional saved whole-cell fallback. It does not treat a failed candidate as zero or discard the cell.

Each complete fixed partition is independently recomputed at 256 and 384 bits. The final assembly checks both exact partitions, both full sums, matching exact definitions, the safe normalized-method condition, and every selected saved-cell fallback. It retains the larger of the two direct whole-spin bounds. The two precisions share mathematical kernels; this is a replay and exact-accounting check, not an independent implementation of the physical kernels.

## Reproduction order

Use a fresh reproduction directory containing the required prior archives, this campaign's research sources/reviews, and its saved `output/six_targets_20261007/target_map.json`, but no generated lower output directory. Existing complete outputs intentionally refuse overwrite. Use the saved runtime/environment and set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`.

```text
python research/six_targets_20261007/lower/prepare.py
python research/six_targets_20261007/lower/run_mixtures.py
python research/six_targets_20261007/lower/run_spin4_v2.py
python research/six_targets_20261007/lower/assemble.py
python research/six_targets_20261007/lower/audit_and_report.py
python research/six_targets_20261007/lower/run_spin0_refinement.py
python research/six_targets_20261007/lower/assemble_spin0_refined.py
python research/six_targets_20261007/lower/audit_selected.py
```

The independent head campaigns can execute concurrently, including the separate spin0 refinement. Their complete command arrays, exit codes, elapsed times, thread environment, source hashes, and logs are saved under `output/six_targets_20261007/lower`. The original failed spin4 first attempt is retained for transparency; reproducing that failed attempt is unnecessary. Its records are diagnostics and are not selected by the final assembler. The selected spin4 command records end in `_v2.json`. `finish_campaign.py` is an optional orchestration alternative to the assembly/audit commands: it waits for the head campaigns and source reviews, then records each individual assembly and audit command with its actual exit code and elapsed time.

The new source review is `research/six_targets_20261007/REVIEW_LOWER_MIXTURES.md`, required before complete assembly. Source adaptation ledgers allow reconstruction of the small changes from the frozen previous implementation. The three `adapt_*.py` scripts are optional source-reconstruction tools, not necessary proof commands: run them in order (`adapt_mixture_head.py`, `adapt_spin4.py`, `adapt_spin4_v2.py`) only in a separate scratch copy lacking their generated source targets. They intentionally refuse to overwrite the supplied generated files.

## Additional spin0 partition refinement

The lower15 and lower45 spin0 calculations receive a further refinement: every ordinary cell of the initial exact-mixture partition is bisected at its exact rational midpoint. Gap intervals, gap bounds, exact weights, all 220 multipliers, and every field/support formula are unchanged. The complete refined partition is evaluated at both 256 and 384 bits. These two precision jobs can run concurrently because both use the same predetermined geometry.

The separate `assemble_spin0_refined.py` requires the original complete packet. It checks exact bisection, unchanged definitions and gaps, all selected whole-cell fallback costs, both full sums, and successful recorded commands. It then replaces only the original retained spin0 upper contribution with the larger of the two refined direct bounds, requiring a strict improvement. All other exact components and both infinite domains are preserved. The original complete packets and their initial final audit are not modified.

The strongest final selected paths are `output/six_targets_20261007/lower/complete_lower15_refined_certificate.json`, `complete_lower45_refined_certificate.json`, and `complete_lower60_certificate.json`. The separate review `research/six_targets_20261007/REVIEW_LOWER_SPIN0_REFINEMENT.md` is required for the refinement assembly. `adapt_spin0_refinement.py` can reconstruct the refined evaluator in a scratch source tree; it makes only selector changes and inserts the exact cell-bisection line. The final selected arithmetic/hash audit and quantitative handoff are `final_selected_lower_review.json` and `FINAL_SELECTED_LOWER_RESULTS.md` in the lower output directory.

## Limits of the claim

These are internal computer-assisted continuum certificates under the physical assumptions of the original constraint bank. They are not proof-assistant verification, external review, or optimality certificates. Comparisons with published primal/dual curves use the saved digitized figures; the published authors' exact optimizer data and certified digitization errors are unavailable. The ordinary 192/384-node numerical probes are diagnostic only and are excluded from the proof arithmetic.
