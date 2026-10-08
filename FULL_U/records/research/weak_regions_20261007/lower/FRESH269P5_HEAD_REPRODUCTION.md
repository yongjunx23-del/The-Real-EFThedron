# Fresh 269.5-degree finite-head certificate

This component treats the newly optimized 220-row witness, not the shifted older 268-degree witness. The unchanged rigorous field evaluators integrate all energies for every even spin from 0 through 40. A separate same-witness certificate for every even spin at least 42 is required before this becomes a complete coupling bound.

The completed selected head has status `COMPLETE_REPLAYED_FINITE_HEAD`. Its conservative directed upper is approximately `0.00189186216831086`; its outward nine-decimal ceiling is `0.001891863`. The exact rational is in `head_summary.json`, SHA-256 `7ada1b850abd1ac767483021aaddb4e1ffdeec0cf9912428dbb0e388cf7e3513`. There are24246 ordinary cells and105 analytic gaps across21 spins; every saved cell and gap has a384-bit replay. The complete head campaign, including both source builds and all replays, took999.23 seconds. Spin2 was the only missed quality target; its valid actual upper was retained. These figures exclude optimization, numerical audits, and the separate infinite-spin tail campaign.

## Frozen witness

- File: `output/weak_regions_20261007/lower/search/angle269p5.npz`.
- SHA-256: `f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae`.
- Exact objective: `(-2515252609197481/288230376151711744) g0 + (-9006856288192973/9007199254740992) g2`.
- Every binary64 multiplier is interpreted as its exact dyadic rational. The exact multiplier table is in `output/weak_regions_20261007/lower/LOWER_REFINEMENT_RECHECK.md` and its JSON counterpart.

## Reproduce

Use the runtime and dependencies in the main reproduction archive. Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`. In a fresh checkout where the following output directory does not yet exist, run:

```text
python research/weak_regions_20261007/lower/certify_fresh269p5_head.py
```

The output is `output/weak_regions_20261007/lower/head_certification/angle269p5/`. The driver refuses to overwrite an existing campaign. For rechecking already saved results, use the exact producer and replay commands in `spin*/commands.json`, changing only each new output path to an unused filename. The final `campaign.json` and `jobs.json` retain all successful commands and their elapsed times.

## Source enclosures and integration

The unchanged `research/all_direction_certification_20261007/head/build_moments.py` creates independent `moments128` and `moments128_replay` directories at 256-bit precision, with source order128, ultraviolet order32, and spin extent40. All six source/ultraviolet moment files must match their independent rebuilds byte for byte. The moment metadata and source-degree exactness checks are also checked by the final assembly.

Every production integral uses the unchanged `research/boundary_length_20261007/head/certify_head_constant.py` at 256 bits. Spin0 uses the centered model; spin2 uses the previously reviewed constant-centered model; spins4 through40 use the existing fast model. The adaptive partition includes every ordinary cell and all five analytic gaps, including the threshold and the source endpoints. Exact rational partition endpoints must cover the full velocity interval from0 through1 without omissions or overlaps.

Each complete saved partition is replayed at 384-bit precision by a fresh process. The assembly uses the larger of the original and replay exact dyadic upper bounds for each spin. These are independent executions of the same reviewed implementation, not independently derived mathematical implementations or a formal proof assistant check.

The spin2 negative-constant cap automatically selects velocity width1/262144 for this witness. Its gap contribution is tiny; the ordinary-cell interval estimate is the slower part of its integration. The selected model uses the unchanged `FiniteSpin.threshold_cap` formula. Paired 256/384-bit preflight data for spins0,2,4,6,8,14,40 are in `preparation.json`. Every production and replay worker computes its own analytic gaps, so the other spins do not rely on sampled preflight coverage.

## Quality targets versus validity

Targets are the independent n384 numerical head integrals plus these absolute slacks:

| Spin | Requested slack |
|---:|---:|
| 0 | 0.000020 |
| 2 | 0.000006 |
| 4 | 0.000005 |
| Each remaining even spin through40 | 0.000001 |

The total requested slack is0.000049. Each producer is capped at7000 ordinary cells. A missed numerical target does not invalidate a finite rigorous upper bound: the actual directed upper value is retained and summed. Numerical values only select target precision; they are never substituted for proof contributions.

## Assembly and unused experiments

The selected `head_summary.json` is produced by the unchanged `research/boundary_length_20261007/head/assemble.py`. It checks witness and source hashes, exact directions, every dependency hash, complete partitions, higher-precision replay, source rebuilds, and exact arithmetic sums. Its complete status is `COMPLETE_REPLAYED_FINITE_HEAD`; its scope is still only spins0 through40.

A separate spin2 quadratic-cap experiment is retained under `spin2_quadratic_narrow/`, with its cancellation record. It chose a smaller cap width1/4096 using the unchanged quadratic formula, but offered little improvement in the convergence of ordinary cells. It was stopped and is excluded from the selected head. The first diagnostic report has a Python formatting exception, preserved verbatim; the corrected paired diagnostic is in `quadratic_spin2_preflight_v2.json`. Neither diagnostic is needed by the selected proof. The alternate wrapper and assembly selector are likewise unused in the final selected head.

The mathematical row bank, normalization, unitarity disk, and all220 physical/equality rows are unchanged from the full constraint catalogue. This campaign adds optimization at a new direction and tighter interval estimates, not new physical assumptions.
