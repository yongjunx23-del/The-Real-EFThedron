# Upper 10% target: unchanged upper86 witness

The boundary-length target is upper10%, whose normal is approximately86.4249594898 degrees. This campaign certifies the nearby saved direction86.4927421125 degrees; the final polygon comparison evaluates the resulting half-plane at the requested boundary location.

Witness: `output/weak_regions_20261007/upper/angle86_refined.npz`, SHA256 `14e2b30df27e95ede876a8ad75ca190904547fb4ace535c14c6a11725d2df97f`. Its exact objective is `(8816243287893001/144115188075855872,8990329285288823/9007199254740992)`. All220 exact binary64 multiplier rationals are copied to `output/six_targets_20261007/upper86/head/exact_witness_definition.json`. No optimization or coefficient change is made by this campaign.

All commands below run from the reconstructed workspace root using the archived Python environment. Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`. Reproduce in a disposable extraction; preserve the archived proof files. The literal adapter generator records provenance, but final proof manifests select the supplied final files, including the reviewed compact replay extent204.

## Finite head

```text
python research/six_targets_20261007/upper86/prepare_head.py
python research/six_targets_20261007/upper86/run_head.py --workers 3
```

These build the degree-exact source tables twice and require byte-identical results. All even spins0..40 are integrated over the entire energy range. Producer256bit and fixed-partition replay384bit use the unchanged Centered model for0/2 and Fast model for4..40, with the original guarded source-activity domains. Numerical384-point values only set quality targets: additive slack0.0004 for spins0/2,0.0001 for4,0.00005 for6,and0.00002 for all others (total0.00129). Caps and ordinary cells, with their exact rational endpoints, are saved. Reaching a cell budget without reaching a quality target affects sharpness, not validity; only completed rigorous upper sums enter assembly.

## Compact all-spin coercivity

```text
python research/six_targets_20261007/tail_upper86/prepare_angleupper86.py
python research/six_targets_20261007/compact_upper86/analytic_angleupper86.py --bits 400
python research/six_targets_20261007/compact_upper86/analytic_angleupper86.py --bits 512
python research/six_targets_20261007/compact_upper86/build_moments.py --kind roy
python research/six_targets_20261007/compact_upper86/build_moments.py --kind strips
python research/six_targets_20261007/compact_upper86/run_piecewiseupper86.py
python research/boundary_length_20261007/audit_compact_moments.py --witness output/weak_regions_20261007/upper/angle86_refined.npz --compact-moments output/six_targets_20261007/upper86/compact/moments --head-moments output/six_targets_20261007/upper86/head/moments128
python research/six_targets_20261007/compact_upper86/assemble_angleupper86.py
```

The exact17 adjacent energy intervals cover4<mu≤49. Analytic Jensen/Bernstein bounds cover each interval above its saved minimum spin. Every complementary even spin from42 is checked by interval subdivision at512bits and fully replayed at640bits. The largest finite spin is204; the fresh moment source degree is204 and energy expansion order256. The normalized compact formulas are regular atmu4, with the physical lower endpoint excluded after taking its valid enclosing limit. No threshold or source endpoint is omitted.

## Ultraviolet and real-norm tail

```text
python research/six_targets_20261007/tail_upper86/certify_angleupper86_large_tau_moments.py --precision 320
python research/six_targets_20261007/tail_upper86/run_angleupper86_finite_queue.py
python research/six_targets_20261007/run_uniform_tail_queue.py --case upper86 --workers 4
python research/six_targets_20261007/run_tail_auxiliary.py --case upper86
python research/six_targets_20261007/assemble_uv_peak.py --angle upper86
python research/six_targets_20261007/assemble_uv_peak.py --angle upper86 --verify-only
python research/six_targets_20261007/banded_upper86.py --angle upper86
python research/six_targets_20261007/review_upper_tail.py
```

Finite UV spins42..448 are partitioned into eight consecutive groups; every400bit cell is replayed at480bits. The analytic eta representation covers every even spin450 and above, including its limit, using16 adjacent eta intervals and tau0..256. Its400bit partitions are fully replayed at480bits. The shifted crossing polynomial plus signed source2 moment estimate covers all tau≥256 with paired320/400bit analytic bounds and explicit geometric remainder. The overlap with low spins is checked exactly. The UV acceptance margin is1/200; actual local paired margins are retained.

The real-source-square integrals use4096 cells at256/384bits, both for finite spins42..448 and for the analytic infinite-spin remainder. The banded accounting takes the minimum valid local coercivity over every touching closed energy interval and all204 finite spins; reciprocal weights are rounded outward and the final weighted sum is exact rational. A separately implemented ledger reconstruction checks interval coverage and reproduces this cost.

## Complete assembly

```text
python research/six_targets_20261007/assemble_upper86.py --out output/six_targets_20261007/complete_upper86_global_certificate.json
python research/six_targets_20261007/assemble_upper86_banded.py
```

Both assemblers require the same witness in every component and recursively check bound source hashes, exact partitions, paired precision results and all support sums. They refuse to overwrite complete results. The banded packet preserves the first global packet and substitutes only its finite-spin tail contribution. Exact support and component fractions are authoritative; decimal bounds round outward before normalization and again afterward.

These are directed interval computations and analytic estimates using shared implementations, internally reviewed. They are not formal proof-assistant certificates, independent external verification, or optimality/attainability claims. Comparisons to published curves use the supplied PDF's digitized figures, not author-supplied exact dual/primal data.

## Selected scheduling and assembly version

During this run, the three-worker head scheduler was replaced by `resume_parallel_head.py`, using six workers and exactly the same quality targets, coefficients, kernels and per-spin cell limits. Active producer10/12 and replay8 were preserved. Opened Windows process handles retained their actual exit codes; the completed producer8's exit was separately identified as inferred from the original scheduler's checked-return-before-replay control flow. The four already completed head pairs were reused without rewriting them. `head/orchestration_switch.json` and the per-spin command records distinguish these cases.

The original scheduler's interruption and a preliminary syntax error that occurred before any action are retained as diagnostics. Neither is a failed mathematical certificate. The idle completion scheduler was separately replaced, preserving its successful independent moment audit and all original logs. No proof-producing child was stopped or restarted by these changes.

The canonical selected complete assembly therefore uses the new provenance-aware versions:

```text
python research/six_targets_20261007/assemble_upper86_v2.py --out output/six_targets_20261007/complete_upper86_global_v2_certificate.json
python research/six_targets_20261007/assemble_upper86_banded_v2.py
```

These add scheduler-ledger bindings to the same support calculation. The final canonical path remains `output/six_targets_20261007/complete_upper86_certificate.json`. `REVIEW_UPPER86_PARALLEL_ASSEMBLY.md` records the independent internal review. For a new construction from scratch, the unchanged ordinary `run_head.py` also accepts `--workers 6`; use the original ordinary-head assembler with that newly generated campaign metadata. The mathematical bounds are determined by each saved integral partition and its replay, not by scheduler width.

## Final result and observed runtime

The canonical certificate passed, with SHA256 `dd31c9e320b1559a514172226af8992ec8aaf220f8160ef9c2fb6b39dfef782f`. Its outward unit-normal support is **0.578462767** at 86.4927421125 degrees. The finite head contributes 0.5780552723853541, the finite omitted-spin tail contributes 0.0004074023123775403, and the remaining infinite-spin tail contributes approximately 0.0000000912419188564. The exact fractions in the certificate are authoritative.

The complete calculation uses 39,940 ordinary head cells over all 21 even spins from 0 through 40, with every saved cell replayed at higher precision. The compact proof has 17 adjacent energy rectangles. The independent finite-tail accounting checks 167,008 paired UV cells against all 4,096 real-norm cells and all 204 finite tail spins. Source endpoint caps, the threshold cap, the energy endpoint at infinity, and the infinite spin domain are included by the separately bound analytic estimates.

The observed interval from the original head controller launch to final complete assembly was **2,918.91 seconds (48.65 minutes)** on this shared machine. This includes concurrent jobs and scheduling gaps, and excludes the earlier numerical optimization and input preparation. It is not a single-core runtime or an independently repeated benchmark. `output/six_targets_20261007/upper86/result_and_timing.json` preserves the timestamps, per-job timings, cell counts and source bindings. Parallel job durations must not be added to estimate elapsed wall time.

Spins 2, 8, 10, 12, 14 and 16 reached their cell limits before the requested quality targets. Their rigorous upper bounds remain valid, and the actual conservative values were used. No quality target was substituted for an integral bound.

The subsequent independent internal accounting review passed: `output/six_targets_20261007/upper86_peer_final_review.json` checks the exact final sum, source hashes, head and compact partitions, infinite UV coverage, and the exact substitution of the already authenticated local tail ledger. It authenticated 302 files. This is an additional internal check, with the shared-code limitations stated above.

At this direction, the digitized published dual support is approximately 0.594477234 and the digitized primal support approximately 0.516705903. Thus this certified upper bound improves the plotted dual support while remaining above the plotted primal value. The comparison is to figure-derived data; it does not establish the precision of the authors' underlying optimization.
