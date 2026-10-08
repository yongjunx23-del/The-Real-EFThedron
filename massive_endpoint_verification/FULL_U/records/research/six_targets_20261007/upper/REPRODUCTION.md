# Upper55 and upper95 campaign: exact witnesses and complete checks

The labels55 and 95 mean percentages of normalized upper-boundary length, with each axis divided by its displayed range. They are not angles. The corresponding physical normal angles are84.11888521120373 and 82.32092080175343 degrees. The exact binary64 normal vectors are recorded in `output/six_targets_20261007/target_map.json`.

All earlier sealed files are inputs only. New research and output files belong to `six_targets_20261007`. A running numerical probe or component proof is not a complete directional certificate. Only a final `PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND` packet supplies a new certified half-plane.

## Exact mixtures

Upper55 is an exact positive combination of the completed83.1759644845-degree and 85-degree witnesses. Solve the exact rational two-by-two linear system for the weights so their objective vectors sum to the requested target normal. All 220 multiplier sums are exact rational numbers; they are not rounded into a new floating-point witness.

Upper95 has a separate complete-mixture fallback using the completed83.1759644845-degree witness and the completed b6 witness. The latter is the original 150-row witness with the exact increment107/5 in the X4 multiplier and objective(1,6). The JSON definition retains each base's complete coefficient list, its exact weight, and the separate X4 increment. Repeated physical rows across the two bases add linearly; they are not additional independent physical assumptions.

Upper55 reintegrates spins 0 and 2. The selected upper95 mixture fallback reintegrates only spin 0; its attempted spin 2 integrations are excluded, as explained below. Every untouched spin and every omitted-spin energy domain inherit the positive weighted sum of the already complete base costs, using positive homogeneity and subadditivity of the disk support. Before subtracting any base head contributions, the assemblers reconstruct the exact original additive decomposition. The selected new contribution is the smaller of the old convex component and the larger of the 256-bit and 384-bit direct integrations. The untouched remainder is never subtracted or optimized numerically.

Upper55 has unequal spin 2 threshold caps. Its83 base has a gap ending at beta=1/1024; the 85 base ends at 1/256. Exactly64 already certified83 cells fill the difference. Their nonnegative bounds are moved into the 83 gap sum, retaining the same overall base bound and making both gap geometries identical. No new cap formula is evaluated over the larger interval.

For the 83+b6 mixture, all source and threshold gap intervals are unioned. For each base, its old aggregate gap bound plus every whole nonnegative ordinary-cell bound whose interior intersects the union bounds the enlarged gaps. Ordinary cells outside the union are integrated directly. Whole-cell bounds can be conservatively counted again on their remaining ordinary portions; no unproved proportional allocation is used. Both precisions use the same exact partition. The final assembler independently reconstructs these gap sums and every selected ordinary convex fallback.

All direct field and derivative evaluations use the existing physical kernels. The normalized spin 2 shortcut is permitted only when beta squared is at most79/119, equivalent to mu at most119/10, where the Roy source remains active. The 83+b6 spin 0 mixture uses no normalized shortcut. The legacy N64 source cache was independently rebuilt from exact rational coefficients and matched in every nested field before its original bytes were read by the new integration.

Commands, all from the workspace root:

```text
PYTHON research/six_targets_20261007/upper/probe_strategies.py
PYTHON research/six_targets_20261007/upper/prepare_upper55.py
PYTHON research/six_targets_20261007/upper/prepare_upper55_assembly.py
PYTHON research/six_targets_20261007/upper/run_upper55.py
PYTHON research/six_targets_20261007/upper/assemble_upper55.py --out output/six_targets_20261007/complete_upper55_certificate.json
PYTHON research/six_targets_20261007/upper/prepare_upper95.py
PYTHON research/six_targets_20261007/upper/run_upper95.py
PYTHON research/six_targets_20261007/upper/prepare_upper95_spin0_assembly.py
PYTHON research/six_targets_20261007/upper/assemble_upper95_spin0.py --out output/six_targets_20261007/complete_upper95_mixture_certificate.json
```

Run regenerating commands in an empty disposable reconstruction. Existing outputs are deliberately preserved. The original `run_upper95.py` attempts all four jobs and is expected to exit after the historical spin 2 exceptions; the selected complete fallback reads only its successful spin 0 records. To reproduce the selected bound without repeating failed diagnostic work, run its two saved spin 0 commands directly and retain their command records, then run the spin 0 selection preparation and assembly. Scripts ending in `prepare` or `adapt` record literal source transformations; supplied generated sources and peer review records are the authoritative versions. Per-job command records include runtime, exit code and log path. Use one OpenBLAS, OMP and MKL thread per process and `PYTHONDONTWRITEBYTECODE=1` so imported old source trees remain unchanged.

## Stronger directly optimized upper95 witness

The new 220-row search uses precisely the physical rows and robust search restrictions documented in `research/weak_regions_20261007/upper/UPPER_NUMERICAL_RECHECK.md`. They include the harmonic, crossing, Roy, tapered-strip, transfer-derivative and resummed rows; threshold samples; dense spin 10..40 guards; finite-spin/UV sampling guards; exact-polynomial crossing optimization guards; and the Roy Bernstein guard bank. There are no omitted extra physical constraints. These optimizer guards restrict the search but do not certify a continuum bound.

The 20-iteration optimization used the saved n192/order96 cache, with source margin0.0001, tail margin0.01, dense guard margin0.002 and the same spin 2/4 near-threshold guard locations. The complete invocation, environment and 245.09376549720764-second elapsed time are in `output/six_targets_20261007/upper95_search/command.json`. The numerical optimizer stopped at its iteration limit; no optimality claim is made.

The frozen witness is `output/six_targets_20261007/upper95_direct/witness.npz`, SHA256 `767b9d4aa9247f505004c392031f60ed82a3aedc7b1354e7a21ee9657e7a24ec`. Its bytes exactly match the optimizer output. Its objective is

```text
(1203580983497459/9007199254740992) g0
+ (4463211602388183/4503599627370496) g2.
```

Every exact dyadic multiplier is in `upper95_direct/exact_coefficients.json`. The 192 and 384 numerical head values are0.4235337050260135 and 0.4235337190912312. Adding the freshly sampled source and UV spins 42..100 gives 0.4235359019109998. This remains a numerical partial support until the full packet passes.

The initial one-panel Bernstein sufficient test for the higher-crossing polynomial was inconclusive. Four exact panels `[3i/32,3(i+1)/32]`, for i=0,1,2,3, prove the sign coefficient by coefficient for all z in[0,3/8] and every k>=10. Eight and 16 panels also passed. Seven exploratory tiny coefficient blends did not satisfy the coarse one-panel condition and were not selected. The full candidate is unchanged: a stronger exact polynomial enclosure resolved the issue. The original failed coarse-test record and exploratory files remain separate diagnostic evidence.

For the compact analytic bound, each coefficient's global upper bound is the maximum of its exact Bernstein upper bounds across those four panels. The original Jensen and spin-monotonicity inequalities then apply unchanged. The 400/512-bit analytic plans agree on17 adjacent energy intervals covering4<mu<=49. Their finite complements end no later than spin 132; newly contracted source moments therefore use jmax=skmax132 and kmax256. Finite compact intervals are produced at 512 bits and replayed at 640 bits. Their required common margin is kappa=1.

The fresh head rebuilds N128/Nuv32 source tables twice at 256 bits, with jmax74, and checks exact byte agreement. All even spins 0..40 are integrated over every energy and fully replayed at 384 bits. Numerical quality targets are head0 plus0.00025, head2 plus0.00008, head4 plus0.00003, head6 plus0.00002, and each remaining even spin plus0.00001; actual outward enclosures determine the bound even if a target is missed. Models are centered0/2 and fast4..40. No normalized-source applicability window is enlarged.

```text
PYTHON research/six_targets_20261007/upper/run_upper95_search.py
PYTHON research/six_targets_20261007/upper/audit_upper95.py output/six_targets_20261007/upper95_search/candidate.npz
PYTHON research/six_targets_20261007/upper/freeze_upper95_direct.py
PYTHON research/six_targets_20261007/upper/adapt_upper95_direct.py
PYTHON research/six_targets_20261007/upper95_direct/prepare_head.py
PYTHON research/six_targets_20261007/upper95_direct/run_head.py --workers 2
PYTHON research/six_targets_20261007/upper/prepare_upper95_parallel.py
PYTHON research/six_targets_20261007/upper95_direct/resume_parallel_head.py
PYTHON research/six_targets_20261007/compact_upper95/analytic_angleupper95.py --bits 400
PYTHON research/six_targets_20261007/compact_upper95/analytic_angleupper95.py --bits 512
PYTHON research/six_targets_20261007/compact_upper95/build_moments.py --kind roy
PYTHON research/six_targets_20261007/compact_upper95/build_moments.py --kind strips
PYTHON research/six_targets_20261007/upper/prepare_upper95_refinements.py
PYTHON research/six_targets_20261007/compact_upper95_v2/run_piecewiseupper95.py
PYTHON research/boundary_length_20261007/audit_compact_moments.py --witness output/six_targets_20261007/upper95_direct/witness.npz --compact-moments output/six_targets_20261007/upper95_direct/compact/moments --head-moments output/six_targets_20261007/upper95_direct/head/moments128
PYTHON research/six_targets_20261007/compact_upper95_v2/assemble_angleupper95.py
PYTHON research/six_targets_20261007/upper/adapt_upper95_assembly.py
PYTHON research/six_targets_20261007/assemble_upper95_direct_v2.py --out output/six_targets_20261007/complete_upper95_direct_global_v2_certificate.json
PYTHON research/six_targets_20261007/assemble_upper95_direct_banded_v2.py
```

The same-witness UV/real/tail preparation must precede the compact analytic commands, and its completed records must precede full assembly. Their separate `tail_upper95` code and recorded queue commands provide the finite42..448, infinite450-up, large-tau, source-remainder, real-norm and independently reconstructed local-coercivity checks. All complete packet manifests recursively bind those inputs. The global-margin packet is preserved before any tighter local-margin substitution.

The head began with two workers. Once other queues freed resources, the separately reviewed `resume_parallel_head.py` adopted the live ordinary spin 10 replay and spin 14 producer, retaining opened Windows process handles and their genuine exit statuses. It terminated only the original task-owned scheduler, preserved all completed records, and continued with four workers. Parent-process matching excludes the separate threshold-refinement process. The source witness, source moments, quality targets, max-cell counts and all producer/replay arguments are unchanged. `head/orchestration_switch.json` records the deliberate scheduler exit and both source identities; this is not a failed scientific proof. The v2 complete assembler authenticates this provenance in addition to every original mathematical check. In a fresh reconstruction the adoption script must be run while the original controller is active; it is not a sequential command after the original controller finishes. Reusing the retained final interval records and running the final verifiers needs no live process adoption.

The direction norm is bounded below by an exact 256-bit rational square-root enclosure. Reported unit-normal supports are rounded outward after division. These certificates combine directed intervals, exact rational identities, analytic infinite-domain estimates and internal peer review, with shared kernel implementations. They are not formal theorem-prover verification or an external derivation. Published dual/primal curves are graphical comparisons; missing-cap sensitivity is explicitly retained in each benchmark record.

The original compact95 finite adapter retained the 86 objective fractions in a guard assertion. Its independent point comparison stopped before emitting any finite interval proof. The failed `compact/point_check.log` and reviewed original sources remain preserved. `compact_upper95_v2` changes only that assertion to the frozen95 objective, uses a new diagnostic log, and binds the narrow correction review and adaptation ledger. All physical formulas and source moments are unchanged. The corrected point comparison passed. An independent exact moment audit checked399 factorial/Beta identities and 58,284 overlaps with the independently rebuilt head tables. Historical fractions occur only in the declared adaptation history, not in selected95 evaluators.

Both 83+b6 mixture spin 2 attempts stopped near the final ultraviolet cells because an interval-series denominator could not be proven nonzero at the requested precision. They produced no final enclosure records and contribute no partial value. The selected complete fallback uses only the two successful spin 0 records; its assembler subtracts only the exact old spin 0 component and preserves the positive combination of all remaining head and tail components. `spin0_run_manifest.json` and `spin0_assembly_adaptation.json` explicitly record this selection. Reproduce this fallback with:

```text
PYTHON research/six_targets_20261007/upper/prepare_upper95_spin0_assembly.py
PYTHON research/six_targets_20261007/upper/assemble_upper95_spin0.py --out output/six_targets_20261007/complete_upper95_mixture_certificate.json
```

For the stronger direct witness, a separate spin 2 threshold-cap refinement uses the byte-identical prior83 implementation. It preserves every ordinary cell of the original completed centered producer, replaces only its beta<=1/256 threshold treatment by the unchanged constant-centered cap at beta<=1/1024, and integrates the newly exposed interval using64 cells. All other source and endpoint gaps remain unchanged. The entire refined partition is recomputed at 384 bits. The final selected head takes the larger producer/replay value independently for each spin; it never subtracts an unverified estimate. The original ordinary global and banded certificates remain intact before the refined head is substituted in a separate complete packet.

```text
PYTHON research/six_targets_20261007/upper95_direct/run_threshold_refinement.py
PYTHON research/six_targets_20261007/review_upper95_cap.py
PYTHON research/six_targets_20261007/upper95_direct/assemble_selected_head.py --witness output/six_targets_20261007/upper95_direct/witness.npz --base output/six_targets_20261007/upper95_direct/head --spin2-directory spin2_refined
PYTHON research/six_targets_20261007/assemble_upper95_direct_refined.py
```

Completed mixture results: upper55 has unit-normal outward support 0.412936597, from exact retained components approximately 0.39299476891250973 (spin 0),0.013589425323210622 (spin 2), and 0.006352401365492951 (all remaining spins and tails). The four low-spin jobs took 1483.897987 seconds of parallel wall time; individual command records retain their own elapsed times. The upper95 spin 0-only fallback has outward support 0.434424025. Its successful producer took 229.114175 seconds including process overhead; the failed spin 2 attempts are separately preserved with exit_code1 and contribute no bound.

The direct95 original spin 2 producer exhausted its 6500-cell quality budget after 824.081428 seconds and returned a valid upper enclosure0.001743636684004362 with target_met=false. The cap refinement took 8.127383 seconds and reduced that producer to 0.0010734015175548252, with aggregate gap bound4.951110526e-9. These are component values; only the final complete packet, after the mandatory384-bit replay and all remaining spin/energy proofs, is a certified directional support.


## Final selected direct95 result and timing

The strongest complete direct95 packet is `output\six_targets_20261007\complete_upper95_direct_refined_certificate.json`, SHA256 `fc4ac8aefa201f4e651877487f45b66317c281839f5725385011acfc3211162b`, with unit-normal outward support **0.424118528**. Its exact rational components are authoritative; their decimal approximations are head_0_40=0.42405143755572683, infinite_450_up_compact=9.6434284757169749e-20, infinite_450_up_UV_and_sliver=4.0648944039449518e-08, finite_42_448_banded_all_energy=6.7048660293973035e-05. Every even spin and every energy above threshold is covered.

The preserved global precursor is 0.424914883; banding alone gives 0.424788764. The final packet additionally substitutes the fully replayed threshold-refined head. The graphical dual support is 0.42524412340588008 and graphical primal support is 0.41800770787645547; the benchmark's extraction and missing-cap caveats remain in every packet.

Recorded numerical optimization:245.093765s; numerical audits:184.556410s; head source preparation:99.918140s. The full original-to-adopted head integration/replay wall time is 2436.603760s (excluding source preparation), derived from saved process timestamps. The resumed four-worker phase took 989.607925s. Other proof branches ran concurrently, so individual elapsed times must not be added as an end-to-end wall time. `upper_consolidated_results.json` and all command records preserve the detailed stage timings.

The final assembly chain is recorded by `upper95_direct/finalize.py`; its four command/log records are under `upper95_direct` and its finalization record names the selected packet. Source and arithmetic peer reviews are linked by the proof manifests. Final mixture arithmetic was also independently checked in `upper_mixtures_peer_final_review.json` across 713 authenticated files. These are internal shared-code checks, not external or proof-assistant verification.

Final tail commands, exact coverage, paired precisions, source hashes and the two preserved metadata-reader diagnostics are detailed in `research/six_targets_20261007/upper95_direct/TAIL_REPRODUCTION_FINAL.md`.
