# Finite-head refinement at 261.3814° and 268°

These campaigns retain the exact saved 220-row dual witnesses and certify
their all-energy integrals for every even spin 0 through 40. The complete
directional bound additionally requires the separately computed same-witness
tail for every even spin 42 and above. Numerical references and intermediate
checkpoints are not complete coupling bounds.

## Frozen inputs and directions

The two witnesses are under `output/boundary_length_20261007/search/`:

| Alias | File | SHA-256 |
|---|---|---|
| 261p3814 | `angle261p3814.npz` | `5bd1b99e5fd4e237fe7e98583aa683ab076d4d76319b05b06c023e583f48383e` |
| 268 | `angle268.npz` | `1a0754b19846867214d914c595cdc744f433abfcdaff309d645677bbe2602044` |

The objective directions are the saved binary64 values interpreted as exact
dyadic rationals; the integration does not reevaluate trigonometric functions:

- 261p3814: `(-5399142818899535/36028797018963968,
  -8905488025140847/9007199254740992)`.
- 268: `(-5029547531033609/144115188075855872,
  -9001712312321383/9007199254740992)`.

All 220 multipliers and their row names are retained unchanged. The bank and
all-energy interval method are documented in `REPRODUCTION.md` in this same
directory. The source/kernel dependencies are unchanged from the completed
85°, 255°, 265°, and 259.2941° certificates.

## Numerical references and quality targets

| Alias | Numerical reference JSON under `output/boundary_length_20261007/` | Numerical head | Total requested slack |
|---|---|---:|---:|
| 261p3814 | `head384/angle261p3814.json` | 0.29865669703701664 | 0.0015 |
| 268 | `head384_268/angle268.json` | 0.022969543923764126 | 0.0008 |

The numerical cache is
`output/full_direction_20261007/right/dense220_n384k128_cache.npz`, SHA-256
`315f478de65c5d87da9bc3db646fd4ec061191922d60bb5547a8ef97dc6edbca`.
It uses source/angular order128 and 384 nodes on ordinary velocity panels;
spins26–40 use only20 nodes on the first panel. Numerical integration is
used solely to choose targets, not as certification evidence.

For 261p3814, spins 0,2,4 receive slacks 0.0008,0.00022,0.00021;
each of the remaining 18 even spins receives 0.000015. For 268 the corresponding
slacks are 0.0004,0.00014,0.00008 and 0.00001. The actual interval upper bound
is always used, even if a quality target is missed.

## Source rebuilds and threshold preflight

Each campaign has its own `moments128/` and `moments128_replay/` directories
under `output/boundary_length_20261007/head/angleALIAS/`. The six production
moment/UV files match their independent rebuilds byte for byte. The tables
use N128, UV order32, spin extent40, moment index extent74, and256-bit Arb
precision. The degree-exact interval source quadratures and explicit
nonpolynomial remainder bounds are the unchanged method in `REPRODUCTION.md`.
Preparation commands, logs, hashes, and elapsed times are saved in each
campaign's `preparation.json`.

The initial screened spins are0,2,4,6,8,14,40, including all five analytic
gaps at each listed spin. Every production worker independently constructs
and bounds its own gaps, so the other spins do not rely on this preflight
sample. The production and full replay cover every even spin through40.

The 268 preflight passes without a cap change. Its spin4 constant coefficient
is close to zero but strictly negative: the directed upper is approximately
−9.609690618747679e−8. The beta-squared coefficient upper is approximately
−31.054641796935606, giving a valid finite-width quadratic cap on
`[0,1/256]`; the cap contribution upper is approximately2.246521e−12.

The 261p3814 initial preflight fails only the spin2 sufficient quadratic-cap
criterion. The selected model reuses the already reviewed
`certify_head_constant.py` and unchanged `FiniteSpin.threshold_cap` with
unchanged centered ordinary-cell evaluation. Its fresh 256/384-bit paired
checks pass over `[0,1/256]`. The exact directed V/p enclosure is approximately
`[-0.14627930537228234,-0.11574923501013541]`, strictly negative.
`constant_cap_exact_sign.json` records both endpoints as exact dyadics.
The diagnostic profile hook observes the local interval at return from the
unchanged method and never changes that interval or formula.

The initial failure is preserved in `preflight.json`; it is not overwritten.
The selected row models are in `budgets_selected.json`, and
`model_selection.json` binds the original preflight, fresh cap evidence,
selected budget file, and dependency hashes. No mathematical kernel or cap
implementation was edited for these campaigns.

The independent orchestration review is `REVIEW_LOWER_REFINEMENT.md` in
this directory. It confirms that the selector and scheduler reuse the same
reviewed cap and driver, and that the endpoint observer only reads the
actual interval. The complete-certificate assembly binds this review and
the fresh selected-cap evidence separately from the unchanged head driver.

## Production, replay, and stored data

New orchestration scripts in this research directory:

- `prepare_lower_refinement.py`: fixed-witness verification, fresh source
  builds, independent source rebuilds, initial preflight, and numerical budgets.
- `select_lower_refinement_models.py`: selection of existing methods and
  fresh paired cap checks with exact negative endpoints where required.
- `run_lower_refinement.py`: up to10 single-thread workers, at most5 per
  direction; production then full384-bit replay, using the unchanged
  `certify_head_constant.py` and unchanged `assemble.py`.

The scheduler records complete commands, worker PIDs, log paths, elapsed
times, and exit codes in `output/boundary_length_20261007/head/lower_refinement_jobs.json`.
Its status file is `lower_refinement_status.json` in the same directory.
All workers explicitly set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`,
and `MKL_NUM_THREADS=1`.

Each `spinELL/` directory stores `head_lELL.json`, `head_lELL.log`,
`replay_lELL.json`, and `replay_lELL.log`. Production uses256-bit arithmetic
and an initial10000-cell quality budget. Every final partition is fully
recomputed at384 bits. The generic assembler checks exact partition
coverage, source rebuilds, immutable final records, all dependency hashes,
and strictly higher replay precision; it takes the larger production/replay
upper per spin before exact summation and outward decimal rounding.

## Reproduction commands

To reproduce a complete campaign from scratch, run the following in a fresh
copy of the workspace where the respective output directory does not yet
exist. The preparation scripts refuse to overwrite a frozen campaign.
The current environment loads **python-flint 0.9.0** from
`tmp/certification_deps/`; this version is recorded by that directory's
`python_flint-0.9.0.dist-info/METADATA`. Preserve the dependency setup from
the accompanying full reproduction bundle when replaying the interval code.

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$py = 'C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'

& $py research/boundary_length_20261007/head/prepare_lower_refinement.py --campaign 261p3814
& $py research/boundary_length_20261007/head/prepare_lower_refinement.py --campaign 268
& $py research/boundary_length_20261007/head/select_lower_refinement_models.py --campaign 261p3814
& $py research/boundary_length_20261007/head/select_lower_refinement_models.py --campaign 268
& $py research/boundary_length_20261007/head/run_lower_refinement.py --total 10 --per-angle 5
```

To replay an existing frozen partition without regenerating or overwriting
it, choose a new output path. For example:

```powershell
& $py research/boundary_length_20261007/head/certify_head_constant.py `
  --witness output/boundary_length_20261007/search/angle261p3814.npz `
  --moments output/boundary_length_20261007/head/angle261p3814/moments128 `
  --replay output/boundary_length_20261007/head/angle261p3814/spin2/head_l2.json `
  --out output/boundary_length_20261007/independent_replay_261_spin2.json `
  --bits 384
```

The driver reads the saved model key from the production record. Replace
the witness, moments path, spin, and output path for other records. Preserve
the exact witness/source/code bytes and the source manifest; hash mismatch
is an error, not something to suppress. A second384-bit replay is a repeat
of the same reviewed implementation, not formal proof-assistant verification
or an independent mathematical implementation.

## Sharper 268° spin-2 threshold cap and final selection

The original 268° preflight was valid. During integration, however, the
quadratic threshold estimate alone cost 0.00101997609802265, already above
the full spin-2 target 0.0008304318542617719. Ordinary-cell refinement cannot
remove that fixed cap cost. The original spin-2 process was therefore stopped,
preserving its incomplete `spin2/head_l2.json`, complete log, command, and
checkpoint hash. It is not a selected final record. The original scheduler
was also stopped; other live proof workers continued unchanged.

The replacement reuses the same reviewed constant-cap wrapper, with no
kernel or cap algebra changes. Fresh 256/384-bit checks select the exact
threshold interval `[0,1/1024]`. The directed V/p enclosure is approximately
`[-0.0018475906243137763,-0.0007586945746906315]`; its exact dyadic endpoints
are saved. The threshold cost is approximately 2.67e-9. Its production and
replay paths are `angle268/spin2_constant/head_l2.json` and
`angle268/spin2_constant/replay_l2.json`.

The final selection is `angle268/model_selection_final.json`, with
`selected_budgets_path` pointing to `angle268/budgets_final.json`.
It binds `constant_cap_preflight_final.json` and
`constant_cap_exact_sign_final.json`. All original model selection, budget,
and preflight files remain available. The replacement was chosen for
quality, not because the original sign test failed. The initial generated
final-selection metadata used the old selector's failure-field convention;
`finalize_lower_model_selection.py` retains that draft as
`model_selection_final_initial_metadata.json` and records the correct
distinction in the selected final metadata. No interval result is changed.

`run_lower_refinement_final.py` adopts surviving workers through read-only
Windows process handles and records new work under
`lower_refinement_final_jobs.json` and `lower_refinement_final_status.json`.
It allows either direction to use available slots while never exceeding
ten simultaneous proof workers. The original scheduling and stop evidence
are retained in the original job/status files and
`lower_refinement_handoff.json`. `adopt_worker.py` only monitors process exit
status; it does not evaluate a kernel or modify a proof worker.

`assemble_lower_final.py` is a separate copy of the reviewed generic
assembler. It substitutes only the 268° spin-2 selected directory and
output-summary filename; every source/hash/coverage/replay/sum check remains.
The selected 268° packet is `head_summary_final.json`. The 261.3814° packet
keeps `head_summary.json`. Every selected record path is explicit in these
packets, so the abandoned quadratic-cap checkpoint cannot be selected
implicitly by a glob.

For the final method from a fresh workspace, run the two preparation and
original model-selection commands above, then use these commands instead
of the original scheduler:

```powershell
& $py research/boundary_length_20261007/head/select_lower_refinement_models_final.py --campaign 268
& $py research/boundary_length_20261007/head/finalize_lower_model_selection.py
& $py research/boundary_length_20261007/head/run_lower_refinement_final.py --campaigns 268 261p3814 --total 10 --per-angle 10
```

The final scheduler also works without any historical job file, so the
abandoned cap experiment need not be repeated when reproducing the selected
certificate. The final replay command uses `spin2_constant/head_l2.json`
for 268° spin2, still through the unchanged constant-cap driver.

## Frozen completed 268° head

The selected packet is
`output/boundary_length_20261007/head/angle268/head_summary_final.json`,
SHA-256 `d403882a2e3ec448e3fecbb017fed96eb645bfd52a43465ed74e14a7218f813e`.
Its status is `COMPLETE_REPLAYED_FINITE_HEAD`. The exact rational sum is
stored in `conservative_max_of_original_replay_upper_exact`; its approximate
value is **0.02376831517256743**, with outward decimal ceiling
**0.023768316**.

This packet includes all 21 even spins through 40, all 21 full 384-bit
replays, and six matching source rebuilds. The selected partitions contain
**22,001 ordinary cells plus 105 analytic gap intervals**. All production
quality targets were met. The spin-2 row explicitly selects the new
`spin2_constant` record and replay; the stopped quadratic-cap checkpoint
is excluded.

The sum of selected production elapsed times is 3,628.7067 seconds and the
sum of replay elapsed times is 2,250.0796 seconds. These jobs overlap, so
these sums are not wall time. They also exclude the stopped original
spin-2 attempt, which remains recorded separately. Scheduler phases and
handoff timestamps are preserved in `lower_refinement_run_context.json`,
the two job/status ledgers, and the stop record.

The packet is an all-energy finite-head bound only. The complete coupling
bound must add the separately certified same-witness infinite-spin tail.

## Frozen completed 261.3814° head

The selected packet is
`output/boundary_length_20261007/head/angle261p3814/head_summary.json`,
SHA-256 `71ea32d7de96403efcff92e13a7236bfe892c287695ddfdc78ba37c619ffa02f`.
Its status is `COMPLETE_REPLAYED_FINITE_HEAD`. The retained exact sum is

`284717582548998713024733620566418530998600158596136666478607718939850045982489445/948568795032094272909893509191171341133987714380927500611236528192824358010355712`.

This is approximately **0.3001549113149621**, with outward decimal ceiling
**0.300154912**. The packet contains all 21 even spins through 40, all 21
full 384-bit replays, and six matching source rebuilds. The partitions
contain **39,740 ordinary cells plus 105 analytic gap intervals**. All
production quality targets were met. The larger production/replay upper
is retained separately for each spin before exact summation.

The sum of selected production elapsed times is 8,001.3111 seconds and the
sum of replay elapsed times is 4,468.0306 seconds. These parallel sums are
not wall time. The two-direction final scheduler phase lasted 1,121.6190
seconds and exited successfully. Including the preceding scheduler phase
and handoff interval gives approximately 35 minutes 49 seconds from the
first scheduler's process start to final completion, excluding source
preparation and earlier setup. Surviving proof workers continued computing
during the scheduler handoff. The exact process timestamps and phase ledgers
are preserved in the run-context and job/status files.

Both new finite-head packets are now frozen. No proof workers remain
active. Same-witness all-spin tail assembly is a separate required step
before either number is presented as a complete coupling bound.
