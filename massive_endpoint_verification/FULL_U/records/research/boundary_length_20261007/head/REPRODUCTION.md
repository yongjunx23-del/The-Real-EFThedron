# Finite-head certificates for the 85° and 255° boundary witnesses

This directory certifies the integral over **every energy for each even spin
0, 2, …, 40**. It does not by itself certify the omitted spins 42, 44, ….
The complete directional bound must use a separately completed tail for the
same exact witness and objective. In-progress output and numerical quality
targets must not be reported as completed coupling bounds.

## Frozen input witnesses

| Nominal angle | Witness | SHA-256 |
|---|---|---|
| 85° | `output/full_direction_20261007/right/angle85_220.npz` | `58fafda2eafa3a741610175f416a009c84f9a18f309d36e74d0e4007f65fd6d8` |
| 255° | `output/full_direction_20261007/left/angle255.npz` | `5b7febef7cb467ff224cd7162a82d76a3534fd3439102f459733e513b4f31615` |

The objective is the saved binary64 direction interpreted as exact dyadic
rationals, not a newly evaluated trigonometric direction:

- 85°: `(24532160660097/281474976710656,
  4486462071114449/4503599627370496)`.
- 255°: `(-2331234710160199/9007199254740992,
  -8700286382685973/9007199254740992)`.

All 220 saved multipliers are used. No coefficients are refitted in this
campaign. The row bank comprises 39 H rows, seven higher crossing rows,
60 Roy rows, 64 tapered strip rows, 36 derivative rows, and 14 rational
rows. The exact names and parameters are supplied by the witness and the
unchanged shared kernel implementation. They are the bank documented in
`research/all_direction_certification_20261007/head/LEFT180_FINITE_HEAD_FINAL.md`;
the direction and multipliers here are those of the two witnesses above.

## Analytic and interval method

The integration variable is `beta = sqrt(1 - 4/mu)`, with the exact range
`[0,1]` corresponding to all energies `4 <= mu <= infinity`. The support
function is `H(U,V) = sqrt(U^2 + V^2) + V`, with the reviewed normalization
and measure `4(2 ell + 1)/pi d beta` after scaling the fields by
`(1 - beta^2)^2`.

The new `certify_head.py` calls the unchanged reviewed classes under
`research/all_direction_certification_20261007/head/`:

- `FiniteCentered` at spins 0 and 2;
- `FiniteFast` at spins 4 through 40.

These evaluate the exact fixed dual witness with Arb enclosures, polynomial
source moments, explicit geometric remainders, continuous Taylor bounds,
threshold caps, source endpoint logarithmic caps, and the analytic UV
expansion. They do not infer a certificate from agreement of numerical
quadrature grids. The numerical values in `budgets.json` only select
desired error budgets.

The source tables use moment order 128, UV order 32, spin extent 40,
moment index extent 74, and 256-bit Arb precision. Roy Gauss order 109 is
sufficient for degree bounds 216 in energy and 207 in angle. The strip
tables use order 103 for degree bounds 136 and 204. The polynomial moments
are therefore enclosed by degree-exact Gauss rules using interval nodes
and weights; the nonpolynomial geometric remainder is separately bounded.
Every one of the six moment/UV files is independently rebuilt in
`moments128_replay/` and compared byte for byte with the production table.
The manifest also binds each UV table to its source table hash.

The beta partition initially contains the 1/256 grid, all analytic gap
endpoints, and an outward rational UV-start knot. Every cell outside an
explicit gap is included, and each gap is bounded analytically. Adaptive
bisection uses a floating estimate only to choose the next cell. All bounds
and sums use directed interval arithmetic. Each saved checkpoint contains
a complete exact rational partition of `[0,1]`; no unsuccessful or missing
cell is set to zero. A nonfinite cell is subdivided or raises an error.

Production uses 256-bit arithmetic and at most 10,000 ordinary cells per
spin. Every final partition is recomputed at 384 bits, including every
ordinary cell and every analytic gap. The source tables remain valid
256-bit enclosures. The assembler selects the larger original/replay upper
endpoint for each spin, sums exact dyadic endpoints, and rounds the final
sum outward to nine decimal places. Replay is an internal independent run
at higher precision, not an independent mathematical implementation or a
formal proof assistant verification.

## Quality budgets

| Angle | Dense numerical head, through spin 40 | Requested total slack | Sum of per-spin targets |
|---|---:|---:|---:|
| 85° | 0.4485298646067195 | 0.00187 | 0.4503998646067195 |
| 255° | 0.993606804176465 | 0.00594 | 0.999546804176465 |

At 85°, the slacks for spins 0, 2, 4 are respectively 0.0008, 0.00025,
0.0001; every even spin 6–40 receives 0.00004. At 255°, the first three
slacks are 0.003, 0.001, 0.0005, and every even spin 6–40 receives 0.00008.
These targets are not substituted for the resulting enclosures. A final
record that misses its target is still usable only at its actual bound.

## Reproduction commands

Run from the workspace root. Always use a new output directory when
reproducing a frozen record. The example below reproduces 85° spin 2;
replace the witness, output base, spin, model, and target with the exact
entries in the corresponding `budgets.json` for other components.

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$py = 'C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$witness = 'output/full_direction_20261007/right/angle85_220.npz'
$run = 'output/boundary_length_20261007/reproduction_angle85'

& $py research/all_direction_certification_20261007/head/build_moments.py --witness $witness --out "$run/moments128" --lmax 40 --N 128 --Nuv 32 --bits 256
& $py research/all_direction_certification_20261007/head/build_moments.py --witness $witness --out "$run/moments128_replay" --lmax 40 --N 128 --Nuv 32 --bits 256
& $py research/boundary_length_20261007/head/certify_head.py --witness $witness --moments "$run/moments128" --out "$run/spin2/head_l2.json" --spin 2 --model centered --target .02563318705806788 --max-cells 10000 --bits 256
& $py research/boundary_length_20261007/head/certify_head.py --witness $witness --moments "$run/moments128" --out "$run/spin2/replay_l2.json" --replay "$run/spin2/head_l2.json" --bits 384
& $py research/boundary_length_20261007/head/assemble.py --witness $witness --base $run
```

The last command deliberately reports an incomplete head when only spin 2
has been reproduced. Repeat the production/replay pair for every even spin
0–40 before expecting `COMPLETE_REPLAYED_FINITE_HEAD`.

`assemble.py` checks witness/objective hashes, current dependency hashes,
source metadata and six rebuilds, immutable `final:true` production records,
exact complete partitions, cell/gap sums, replay linkage, and strictly
higher replay precision. Its `head_summary.json` includes every selected
record/replay path and hash, exact sums, missing components, and scope.

## Scheduling and records

The original scheduler was stopped without stopping any proof worker when
Windows committed memory became scarce. Its records are retained. The
replacement `run_campaign_limited.py` set all three thread environment
variables explicitly for every new child and limited the combined campaigns
to six head workers. Once committed-memory capacity recovered, a new
`run_campaign_resumed.py` raised the limit to ten total, six per angle, with
the same explicit thread limits. Both prioritize 255°. The current live
records are `output/boundary_length_20261007/head/campaign_resumed_status.json`
and `campaign_resumed_jobs.json`; the latter records complete commands, PIDs,
logs, and scheduler-relative start/finish times. The earlier
`campaign_limited_*` records remain intact. The resumed manager inherited
the earlier active tags even when a child's first checkpoint had not yet
been written, preventing duplicate launches. Early workers were launched
before either scheduler and record their own elapsed times in proof files.

Production checkpoints are mutable while `final:false`. Only final records
are scheduled for replay. A scheduler crash does not invalidate an already
written enclosure, but a retained nonfinal file is not a completed campaign.
Per-job elapsed times overlap and must not be summed as wall-clock time.

The independent orchestration review is `REVIEW_CERTIFY_HEAD.md`. This
campaign does not modify any earlier g2, negative-g0, objective-family,
or tail proof artifact.

## Completed finite-head packets

Both campaigns completed all 21 production records, all 21 higher-precision
replays, and all six source-file rebuild comparisons. All requested quality
targets passed. These are the selected **finite-head** packets; a complete
coupling inequality still requires the separate same-witness tail.

| Angle | Retained head upper, approximately | Outward ceiling | Ordinary intervals | Analytic gap caps |
|---|---:|---:|---:|---:|
| 85° | 0.4503946446346439 | 0.450394645 | 35,112 | 105 |
| 255° | 0.9995372645459412 | 0.999537265 | 43,918 | 105 |

The authoritative exact dyadic sums and every record/replay hash are in:

- `output/boundary_length_20261007/head/angle85/head_summary.json`, SHA-256
  `0712fa77877ff76ab01a83adaeef979c72e2546365802046d48f6b2c186070f5`.
- `output/boundary_length_20261007/head/angle255/head_summary.json`, SHA-256
  `55efc63dbb352d40ec349067bf270ae2ec579ea8f950ba896f292e604b21b6c2`.

The long 85° spin-2 replay was additionally divided into four disjoint
index chunks by `parallel_replay.py`. Its combined result covers the same
7,965 cells and five gaps, agrees with the serial replay through more than
100 decimal places, and is preserved separately as
`angle85/spin2/replay_l2_parallel.json`. The selected summary above retains
the original serial replay. The parallel check took a 232.875 s wall-time
envelope; the serial replay took 589.590 s. They ran concurrently for part
of that period. `REVIEW_PARALLEL_REPLAY.md` reviews the exact coverage and
summation. No overwritten or partially replayed record was selected.
