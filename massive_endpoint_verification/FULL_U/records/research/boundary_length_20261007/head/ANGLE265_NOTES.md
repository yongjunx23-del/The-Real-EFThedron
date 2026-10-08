# 265° finite-head campaign

The fixed witness is
`output/boundary_length_20261007/search/angle265.npz`, SHA-256
`5aae4902c2db47f5e3a5f5a37d87c136cd7ec3cb191752d8f3b2f7ee04d6d8c3`.
The direction and all 220 multiplier values remain exact binary64 dyadics.
No reoptimization occurs in this head campaign.

The dense numerical head in `output/boundary_length_20261007/head384/angle265.json`
is 0.1196042268215931. This is only a quality reference. The desired rigorous
head slack totals 0.00194: spins 0, 2, 4 receive 0.0008, 0.0004, 0.0002;
each even spin 6–40 receives 0.00003. Thus the intended sum of targets is
0.1215442268215931. The actual retained interval sum, including every spin
and higher-precision replay, must be used in a complete certificate.

`prepare_angle265.py` built N128/UV32 source tables at 256 bits and independently
rebuilt all six moment/UV files. All six hashes match byte for byte. The
initial production/rebuild times were 44.25 s and 45.06 s; preparation plus
the first cap screen took 102.58 s. These elapsed times were measured while
the two earlier campaigns were running, and are not a claim about isolated
runtime or total optimization time.

The original `preflight.json` intentionally remains marked failed for spin 2.
The old `FiniteTaylor.threshold_cap` requires both a negative constant and a
negative quadratic upper coefficient. This sufficient condition does not
apply to this spin-2 witness. It must not be treated as a failure of the
constant-sign condition itself or concealed by replacing the old record.

The separately hashed `certify_head_constant.py` driver adds a
`constant_centered` model. Its only change is that a `FiniteCentered`
subclass delegates `threshold_cap` to the unchanged and previously reviewed
`FiniteSpin.threshold_cap`. Ordinary-cell evaluation remains centered.
All original models, the original driver, and the completed 85°/255°
dependencies remain unchanged.

The existing constant-sign cap encloses V/p over the whole beta interval
`[0,1/256]`, including its source remainders, and proves negativity. It then
uses `H <= U^2/(-2V)` and integrates the explicit beta power analytically.
There is no requirement on the sign of a Taylor quadratic coefficient.
The cap was independently evaluated at 256 and 384 bits:

- the V/p enclosure is displayed as `[-0.4 +/- 0.0346]`, strictly negative;
- the threshold integral upper is 1.099415773639008e-7 at 256 bits;
- the 384-bit result is 1.099415772037823e-7;
- all source endpoint gaps also have finite nonnegative upper bounds.

The full exact dyadic bounds, gap endpoints, implementation hashes, and
source manifest hash are in `constant_cap_preflight.json`. This is a cap
component only, not a full head or coupling certificate. The other screened
spins 0, 4, 6, 8, 14, 40 passed the original cap method.

`budgets_constant.json` retains the original numerical targets and changes
only the spin-2 model key to `constant_centered`. The scheduled integration
driver is `run_campaign265.py`; it uses the new wrapper for both production
and replay, explicitly sets OPENBLAS_NUM_THREADS, OMP_NUM_THREADS, and
MKL_NUM_THREADS to 1, and defaults to six workers. It must be launched only
when the ongoing priority campaigns have capacity. Its logs are
`campaign265_status.json` and `campaign265_jobs.json` under the common head
output directory. The generic assembler checks the unchanged output schema
and adds all wrapper hashes appearing in the selected proof records.

The actual campaign initially shared a ten-worker cap with the last
85°/255° replays by reading their live status. After both earlier campaigns
closed, `run_campaign265_expanded.py` inherited all six active worker tags
and filled the remaining four slots. The original manager alone was stopped;
all proof workers and earlier logs were preserved. Its current records are
`campaign265_expanded_status.json` and `campaign265_expanded_jobs.json`.
The expanded scheduler uses ten single-thread workers and does not modify
the mathematical drivers or the frozen records.

The command, after earlier head campaigns have freed capacity, is:

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
& 'C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' research/boundary_length_20261007/head/run_campaign265.py --per-angle 6 --total 6
```

The same all-energy source, partition, dependency, and 256/384-bit replay
requirements documented in `REPRODUCTION.md` apply. A separately completed
same-witness tail for all even spins 42 and higher remains necessary.

## Completed packet

All 21 production integrals, all 21 higher-precision replays, and all six
source rebuild comparisons passed. The retained finite-head upper is
0.12153957717774937 (approximately), with outward ceiling **0.121539578**.
It covers 30,914 ordinary beta intervals and 105 analytic gap caps. Exact
dyadic sums and all selected record/replay hashes are in
`output/boundary_length_20261007/head/angle265/head_summary.json`, SHA-256
`ef7c221e436388a95ddda09becded2fa5a377fe100f4ec502ae07f3bb06c10e5`.
The scheduler exited successfully. The constant-cap implementation and
selected proof dependencies are frozen.
