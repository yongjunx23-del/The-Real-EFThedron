# 259.2941° finite-head campaign (directory alias angle259)

The frozen witness is
`output/boundary_length_20261007/search/angle259p2941_from255.npz`, SHA-256
`321695c10783e2e0f6e6e28589f132fdb7a68497314b1cf92eb71a07612c704f`.
The saved exact objective is

`(-6692978864176391/36028797018963968,
-4425208767245823/4503599627370496)`.

The directory alias `angle259` is only a name; the integration uses this
saved direction and never substitutes the trigonometric direction at 259°.
All 220 exact dyadic multipliers are unchanged.

The dense numerical head reference is 0.4548106193881146, from
`output/boundary_length_20261007/head384_259/angle259p2941_from255.json`
(384 nodes per ordinary velocity panel, source/angular order 128). For
spins 0–24 the cache contains 3,457 energy points, including one zero-weight
threshold probe; for spins 26–40 it contains 3,092 points, using only 20
nodes on the first velocity panel. The cache and its metadata are
`output/full_direction_20261007/right/dense220_n384k128_cache.npz` and `.json`.
Its SHA-256 is `315f478de65c5d87da9bc3db646fd4ec061191922d60bb5547a8ef97dc6edbca`.
Its role is to
choose error budgets, not establish the certificate. The requested total
head slack is 0.006: spin 0 receives 0.003, spin 2 receives 0.0015, spin 4
receives 0.0006, and each even spin 6–40 receives 0.00005. The intended sum
of per-spin targets is 0.4608106193881146; actual retained bounds must be used.

`prepare_angle259.py` produced N128/UV32 source tables at 256-bit precision
and independently rebuilt all six moment/UV files. All hashes match. The
entire source preparation and initial cap screen took 77.786 s while the
265° campaign was finishing. It waited for a free slot before starting.

The original quadratic-cap preflight is retained as failed for spin 2.
The already reviewed `certify_head_constant.py` wrapper resolves this by
using the unchanged `FiniteSpin.threshold_cap` with unchanged centered
ordinary-cell evaluation. Its 256- and 384-bit checks pass over the entire
beta interval `[0,1/256]`; the threshold support upper is at most
7.591759489543e-7. All other screened spins 0,4,6,8,14,40 pass their original
cap method. `budgets_constant.json` changes only the spin-2 model key.

The exact spin-2 V/p enclosure has endpoints approximately
`[-0.050314321505333374, -0.01671988898558457]`, strictly below zero.
Arb's default string prints the wider coarse representation `[+/-0.0504]`.
To avoid ambiguity, `constant_cap_exact_sign.json` preserves the actual
directed dyadic endpoints at both precisions. Its supplemental observation
script reads the local interval at the unchanged cap method's return;
it neither modifies the interval nor substitutes a new formula.

Production and replay are scheduled by `run_campaign259.py`, using the
unchanged constant-cap wrapper and generic integration/assembly schema.
All worker environments explicitly set OPENBLAS_NUM_THREADS,
OMP_NUM_THREADS, and MKL_NUM_THREADS to 1. The scheduler reads the preceding
campaign's active-worker count and enforces a total of ten head workers.
The records are `campaign259_status.json` and `campaign259_jobs.json` under
`output/boundary_length_20261007/head/`. Standard source, exact partition,
dependency, and 256/384-bit replay checks from `REPRODUCTION.md` apply.

This is a finite-head component campaign. A complete same-witness proof
for every even spin 42 and above is separate. No numerical head estimate,
preflight, incomplete checkpoint, or requested target is a complete bound.

## Frozen completed finite-head packet

The campaign is complete. The selected packet is
`output/boundary_length_20261007/head/angle259/head_summary.json`, SHA-256
`9db1acf65efd7e2c6aa99491257e2ec301472f6fa034b6bd28b2836c0194ffc7`.
Its status is `COMPLETE_REPLAYED_FINITE_HEAD`; it includes every even spin
0 through 40, all 21 full 384-bit replays, and all six independent source
file rebuilds. No production record is an unfinished checkpoint.

The exact retained sum, taking the larger production/replay upper bound
for each spin, is

`874200249057541354479767200908094704335567179366139075589042277330774729862950855/1897137590064188545819787018382342682267975428761855001222473056385648716020711424`.

This is approximately **0.4607996033793012**, with outward decimal ceiling
**0.460799604**. All individual quality targets were met. The rigorous
finite-head overhead relative to the stated numerical reference is about
0.005988984; the numerical reference is not used as proof evidence.

The partition contains **40,643 ordinary cells and 105 analytic gap
intervals** across the 21 spins. The gap intervals cover the threshold and
four finite source endpoints for each spin. Each production partition
covers the full beta interval `[0,1]`; each higher-precision replay uses the
same exact partition. The head summary binds the selected production and
replay files, their hashes, source inputs, exact directions, and code
dependencies. The constant-cap sign supplement remains separately bound
by the complete-certificate assembler.

The scheduler reported completion after **1,136.2186 seconds** (about
18 minutes 56 seconds), excluding the earlier source preparation. The sum
of individual production elapsed times is 5,671.9518 seconds and the sum
of replay elapsed times is 3,165.7938 seconds. These sums overlap across
parallel workers and must not be described as campaign wall time. The
scheduler exited successfully and has no remaining active workers.

This frozen packet proves the all-energy contribution from spins 0–40.
The separate same-witness proofs for spins 42 and above must be included
before reporting a complete directional bound.
