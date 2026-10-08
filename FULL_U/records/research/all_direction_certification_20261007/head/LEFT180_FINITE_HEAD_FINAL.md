# Final interval head for the competitive left endpoint

The selected finite-head certificate covers every even spin0,2,...,40
at every energy mu>4. Its conservative cost is at most **7.747533**.
The exact sum before outward decimal rounding is in
`output/all_direction_certification_20261007/head/left/head_refined_summary.json`.
It is the sum, over spins, of the larger of each adaptive result and its
partition replay, and displays7.7475328978105.

This is a finite-head result. The parent complete certificate adds a
same-witness all-spin tail of at most.003968804862404489 and gives
**g0 >= -7.751502**. The authoritative complete combination is
`output/competitive_20261007/left/complete_left_certificate.json`;
its final tail and assembly provenance are owned by the parent report.

## Frozen input and selected records

The exact objective is(-1,0). The unchanged witness is
`output/full_direction_20261007/left/angle180_final.npz`, SHA256
`e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
Every stored binary64 multiplier is interpreted as an exact dyadic rational.
The canonical row parameters and C=4/3 are exact rationals.

The original completed head uses21 adaptive records and21 full replays.
Its conservative max-per-spin sum is7.757833840054715, rounded up to7.757834.
It is preserved independently in`left/head_summary.json`, SHA256
`281ac1dfe687df4a25576dafb38e89ccd405662f93d48917283ec5136a582d66`.

The refined head replaces only spin0 with the separately computed and
replayed`left/low0_refined/`record. All other records and the witness are
unchanged. The refined summary SHA256 is
`a42e330cf92342de5de98e91ce0ad438fc444dd66463a19bdc051201435fb550`.
It identifies every selected adaptive/replay path and SHA256, exact dyadic
upper sums, code dependencies, source tables, and independent table rebuilds.

## Constraints retained in this220-row witness

These are the same rows used in the numerical search; certification does
not introduce an unrecorded multiplier or assume a omitted row vanishes.

| Row family | Count | Retained parameters |
|---|---:|---|
| Forward H rows |39|H2 through H40|
| Crossing nulls |7|X4,X5,X6,X7_1,X7_2,X8_1,X8_2|
| Roy source rows |60|B=119/10; J=0,2,4,6; n=0,...,14|
| Tapered strip sources |64|Energy windows[12,24]and[24,48]; T=4/3; j=0,...,3; n=0,...,7|
| Transfer derivatives |36|Derivative orders1,2,3, each with H2,...,H12 and one regularized G row|
| Resummed rational rows |14|Cosine/sine rows, r=4/5, c=-9/10,-3/4,-1/2,0,1/2,3/4,9/10|

Their exact analytic implementations are in`forward220.py`,
`research/crossing_kernels.py`, the Roy/strip moment code, and the frozen
numerical row builders. This inventory specifies which constraints are
present; it does not claim that the finite list is a complete description
of all physical amplitudes. Each added identity enlarges the available
dual-witness combinations, while the disk-support inequality gives a valid bound
for every amplitude satisfying the underlying dispersive assumptions.

## Continuous-energy enclosure

Set beta=sqrt(1-4/mu) and e=1-beta^2. The complete energy domain becomes
0<beta<1. The code evaluates the scaled fields U/e^2 and V/e^2 and
integrates the exact disk support H(U,V)=sqrt(U^2+V^2)+V with the
beta-space measure4(2ell+1)/pi. The measure/field conventions are inherited
from the reviewed base head implementation.

Each selected spin has an exact rational partition of[0,1]. Ordinary
intervals and five explicitly bounded gap intervals have disjoint
interiors and meet exactly; every full replay verifies this union.
The refined selection contains61854 ordinary cells and105 gap intervals.
The cells are an integration partition, not a grid-only feasibility test.

For ordinary cells, the code chooses among valid outward Arb enclosures:
direct interval support, convex endpoint interpolation with enclosed
second derivatives and source-error allowances, and a normalized
negative-V estimate when its sign is proved. The latter uses
H(U,V)<=U^2/(-2V) for V<0. Centered second- and fourth-order Taylor
enclosures reduce interval overestimation. Source remainders are uniform
value errors; the code never differentiates an error ball as though it
were an exact constant.

Near threshold, even-analytic quadratic estimates handle the especially
small spin2/4 margins. Separate threshold caps enclose the removed
neighborhood; no threshold interval is discarded. At finite source
endpoints, polynomial/log formulas are enclosed directly or by explicit
integrable log-gap estimates. Roy endpoint terms above Jmax=6 vanish
identically. Strip endpoint factors include the tapered zero that makes
the relevant x*log(x) bound integrable. Active and inactive source portions
are both enclosed when an interval touches an endpoint.

Principal-value approximants are selected from the full interval and held
fixed for its midpoint and range Taylor jets. This avoids combining
derivatives of one approximant with point values from another.

## Source tables, ultraviolet truncation, and arithmetic

All calculations use Arb at256bits. Finite-energy source moment order is
N=128; ultraviolet order is32. The ultraviolet expansion is used for
mu>=128 and includes beta=1 with an explicit geometric remainder.
No cutoff of the energy integral is treated as zero.

Roy source tables use enclosed Gauss-Legendre nodes/weights of order109.
The retained integrands have energy degree at most216 and angular degree
at most207, within the exactness degree217. Each strip uses order103,
with degree bounds136 and204, within the exactness degree205.
This polynomial exactness plus enclosed node/weight evaluation certifies
the retained moments; it is not an empirical quadrature-convergence claim.

Three finite moment tables and three UV tables were independently rebuilt
and match bit for bit. Their manifest is`left/moments128/manifest.json`,
SHA256`2ff5d50c715f4b2a091db0e38d90ad51326d640a62eb1412b65af01bb43e5601`.
The rebuilt tables are in`left/moments128_replay/`. Source norm bounds and
geometric remainders are retained in the table metadata and used by the
field evaluators. The complete source and UV moment formulas are in
`research/certified_roy_moments.py`,`certified_strip_moments.py`, and
`certified_uv_moments.py`.

Numerical quadrature was used only to set convenient error budgets and to
guide optimization. The numerical all-energy head estimate7.73306457817
is not substituted for any certified interval contribution.

## Replays and exact-polynomial acceleration

Every original adaptive partition has been replayed without trusting its
stored cell bounds. Spins0/2 use the centered evaluator,4/6 the quartic
evaluator, and8..40 the gap-aware evaluator for adaptive integration.
Some later replays use`FiniteFast`, which exactly contracts finite source
polynomials before interval evaluation. The independent formula review is
`../REVIEW_FINITE_FAST.md`. Before adoption,252 interval/jet comparisons
passed for spins14/40 across representative compact and UV cells, with
identical range-selected PV modes. The diagnostics are separate from proof.

Algebraically equivalent interval evaluations can have different widths.
Initial fast replays12/14/16 exceeded tight quality budgets by about
9.5e-6,7.1e-6,2.6e-6. Those attempts emitted no selected record. Spin12
was replayed with the original evaluator. Spins14/16 were rerun with
slightly looser explicit ceilings, and their actual upper values are kept.
The final conservative sum takes max(original,replay) separately for each
spin, including the small widening at14/16/18. Failure records are retained
in`left/fast_replay12_failure.json`and`left/fast_replay14_16_failures.json`.

`assemble_left180.py` verifies the original complete head.
`assemble_left180_refined.py` creates the separately selected refined head.
The latter never overwrites the former. The refined spin0 upper is
6.83685814962010641596, compared with6.84715909186432157129 previously.
The reduction is0.0103009422442151553 and its full replay passed.

## Fixed negative-g2 family

`replay_left_family.py`and the separately hashed fast/ten-slope variants
retain the same multipliers and source rows, replacing the exact objective
with(-1,-b). The identity is U_b=U_0,V_b=V_0-b*p with p>0.
Monotonicity of H in V justifies reusing the base gap and complete tail
bounds. The selected shifted interval cells are recomputed, not inferred
from a sampled curve. Parent assembly performs outward normalization by
sqrt(1+b^2).

All21 finite spins are complete for b=1/4,1/2,1,2,4,8. Their unnormalized
head upper costs are approximately7.68637540723,7.62783473620,
7.53154597174,7.39725299664,7.23181003875,7.04339864362.
The exact sums and provenance are in`left/family_packets.json`, SHA256
`306c5897c6b8c11c34646e024824d778eefa5bd8b673d4cbec6b2a96541d3f16`.
Four larger shifts remain explicitly incomplete and are not selected.
These fixed-family bounds are not independently optimized and become
weaker farther from180degrees; competitiveness must be assessed separately
from validity. See`LEFT_FAMILY_CERTIFICATION.md`and
`REVIEW_LEFT_FAMILY_ASSEMBLY.md`.

## Reproduction and timing

Run from the repository root with the configured Python containing NumPy,
SymPy and python-flint0.9.0. The scripts also add`tmp/certification_deps`
to the import path. Rebuild source tables with`build_moments.py --witness
output/full_direction_20261007/left/angle180_final.npz --out NEW_DIRECTORY
--lmax 40 --N 128 --Nuv 32 --bits 256`.

Exact adaptive commands for4..40 are in`left/launches_4_40.json`.
The original0 budget is numerical*1.002+.0001;2..40 use
numerical*1.01+.0001. The refined0 command uses`head_centered_batch.py`,
`--start 0 --stop 0 --relative-slack .0005 --absolute-slack .00005
--max-cells 12000`, with the same witness/moment directory and numerical
budget file`output/full_direction_20261007/left/head180_n160/angle180_final_audit.json`.
Always use a new output directory when recomputing adaptive artifacts.

To replay a selected ordinary component, use`replay_gap.py --witness W
--moments M --record R --out NEW_REPLAY --ceiling C`, with the exact
record and ceiling from the frozen summary/replay. Family replay variants
take`--moments M --record R --out NEW_REPLAY --family-out NEW_FAMILY
--ceiling C`; the base witness is hash-fixed inside those drivers. Select
the driver identified by the recorded dependency hashes. For fast14/16
the replay ceilings are.00845/.00636; fast18 uses.00350. Failed attempts
must not replace successful selected records.

The selected refined adaptive jobs sum to38188.272 elapsed job-seconds;
their selected replays sum to12645.372 job-seconds. These overlapping
parallel durations are neither total wall time nor measured CPU time.
They exclude discarded attempts and some source-table/diagnostic work.
The optional refined0 alone took320.464seconds to integrate and
130.443seconds to replay. Full failed-fast timings were not retained;
bookkeeping records only the last printed elapsed time and marks the
unknown total explicitly.

The proof is an interval-enclosed computation under the stated dispersive
and partial-wave assumptions, with internal independent formula review
and reproducible hash-bound records. It is not a proof-assistant
formalization or an external audit. The finite-head report alone never
certifies omitted spins; the separate complete parent tail is essential.
