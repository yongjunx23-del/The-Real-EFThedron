# Left180 fixed-witness finite-head continuation

This branch reuses the reviewed generic220-row interval head machinery,
without changing or reoptimizing the left witness.

Witness: `output/full_direction_20261007/left/angle180_final.npz`.
SHA256: `e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
Direction is exactly(-1,0). The allenergy numerical head through40 is
7.7330645781685865; this is a convergence diagnostic, not a certificate.

Exact source moments N128 and UV moments N32 for this same hash were
built in `output/all_direction_certification_20261007/head/left/moments128/`.
Their manifest SHA256 is
`2ff5d50c715f4b2a091db0e38d90ad51326d640a62eb1412b65af01bb43e5601`.
The separate threshold-cap inventory already passed every even2..40.
Spins2/4 have weak threshold limits, so the later even-analytic quadratic
cap is used by the actual integration driver instead of relying on a
single extremely short constant-margin cell.

Spin0 completed with upper
6.847159091864321571288553923423997494813397628507400578644404614008385313735,
using2950 ordinary cells and5 gaps,172.27276seconds. The full rational
partition replay also passed. Records are `left/low0/head_l0.json` and
`left/low0/replay_l0.json`. Requested budget was numerical*1.002+.0001.

Spin2 is running in `left/low2/` with requested budget
0.519411992983221 = numerical*1.01+.0001. A requested budget is not a bound;
only a completed interval record and its replay may be selected.
All other left spin components remain to be integrated. A full left
coupling bound additionally needs a same-hash continuum spin42..infinity
certificate; the g2 tail certificate must not be reused for this witness.

Numerical tail diagnostics, useful only to plan the proof:

* `audit180/source.json`: even42..100, mu4.02..49, maximum sampledV/p
  -51.85771613; numerical support1.22161836e-5.
* `audit180/uv.json`: mu>=49, maximum sampledV/p -0.0011634058;
  numerical support0.00044207208.

The compact margin appears generous. A uniform UV margin0.01 cannot be
assumed for this exact witness; even the sampled minimum negativity is
weaker. The actual continuum UV proof and integrable real-tail bound
remain separate tasks. These diagnostics do not certify omitted spins.

Reproduction uses the same scripts as the robust g2 head, with the left
witness, left moment directory, and numerical budgets file
`output/full_direction_20261007/left/head180_n160/angle180_final_audit.json`.
`head_centered_batch.py` handles spins0/2; every chosen component must
be replayed with `replay_gap.py` using its own witness and moment directory.
The full method and inherited mathematical assumptions are documented in
`HEAD220_STATUS_AND_REPRODUCTION.md` and
`ROBUST_STAGE1_HEAD_CERTIFICATION.md`.

Spin2 has now completed: upper
0.5194092127896018376002879640947783170068498303614224852552575891827281933516,
419.569seconds, requested budget passed. Its dedicated full-partition
replay is running. The rigorously integrated0+2 sum is7.3665683046539234
(decimal display); this is only two finite-spin components.

Continuation now covers every remaining even spin4..40 in separate
`left/spin{ell}/` folders. Spins4/6 use `head_quartic_batch.py`; spins8..40
use `head_gap_batch.py`. Each target is numerical*1.01+.0001 with at most
8000cells. The exact launch commands are retained in
`left/launches_4_40.json`. Every completed result is replayed separately.

All six left source tables were rebuilt in `moments128_replay/` and match
bit for bit. `assemble_left180.py` checks both copies, the polynomial
Gauss exactness degrees, source/witness/implementation hashes, completed
per-spin certificates, and their full partition replays. Its authoritative
output is `left/head_summary.json`; missing spins and missing replays are
listed explicitly. The planned complete-head ceiling is7.75783, but that
ceiling is emitted only after all21 adaptive records finish and their
exact directed dyadic sum is verified below it. No absent spin is set to0.

Spin2's full partition replay has also passed (4836ordinary cells).
Therefore the completed0+2 contributions both have replayed certificates;
remaining4..40 and the infinite-spin tail are still separate work.

After the adaptive jobs were already running, the parent independently
reviewed the exact source-polynomial contraction in`finite_fast.py`.
`check_fast_intervals.py` then passed252 comparisons of interval fields,
second-order direct jets, and fourth-order normalized midpoint/range jets
against`FiniteGap` for spins14/40 at seven compact and UV cells each. PV
mode selection was identical and fixed by the complete cell. The diagnostic
is retained in`left/fast_interval_diagnostic.json`; it is not a substitute
for the algebra review or a complete integral bound.

The new separately hashed`replay_left_family_fast.py` uses`FiniteFast` only
for not-yet-started higher-spin replays. It leaves adaptive results and all
old replay files unchanged. It recomputes the complete saved partition,
checks its total against the original requested ceiling, and records its
additional fast-source implementation hash. A successful replay therefore
provides a second valid enclosure, regardless of whether its numerical
value coincides with the old enclosure. The family assembly uses the exact
directed bound of the selected replay variant and retains its provenance.

The first fast replays for spins12/14/16 exceeded their originally allocated
quality budgets by approximately9.5e-6,7.1e-6 and2.6e-6, respectively. The
driver rejected those attempts before emitting any selected certificate.
Failure bookkeeping is retained in`left/fast_replay12_failure.json` and
`left/fast_replay14_16_failures.json`. These were enclosure-width changes
under algebraic rearrangement, not failures of domain coverage or source
error bounds. Spin12 is replayed with the original evaluator;14/16 are
rerun under slightly looser explicit ceilings. The original adaptive
records are unchanged. The final head and parent assembly retain the
larger of original and replay upper values per spin and round that actual
sum upward; the old quality budget is not confused with validity.

All21 adaptive components have now completed. Their original exact dyadic
sum displays7.75782173898639. This is a rigorous finite-head result, with
some independent partition replays still pending when this note was added.
The updated`head_summary.json` is authoritative for replay completion.

A separate refined spin0 record in`left/low0_refined/` completed with upper
6.836858149620106415961637841227750496472070209888091649960040267585162605673,
using4326 ordinary cells,320.46446seconds. Its complete partition replay
passed in130.44312seconds. This improves the spin0 contribution by
0.0103009422442151553 relative to the original completed head. It does not
replace the base summary. Once the base replays finish,
`assemble_left180_refined.py` produces a separate same-schema
`left/head_refined_summary.json`, selecting only this new0 component and
retaining every other frozen record. The infinite-spin tail remains the
same, because the underlying witness is unchanged.
