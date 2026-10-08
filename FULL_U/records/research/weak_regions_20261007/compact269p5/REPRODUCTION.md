# Reproduce the fresh269.5 compact-energy proof

This is a coercivity component, not a complete directional bound. It proves `V/p < -50` for every even spin `ell >= 42` and every energy `4 < mu <= 49`, for the frozen269.5-degree witness. The separately constructed head, UV and real-norm components must be assembled before claiming a continuum support bound.

Witness: `output/weak_regions_20261007/lower/search/angle269p5.npz`.
SHA256: `f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae`.
Exact objective: `(-2515252609197481/288230376151711744, -9006856288192973/9007199254740992)`.

Use the pinned Python and python-flint environment in the main reproduction archive. Run from the workspace root, with one thread per process:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
```

The new adapter imports `research/weak_regions_20261007/tail269p5/angle269p5_campaign.py` and validates its prepared source manifest. Prepare that tail/source cache first using the tail campaign instructions. Rebuild into an empty disposable copy; do not overwrite these archived completed proof packets. The moment builder intentionally refuses to overwrite completed files.

```text
PYTHON research/weak_regions_20261007/compact269p5/analytic_angle269p5.py --bits 400
PYTHON research/weak_regions_20261007/compact269p5/analytic_angle269p5.py --bits 512
PYTHON research/weak_regions_20261007/compact269p5/analytic_angle269p5_low.py --bits 400
PYTHON research/weak_regions_20261007/compact269p5/analytic_angle269p5_low.py --bits 512
PYTHON research/weak_regions_20261007/compact269p5/build_moments.py --angle 269p5 --kind roy
PYTHON research/weak_regions_20261007/compact269p5/build_moments.py --angle 269p5 --kind strips
PYTHON research/weak_regions_20261007/compact269p5/run_piecewise269p5.py
PYTHON research/boundary_length_20261007/audit_compact_moments.py --witness output/weak_regions_20261007/lower/search/angle269p5.npz --compact-moments output/weak_regions_20261007/fresh269p5/compact/moments --head-moments output/weak_regions_20261007/lower/head_certification/angle269p5/moments128
PYTHON research/weak_regions_20261007/compact269p5/assemble_angle269p5.py
PYTHON research/weak_regions_20261007/compact269p5/derive_stronger_margin.py
```

The two moment builders may run concurrently. The piecewise driver waits for their completed manifests, checks the numerical kernel against the interval implementation at spins42 and60, and runs at most two producer/replay workers concurrently. `commands.jsonl` stores every finite producer/replay command and elapsed time. The independent moment audit additionally requires the fresh same-witness head moment tables.

The analytic inequalities are the reviewed all-powers Jensen/Bernstein compact envelopes and monotone positive-Legendre denominator bounds. Paired400/512-bit runs choose the same17 rectangles. A separate low proof covers `4 < mu <= 9/2`. Only four analytic rectangles need finite complements:

| Closed energy interval | Finite spins proved/replayed | Analytic even spins |
|---|---|---|
|24..30|42,44,46|48 and above|
|30..36|42,44,46,48,50|52 and above|
|36..42|42,44,46,48,50,52,54|56 and above|
|42..49|42,44,46,48,50,52,54,56,58,60|62 and above|

All other energy rectangles have an analytic proof already from spin42. The exact interval union and finite/analytic spin coverage are checked during assembly.

Fresh compact source moments use `jmax=60`, `skmax=60`, `kmax=256`,512-bit arithmetic and validated degree-exact Gauss integration. Their independent audit checks all183 T moments against exact factorial/beta formulas, and43,328 shared entries against independently built head tables. The finite proof uses exact moments, removable-pole/PV formulas and explicit geometric-remainder bounds. Every accepted512-bit cell is recomputed at640bits. There are828 producer spin cells and828 replay spin cells.

The first completed packet retains the original target `kappa=1`:

`output/weak_regions_20261007/fresh269p5/compact/complete_compact_coercivity.json`

It is unchanged. The separate stronger packet is:

`output/weak_regions_20261007/fresh269p5/compact/complete_compact_coercivity_stronger.json`

SHA256: `d05a251ad4efb576a3811dce27bbc1a22ac89776897d100c6c262129f303f7e9`.

The stronger derivation reads and hash-checks every paired analytic rectangle and every paired finite cell. Every retained upper endpoint is below-50. The weakest analytic margin is53.00566421925409 on[36,42], and the weakest finite margin is197.75009530501174. Consequently the conservative common margin50 follows directly from the saved enclosures. No new integration estimate is introduced. `stronger_margin_ledger.json` records exact rational margins and the full coverage checks; the stronger packet binds62 dependency files.

Observed builder durations were35.782s for Roy moments and88.003s for both strip tables, executed concurrently. Complete finite producer/replay timings are in `commands.jsonl`; their durations must not be summed to claim elapsed wall time for concurrent work.

The adaptation review is `REVIEW_ADAPTATION.md`, including a second agent's source-diff review and all eight adapter hashes. It confirms the old mathematical evaluators are unchanged. This is internal review and directed-interval replay using shared mathematical implementations; it is not proof-assistant verification or an independent external derivation.
