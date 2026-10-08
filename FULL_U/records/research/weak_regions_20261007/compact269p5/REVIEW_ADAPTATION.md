# Compact269.5 adaptation review

This adapts the frozen compact268 proof to the exact new269.5 witness; physical inequalities and mathematical evaluators are unchanged.

Witness: `output/weak_regions_20261007/lower/search/angle269p5.npz`.
SHA256: `f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae`.
Exact objective: `(-2515252609197481/288230376151711744, -9006856288192973/9007199254740992)`.

The changed items are import bootstrapping, witness SHA/objective, source-cache locations, new output paths, adapter module names, and the adaptation-review binding. The original shared analytic85 helper is imported from its original location and bound by its actual file hash. Its old main entrypoint is not executed. The all-power Jensen/Bernstein, positive-Legendre denominator-ratio, low-energy derivative, finite polynomial/PV/geometric-remainder, interval partition, and exact cover logic are unchanged.

Analytic400/512 runs select identical17 rectangles and finite spin complements: even42..46 on[24,30],42..50 on[30,36],42..54 on[36,42],42..60 on[42,49]. All other analytic rectangles begin at42. The separate low proof covers4<mu<=9/2. Fresh compact moments therefore use jmax=skmax=60 and kmax=256 at512bits; these choices are sufficient for every finite complement. Each accepted finite interval is replayed at640bits.

An independent collaborating agent compared analytic269p5, compact_finite269p5, finite_slice, replay_piecewise and assemble to their268 counterparts and confirmed unchanged mathematics and exact objective consistency. This is internal source review using shared mathematical implementations, not proof-assistant verification or a new independent mathematical implementation.

## Independently reviewed file hashes

- `analytic_angle269p5.py`: `acbd9b8ae4f8113254e1a3368b745a8dd988a57d45751a13881ee3d7805d13b5`
- `analytic_angle269p5_low.py`: `fe7a751484e6ae34c928c8c00e90248636b8573e98ecbd367e2b3bafe359db20`
- `compact_finite269p5_combined.py`: `0b8aaab34e3be45eea1b031ce72bec0df4becfde7e3dee6f6f316df87ce16636`
- `certify_finite269p5_slice.py`: `bb36bd9748d372c9334facccd1f37e9cd1fc27dbc1422f1009d0270a4ab9151d`
- `replay_piecewise269p5.py`: `d87615800d4581fff3229ea44083d3e24dca48e0df7ae988c1ed9770ad4f4963`
- `assemble_angle269p5.py`: `7ee31ba929ba47159833d67fcf85ac1645839eb59a6f28b5dfdf2e0dd4c973b1`
- `run_piecewise269p5.py`: `df06c4a7a92f5ceb791dcd8b8d5a66cdf4c82978ba9316b626554d312754a7f8`
- `validate_compact269p5_combined.py`: `c875c8ffec67c0cc83182bd055d86e68b546e3ce4b47ba119c59acb6e45b0625`

The five main mathematical adapters were compared line by line against the frozen 268-degree versions. Differences are restricted to the documented imports, exact witness constants, paths, names and review bindings. The unchanged source geometry is Roy B=119/10 with Jmax=6 and strip T=4/3. The new witness SHA and exact objective were checked directly against the saved NPZ.
