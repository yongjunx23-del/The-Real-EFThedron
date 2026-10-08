# Direct upper95: final tail reproduction record

All tail components are complete for the unchanged witness `output/six_targets_20261007/upper95_direct/witness.npz`, SHA256 `767b9d4aa9247f505004c392031f60ed82a3aedc7b1354e7a21ee9657e7a24ec`. This record concerns the omitted-spin tail; the separate complete directional certificate also includes the finite head.

The finite tail covers every even spin from 42 through 448 and every energy above threshold. Its selected upper contribution is **0.00006704866029397303**. The authoritative exact fraction, all real-norm cells, their retained paired upper bounds, the local coercivity margins, and outward reciprocal weights are in `output/six_targets_20261007/upper95_direct/real/finite_banded_cost.json`, SHA256 `5b7da07048965378927bcc2f780a5c1899240a000664cd6b06acdfbb7420c6e5`.

The independently reconstructed ledger is `output/six_targets_20261007/upper95_direct/tail/independent_coverage_accounting.json`, SHA256 `c30e8dcd6d1831349e1ad07def24d553dd077b3039433ec41196b46a4617e697`. It checks all 100,145 paired finite UV cells, all 4,096 real-norm cells, all 204 finite tail spins, every touching closed energy interval, and the final exact weighted sum. It uses a separately written interval-index reconstruction, while sharing the original physical kernel implementation.

The UV packet is `output/six_targets_20261007/upper95_direct/tail/uv_complete.json`, SHA256 `290cc0ea33105ace41d3b3923fa14483ba43ae97783282ccfdf0166a3aba21cd`. Finite spins use 400-bit interval producers and full 480-bit replay. The 16 adjacent eta panels cover every even spin at least 450 and its infinite-spin limit for the bounded tau region. The paired signed source-moment/Jensen estimate closes the unbounded tau region. The exact source algebra and four-panel crossing certificate are recorded separately and bound by this packet.

The compact proof, owned by the separate head/compact calculation, covers `(4,49]` with 17 adjacent rectangles and the analytic/finite-spin complements recorded in `compact/complete_compact_coercivity.json`. Its finite producers use 512 bits, replays 640 bits, and analytic estimates 400/512 bits. This packet supplies the common compact coercivity value 1 used in the local-cost accounting.

The infinite-spin real-norm remainder is separate. The compact squared-source norm is bounded by approximately `1.928685695143395e-19`, and the UV squared-source norm by `4.064894403944952e-10`. Divide each by twice its corresponding coercivity, respectively 1 and 1/200. The resulting infinite-spin support contribution is approximately `4.064894403954596e-8`. The complete assembler uses the exact fractions and the conservative maximum across 256/384-bit real-norm calculations.

## Selected reconstruction commands

The source and raw producer commands are in `UV_UPPER95_PROOF_AND_REVIEW.md` and the individual queue records. After those raw components and the compact packet exist, the selected acceptance commands are:

```text
python research/six_targets_20261007/assemble_uv_upper95_v2.py --angle upper95
python research/six_targets_20261007/assemble_uv_upper95_v2.py --angle upper95 --verify-only
python research/six_targets_20261007/banded_upper95_v2.py --angle upper95
python research/six_targets_20261007/review_upper95_tail_v2.py
```

Run producers in a fresh reconstruction workspace; do not overwrite frozen packets. The final recorded continuation is `run_tail95_completion_v3.py`, which only verifies the already completed UV packet and constructs the missing cost/ledger. Its command record `tail_completion_commands_v3.jsonl` and summary `tail_completion_campaign_v3.json` contain three zero exit codes. This final verification/accounting phase took 52.00 seconds; that is not the runtime of the preceding proof producers or the full directional certificate.

Two earlier metadata readers stopped before writing their respective result because an authenticated source-review JSON has a list root. The selected readers authenticate that file in full and only inspect `sha256_manifest` on dictionary roots. The repairs do not change any physical coefficient, formula, numerical proof cell, domain assertion or support bound. Original sources, failure logs, exact substitution ledgers, and the independent reviews `REVIEW_UPPER95_UV_METADATA_FIX.md` and `REVIEW_UPPER95_BANDED_METADATA_FIX.md` are preserved. The earlier UV note's v2 completion driver is superseded operationally by v3; its mathematical proof description is unchanged.

This is internal interval certification with exact arithmetic in the coverage/accounting layer, analytic remainder estimates, and higher-precision replay. It is not a formal proof-assistant result or an independent implementation of the physical kernels.
