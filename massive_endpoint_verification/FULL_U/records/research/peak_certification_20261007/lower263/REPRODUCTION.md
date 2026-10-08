# Lower-boundary peak: exact witness and reproduction

The target is near 64.44% of normalized lower-branch length, including the user's marked 65% station. Its physical supporting normal is approximately 263.27473709797 degrees. The horizontal coordinate of the difference plot is boundary length, not that angle. Each physical axis is divided by the observed dual range before measuring length.

Witness: `output/weak_regions_20261007/lower/search/angle263p2747_from261.npz`.

SHA-256: `8931658c2a353a2bbf9b9cb074db490ba5d15f044c80837236d3fd9dd801c541`.

Exact objective:

```text
a = -2109641626282615/18014398509481984
b = -8945221711437375/9007199254740992
```

The 220 multipliers are the exact binary64 values in that NPZ. Their exact fractions, row names, search settings and original numerical audits are already exported in `output/weak_regions_20261007/lower/LOWER_REFINEMENT_RECHECK.md`. This campaign introduces no new physical/null-equation row. It certifies the already audited numerical witness using the existing full 220-row catalogue in `ALL_CONSTRAINTS_HUMAN_READABLE_2026-10-07.md`.

## Finite head

The ordinary mathematical producers remain the frozen `research/boundary_length_20261007/head/certify_head_constant.py` and its imported reviewed models. Spin 0 uses Centered, spin 2 uses ConstantCentered, and even spins 4 through 40 use Fast. The Centered normalized ordinary-cell shortcut remains restricted below beta=1/2. No expanded normalized branch from the previous corrected wrappers is used.

Source moments are built twice with lmax=40, N=128, Nuv=32 and 256-bit arithmetic. Every moment and UV moment file must agree byte for byte between the two builds. All 21 spin caps are preflighted at 256 and 384 bits. The ordinary producers run at 256 bits and every complete partition is replayed at 384 bits. Analytic threshold, source-endpoint and UV gaps remain included. The head assembler retains the larger producer/replay integral for every spin.

The original quality targets were numerical spin values plus 1e-4 for spin 0, 3e-5 for spin 2, 1e-5 for spin 4 and 1e-6 for remaining spins, with maximum 7000 ordinary cells. Those targets were overly expensive. Original completed and active jobs were preserved. Only unstarted even spins 12..40 received revised targets of numerical value plus 1e-5 and maximum 3000 cells. The scheduler then allowed six workers. A missed target is a quality issue, not a failed proof: the actual completed upper and exact complete partition are used, never the requested target.

Records are under `output/peak_certification_20261007/lower263/head/`:

- `preparation.json`, original `budgets.json` and duplicated `moments128*` tables preserve original input preparation.
- `orchestration_switch.json` records the deliberate scheduler replacement; its active spin-6/8/10 workers were verified alive and adopted without restarting or changing their settings.
- `budgets_relaxed.json` records every changed target and cell budget.
- Each `spin*/commands.json` stores the final producer/replay command pair. Adopted workers retain their observed original command line and actual operating-system exit code.
- `jobs.json`, `campaign.json` and `head_summary.json` are written only after all 21 full partitions and replays finish successfully.

To reproduce cleanly, first execute the two source-build command arrays recorded in `preparation.json`, then the per-spin command arrays in the final `jobs.json`, and finally the original head assembler. There is no mathematical need to reproduce process IDs or interrupt a scheduler. Replay the saved exact partitions with the original frozen producer's `--replay` option when independently checking an existing packet. Use new output paths because proof drivers refuse to overwrite frozen records.

## Compact energy, all omitted spins

Adapters are in `research/peak_certification_20261007/compactlower263/`. The full literal adaptation ledger is `lower263/compact/adapter_manifest.json`; the parent independently reconstructed and hash-checked all nine adapters.

Run from the workspace root, setting OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and MKL_NUM_THREADS to 1:

```text
python research/peak_certification_20261007/compactlower263/analytic_anglelower263.py --bits 400
python research/peak_certification_20261007/compactlower263/analytic_anglelower263.py --bits 512
python research/peak_certification_20261007/compactlower263/analytic_anglelower263_low.py --bits 400
python research/peak_certification_20261007/compactlower263/analytic_anglelower263_low.py --bits 512
python research/peak_certification_20261007/compactlower263/build_moments.py --angle lower263 --kind roy
python research/peak_certification_20261007/compactlower263/build_moments.py --angle lower263 --kind strips
```

These require the prepared source-cache packet from `tail_lower263/anglelower263_campaign.py`. The finite compact source tables use jmax=skmax=68 and kmax=256 at 512 bits, sufficient for every finite complement required by the analytic plan. The independent moment audit checks 207 exact factorial/beta T moments and 49,512 entries shared with separately built head tables.

Eight finite complements need fresh evaluation and replay:

| Energy interval | Finite even spins | Analytic even spins |
|---|---|---|
| 13/2..7 | 42 | 44 and above |
| 7..15/2 | 42..46 | 48 and above |
| 15/2..8 | 42..52 | 54 and above |
| 20..24 | 42..46 | 48 and above |
| 24..30 | 42..52 | 54 and above |
| 30..36 | 42..58 | 60 and above |
| 36..42 | 42..62 | 64 and above |
| 42..49 | 42..68 | 70 and above |

Other analytic rectangles already start at spin 42. The separate low-energy proof covers 4 < mu <= 9/2. Required producer/replay commands are recorded in `compact/commands.jsonl` and the two successful extent-correction reruns in `compact/replay_extent_commands.json`. The producer uses 512 bits; replay uses 640 bits.

The first replay wrapper inherited a maximum-spin-60 operational guard. It safely rejected the last two finite complements before writing results. The unchanged lower-spin replays remain valid. The separate `replay_piecewiselower263_extended.py` changes only that guard to 68, matching the validated moment degrees; both last complements were successfully replayed with it. `replay_extent_adaptation.json` and the parent review preserve this correction. For a clean reconstruction, use the extended replay wrapper for the last two rows of the table.

```text
python research/boundary_length_20261007/audit_compact_moments.py --witness output/weak_regions_20261007/lower/search/angle263p2747_from261.npz --compact-moments output/peak_certification_20261007/lower263/compact/moments --head-moments output/peak_certification_20261007/lower263/head/moments128
python research/peak_certification_20261007/compactlower263/assemble_anglelower263.py
python research/peak_certification_20261007/compactlower263/derive_stronger_margin.py
```

The final compact component proves V/p < -3 on every even spin at least 42 and every 4 < mu <= 49. This is checked directly against every paired analytic rectangle and all 4294 producer/replay finite spin cells. The weakest saved margin is 3.2203079282759752. Both the original kappa=1 component and the derived stronger component are preserved.

## UV and real-norm costs

The tail agent's recorded campaigns construct same-witness exact source contractions, large-tau estimates, finite even spins 42..448, and the continuous inverse-spin domain for all spins 450 and above. Paired finite/uniform UV proofs run at 400/480 bits; large-tau proofs run at 320/400 bits. A conservative common UV kappa of 1/50 is used for the infinite domain, while finite costs retain actual local paired margins. The large-tau region has kappa=1.

`real/finite_banded_cost.json` bounds the entire finite 42..448 all-energy contribution by 0.000533703297885395. Each of its 4096 real cells uses the minimum valid margin across all 204 finite spins and every touching closed UV cell, as well as any required compact or large-tau piece. The reciprocal is rounded upward to 40 significant binary bits; all subsequent weighted sums are exact. The separate peer review independently reconstructed these overlaps and checked all 19,015 paired UV cells.

The infinite-spin real cost uses compact kappa=3 and UV kappa=1/50, with the smaller margin on the rational-cutoff sliver. Its total contribution is approximately 0.00000129020729205. These tail costs are components, not complete directional bounds on their own.

## Complete assembly and interpretation

```text
python research/peak_certification_20261007/assemble_lower263.py --out output/peak_certification_20261007/complete_lower263_certificate.json
```

Only successful assembly establishes the completed bound. The assembler checks the same witness everywhere, full head/source partitions and replays, exact compact spin/energy coverage, full UV coverage, every banded real-cost choice, the infinite spin/energy domains, the original and relaxed quality provenance, and the independent component reviews. The exact unnormalized rational inequality is authoritative; its displayed unit-normal value uses a directed lower norm bound and outward decimal rounding.

This remains an internal computer-assisted certificate under the stated physical assumptions. It is not proof-assistant verification or a proof of optimality. Comparisons with the published dual and primal are against digitized figures, without author-supplied exact witnesses or certified digitization error. The baseline catalogue and previous reproduction archive remain required dependencies.

## Completed result

The full assembly succeeded. The authoritative packet is `output/peak_certification_20261007/complete_lower263_certificate.json`, SHA-256 `47860e044d31ab8bf5454c3fe916097bdcd1068a56bfd55dfd9c4e1816fc1dc4`. Its outward unit-normal support upper bound is **0.193047976** at the exact objective direction recorded above (nominal angle 263.27473709797 degrees).

| Component | Retained upper value, approximate |
|---|---:|
| All energies, even spins 0 through 40 | 0.19251298179262752 |
| All energies, even spins 42 through 448 | 0.000533703297885395 |
| Compact energy, even spins 450 and above | 2.3125763711942743e-18 |
| UV and cutoff sliver, even spins 450 and above | 0.0000012902072920523135 |

All 21 finite-head producer/replay pairs completed. Saved producers at spins 2, 6, 8, and 10 reached their maximum cell budget before their requested quality targets; their complete interval enclosures and full higher-precision replays remain valid. The authoritative support uses the conservative maximum of each producer/replay pair, not a requested target or a quadrature estimate.

The separate final packet audit verified all 255 dependency hashes, exact component sums, the complete finite-head spin list, same-witness identities, exact tail accounting, and conservative direction normalization. Its record is `output/peak_certification_20261007/lower263/final_packet_review.json`; reproduce it with `python research/peak_certification_20261007/lower263/audit_complete_packet.py` in a fresh reproduction directory. The audit refuses to overwrite its result, as does the complete assembler.

The digitized published dual support at this direction is approximately 0.19662451062072475, so the new certified support is about 1.82% smaller. This comparison is to the digitized curve; it does not certify the digitization or independently reconstruct the published authors' certificate.
