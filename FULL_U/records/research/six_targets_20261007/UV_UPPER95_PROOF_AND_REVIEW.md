# Direct upper95 same-witness ultraviolet proof

The fixed witness is `output/six_targets_20261007/upper95_direct/witness.npz`, SHA256 `767b9d4aa9247f505004c392031f60ed82a3aedc7b1354e7a21ee9657e7a24ec`, with exact objective `(1203580983497459/9007199254740992,4463211602388183/4503599627370496)`. All coefficients are unchanged from its saved numerical optimization; no rounding repair is introduced here.

The kernel, source expansion, remainder estimates, spin cutoffs, energy cutoffs, eta partition and precision replays are identical to the reviewed upper86 setup described in `UV_TAIL_PROOF_AND_REVIEW.md`, evaluated freshly for this witness. The finite UV covers even42..448, the uniform eta proof all even≥450, and the large-tau Jensen estimate closes tau≥256 for all even≥108. Exact overlap106*107/45<256 closes the domain without an unchecked spin or energy interval.

One sufficient crossing-coefficient test is sharpened. The shifted polynomial at k=y+10 is bounded on four exact adjacent z intervals `[3i/32,3(i+1)/32]`, i=0..3. For each shifted power, the maximum of all Bernstein coefficients over these panels is a uniform upper coefficient over0≤z≤3/8. Every maximum is nonpositive, so the same negative-power Jensen argument applies. The single-panel sufficient test was inconclusive for the fourth shifted power; it is not a counterexample to the inequality. The exact original k64 omitted-tail sign test succeeds on one panel and remains unchanged.

The standalone ledger `upper95_direct/tail/fourpanel_crossing_review.json` preserves every panel endpoint, every Bernstein coefficient, the exact global maximum for each shifted power, the unchanged k64 test and equality with the large-tau coefficient packet. `audit_upper95_fourpanel.py` regenerates this ledger exactly. In `angleupper95_models.py`, the inherited constructor precondition is replaced by the existing exact four-panel certificate; the numerator and all field evaluations remain unchanged. `REVIEW_SIGNED_SOURCE2_UPPER95.md` independently reviews this adaptation and the signed source2 convention/remainder.

The signed source2 estimate retains−m0, bounds higher moments by their absolute coefficients at zmax=3/143 and adds the geometric remainder with ratio(317/60)zmax<1. The equalitym0=Hs is asserted. No source term is discarded. All source moments, source norms, polynomial source contractions, paired interval cells and analytic coefficient bounds are saved and hash-bound. Exact source regeneration, polynomial specialization and independent real-field checks supplement the domain proof.

Commands, with the same single-thread environment as the upper86 reproduction note:

```text
python research/six_targets_20261007/tail_upper95/prepare_angleupper95.py
python research/six_targets_20261007/tail_upper95/certify_angleupper95_large_tau_moments.py --precision 320
python research/six_targets_20261007/audit_upper95_fourpanel.py
python research/six_targets_20261007/tail_upper95/run_angleupper95_finite_queue.py
python research/six_targets_20261007/run_uniform_tail95_queue.py --case upper95 --workers 2
python research/six_targets_20261007/run_tail95_auxiliary.py --case upper95
python research/six_targets_20261007/assemble_uv_upper95.py --angle upper95
python research/six_targets_20261007/assemble_uv_upper95.py --angle upper95 --verify-only
python research/six_targets_20261007/banded_upper95.py --angle upper95
python research/six_targets_20261007/review_upper95_tail.py
```

The finite queue uses two workers. The auxiliary queue records the400bit large-tau run, exact source/algebra checks,4096-cell256/384bit real norms and independent real-field check. Full directional assembly additionally requires the separately certified compact-energy region and finite head. A successful point screen, source sign test or UV-only packet must not be presented as a complete directional bound.

All producer/replay calculations use directed interval arithmetic and exact rational coverage, with shared implementations and internal review. This is not proof-assistant verification, external independent verification, or an optimality claim. Any comparison with published curves is limited by the digitized PDF reference.
## Selected metadata reader

The original UV assembler stopped before writing any component packet because the source-difference review JSON has a list root, while its generic metadata reader assumed a dictionary. All numerical and analytic proof records had already passed. The original source, command and failed log are preserved.

The selected replacement is `assemble_uv_upper95_v2.py`. It authenticates the entire list-root review file by SHA256 and only attempts recursive `sha256_manifest` traversal for dictionary-root JSON. This is the same root-type handling used by the established complete-packet reader. It changes no physical field, interval formula, coverage assertion, coefficient, precision or support cost. `metadata_reader_fix.json` records the exact source substitutions and old/new file hashes; `REVIEW_UPPER95_UV_METADATA_FIX.md` records the separate internal review.

The selected recorded driver is `run_tail95_completion_v2.py`. Its output command log and campaign metadata use a `_v2` suffix to preserve the failed first attempt. For manual reconstruction, use `assemble_uv_upper95_v2.py --angle upper95` and then the same command with `--verify-only`; the banded-cost and independent ledger commands remain unchanged.

