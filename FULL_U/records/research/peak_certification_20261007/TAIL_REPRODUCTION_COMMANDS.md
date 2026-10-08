# Commands for the two complete tail calculations

Run from the reconstructed workspace root, in a disposable output workspace. Preserve the archived evidence. Set `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and `MKL_NUM_THREADS` to `1`. Use the Python environment and library versions recorded in the inherited reproduction handoff. Exact multipliers are read from the frozen witness NPZs; no optimizer needs to be rerun for certification.

For upper83:

```text
python research/peak_certification_20261007/tail_upper83/prepare_angleupper83.py
python research/peak_certification_20261007/tail_upper83/certify_angleupper83_large_tau_moments.py --precision 320
python research/peak_certification_20261007/tail_upper83/run_angleupper83_finite_queue.py
python research/peak_certification_20261007/run_uniform_tail_queue.py --case upper83 --workers 8
python research/peak_certification_20261007/run_tail_auxiliary.py --case upper83
python research/peak_certification_20261007/assemble_uv_peak.py --angle upper83
python research/peak_certification_20261007/assemble_uv_peak.py --angle upper83 --verify-only
```

For lower263:

```text
python research/peak_certification_20261007/tail_lower263/prepare_anglelower263.py
python research/peak_certification_20261007/tail_lower263/certify_anglelower263_large_tau.py --precision 320
python research/peak_certification_20261007/tail_lower263/run_anglelower263_finite_queue.py
python research/peak_certification_20261007/run_uniform_tail_queue.py --case lower263 --workers 1
python research/peak_certification_20261007/run_tail_auxiliary.py --case lower263
python research/peak_certification_20261007/assemble_uv_peak.py --angle lower263
python research/peak_certification_20261007/assemble_uv_peak.py --angle lower263 --verify-only
```

The finite queues cover the exact groups 42:90, 92:142, 144:192, 194:244, 246:294, 296:346, 348:396, and 398:448 in steps of two. Each group runs a 400-bit producer then a 480-bit replay of that producer's complete partition. The upper uniform queue uses `certify_angleupper83_uniform.py` separately for groups `i:i`, i=0,...,15; lower uses `certify_anglelower263_uniform_panels1.py` with `0:0`. Both queues replay every cell at 480 bits. Saved queue logs and JSON preserve each command and exit status where recorded; the finite queue source is the complete deterministic command specification.

The auxiliary driver records all exact commands. It runs the 400-bit large-tau calculation, two exact source regeneration/algebra checks, finite and infinite real-norm integrations with 4096 cells and 256/384-bit replays, and the independent real-field check. Its scientific outputs and producer/replay packets bind the final source hashes.

After the separately documented compact coercivity calculation is complete, run:

```text
python research/peak_certification_20261007/banded_upper83.py --angle upper83
python research/peak_certification_20261007/review_upper_tail.py
python research/peak_certification_20261007/banded_lower263.py --angle lower263
python research/peak_certification_20261007/review_lower_tail.py
```

Upper banding uses the original compact kappa1 packet; lower banding uses its independently checked stronger kappa3 packet. Both read the same-witness real-norm cells and retained local UV margins. The ledger checks independently reconstruct all touching closed intervals and outward weighted sums. Complete directional assembly also requires the separately recorded finite-head certificate and its selected threshold corrections.

The initial point screens, upper absolute-source large-tau failure, unused 16-panel lower probe, and unreplayed four-panel upper probe are preserved diagnostics. They are not required as proof components and must not replace the selected paired partitions or the sharpened upper large-tau producer. Read `UV_TAIL_PROOF_AND_REVIEW.md` and the canonical current-certificate registry before selecting results.
