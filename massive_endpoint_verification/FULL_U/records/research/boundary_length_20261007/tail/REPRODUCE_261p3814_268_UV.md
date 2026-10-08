# Reproduce the 261.3814° and 268° ultraviolet components

These commands reproduce only the ultraviolet coercivity components. They do not by themselves certify the finite head, compact-energy region, real-part norms, or complete physical bounds. The final archive must include those separate same-witness components.

Run from the workspace root in a disposable extracted copy of the archive. Production commands write proof files; preserve the original frozen archive. Before a full regeneration, move the two generated campaign directories `output/boundary_length_20261007/tail/angle261p3814` and `angle268` to backup locations in that disposable copy, so the new runs start with empty output directories. Do not mix freshly rewritten preparation records with old proof records: timing-dependent preparation hashes will differ. Python must have NumPy, SymPy, and python-flint available. The scripts load the archived `tmp/certification_deps` before importing FLINT. The archived final manifests identify exact source bytes; proof timing fields need not reproduce byte for byte.

Set process-local thread limits:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
```

Use `python` below for the interpreter in the archive environment. The development interpreter was `C:\Users\laure\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`.

For read-only verification of the completed archived packets, without regenerating or rewriting anything, run:

```powershell
python research/boundary_length_20261007/tail/assemble_uv261p3814.py --angle 261p3814 --scope all --verify-only
python research/boundary_length_20261007/tail/assemble_uv268.py --angle 268 --scope all --verify-only
```

These commands check every saved paired cell's exact directed sign, full rational coverage, source/objective identity, recursive dependency hashes, and agreement with the frozen assembled records. They do not recompute the interval kernels; the full production/replay commands below do that in the disposable copy.

## 261.3814°

Frozen witness: `output/boundary_length_20261007/search/angle261p3814.npz`.
SHA256: `5bd1b99e5fd4e237fe7e98583aa683ab076d4d76319b05b06c023e583f48383e`.

```powershell
python research/boundary_length_20261007/tail/prepare_angle261p3814.py
python research/boundary_length_20261007/tail/rebuild_source_cache261p3814.py --angle 261p3814
python research/boundary_length_20261007/tail/certify_angle261p3814_large_tau.py --precision 320
python research/boundary_length_20261007/tail/certify_angle261p3814_large_tau.py --precision 400
python research/boundary_length_20261007/tail/certify_angle261p3814_uniform.py --precision 400
python research/boundary_length_20261007/tail/certify_angle261p3814_uniform.py --precision 480 --replay output/boundary_length_20261007/tail/angle261p3814/uniform_p0_3_400_k1_50.json
python research/boundary_length_20261007/tail/run_angle261p3814_finite_queue.py
python research/boundary_length_20261007/tail/assemble_uv261p3814.py --angle 261p3814 --scope all
```

The finite queue runs eight disjoint spin groups with two single-thread workers. Each group first produces a 400-bit partition, then replays every saved cell at 480 bits. The queue log records group runtimes, while per-spin producer and replay logs preserve cell counts and runtimes. A saved successful file causes the queue to skip that stage, supporting continuation of an interrupted regeneration with unchanged preparation bytes.

Claimed UV margin: `1/50`. Infinite spins use four exact eta panels.

## 268°

Frozen witness: `output/boundary_length_20261007/search/angle268.npz`.
SHA256: `1a0754b19846867214d914c595cdc744f433abfcdaff309d645677bbe2602044`.

```powershell
python research/boundary_length_20261007/tail/prepare_angle268.py
python research/boundary_length_20261007/tail/rebuild_source_cache268.py --angle 268
python research/boundary_length_20261007/tail/certify_angle268_large_tau.py --precision 320
python research/boundary_length_20261007/tail/certify_angle268_large_tau.py --precision 400
python research/boundary_length_20261007/tail/certify_angle268_uniform.py --precision 400
python research/boundary_length_20261007/tail/certify_angle268_uniform.py --precision 480 --replay output/boundary_length_20261007/tail/angle268/uniform_p0_0_400_k1_500.json
python research/boundary_length_20261007/tail/run_angle268_finite_queue.py
python research/boundary_length_20261007/tail/assemble_uv268.py --angle 268 --scope all
```

Claimed UV margin: `1/500`. A single exact eta panel suffices. The sampled joint limit approaches approximately `-0.00289268` at tau zero, so this fixed witness cannot have the larger global margin `1/100`.

## Scope, seams, and checks

Both campaigns cover every even spin at least 42 and every energy mu at least 49. Finite spins 42 through 448 cover tau from zero through `min(256,l(l+1)/45)`. Infinite spins at least 450 cover tau zero through 256 and all eta from zero through `1/(450*451)`. The separate all-negative shifted-crossing-power Jensen estimate covers tau at least 256 for every spin at least 108, with margin 1. Spins through 106 need no large-tau piece because their entire mu-at-least-49 domain has tau below 256.

Acceptance requires a finite directed interval upper strictly below the claimed negative margin. The assembler verifies exact rational closed partitions, both precision runs, exact objective components, omitted-crossing signs, all source hashes, and the unbounded-domain seams. It recursively binds source caches and dependency files into the final `uv_complete.json` manifest. The rebuild verifies all 972 exact source moments and all nine remainder coefficients. Higher-precision replays use the same reviewed formula; they are not advertised as independent mathematical derivations.

The screening scripts are optional diagnostics and are not inputs to a valid continuum certificate. The complete archive should retain the exact builder, all cached data, producer/replay JSON files, logs, notes, source reviews, final manifests, and root's separate exact algebra audits.
