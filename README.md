# Three certified massive endpoints

This directory verifies only these three inequalities, in the normalization
`M/(32*pi) = g0 + g2*(sbar²+tbar²+ubar²) + ...`, with `m=1`:

| Certificate | Exact objective | Outward decimal result |
|---|---|---|
| C08 | `g0` | `g0 <= 2.705955` |
| C20 | `g2` | `g2 <= 1.052351` |
| C21 | `-g0` | `g0 >= -7.751502` |

The underlying support constants are stored as exact rational numbers in
`FULL_U/data/selection_registry.json`. These three supports do not determine
an entire curved boundary. They are not a new optimization or an optimality claim.

## Run the checks

Use Python 3.12, with the tested packages in `requirements-verify.txt`.
The original interval runtime was python-flint 0.9.0. Install a wheel appropriate
for your operating system; platform-specific binaries are not bundled.
Do not run Python with `-O`, because the preserved proof scripts use assertions.

After extracting the archive, open a terminal in `massive_endpoint_verification/verification`.
Create the environment outside the extracted package:

```sh
python3.12 -m venv ../../verification-venv
../../verification-venv/bin/python -m pip install -r requirements-verify.txt
../../verification-venv/bin/python -B verify.py --out ../../endpoint-check
```

On Windows, use `../../verification-venv/Scripts/python.exe` in place of
`../../verification-venv/bin/python`.

Choose a new output directory outside `verification`. The driver makes a
working copy, runs the checkers there, keeps reports in the output directory,
and removes that temporary copy. The frozen inputs are not modified.
Allow roughly twice the unpacked verification size in free disk space.
No optimizer or solver is required.

The default check does the following:

1. Reconstructs all 590 exact multipliers from the three frozen binary witnesses:
   150 for C08 and 220 each for C20 and C21.
2. Checks the enclosing release manifest, original dependency hashes and exact saved head partitions, including
   the finite-width endpoint caps for the supported head-record schemas.
3. For C08, checks both refined low-spin records and their replays, retaining the
   larger upper value for each spin 0, 2, 4, 6. It reruns the original directed
   merge of the spin 8–40 real-norm masses with their cellwise coercivity bounds,
   with the objective restricted to slope zero. It checks and assembles the
   original all-energy, spin >=42 tail using its `31/50` coercivity margin.
4. For C20, reruns the original compact-energy, base UV and strengthened UV
   sign assemblers, checking exact partitions, individual negative enclosures,
   and the analytic domains covering infinite spin and energy.
5. For C21, invokes the original `Coercivity` checker: every even spin 42–448,
   exact adjacent UV partitions, paired higher-precision records, and strictly
   negative saved upper enclosures. Its compact-energy and infinite-spin
   components are authenticated and checked by the original assembly route.
6. Runs the original C20 and C21 endpoint assemblers and requires exact rational
   equality with all three selected support constants. C20 uses both
   `--strong-uv` and `--refined-head`. C21 explicitly uses `local_cost2048.json`;
   its historical default `banded_cost.json` gives a different, weaker result.

These are substantive checks of the saved interval certificates, their
partitions and assembly. They are not fresh evaluations of every interval
integral or a machine-checked proof in a theorem prover. Hashes authenticate
which inputs were used; hashes alone are not the mathematical justification.

## A fresh interval calculation

To additionally recompute the C21 finite-spin real-norm contribution at 256 and
384 bits, run:

```sh
../../verification-venv/bin/python -B verify.py --out ../../endpoint-check-fresh --fresh-local
```

This recomputes the weighted integral for every even spin 42–448, using 2,048
energy cells and the saved coercivity proofs. It checks exact agreement with
the saved component. It does not recompute the low-spin heads, the source
moment tables, or the infinite-spin tail. A complete fresh certification
requires regenerating those records with the retained original producers.
The proofs remain conditional on the stated dispersion, crossing, analyticity
and unitarity assumptions.

## What is retained

- `FULL_U/data/multipliers_exact.json`: all 590 multipliers, with row names,
  rational values and original binary64 representation.
- `FULL_U/data/selection_registry.json`: the three exact objective directions,
  supports and selected component ledgers.
- `FULL_U/records/`: substantive source moments, interval cells, endpoint caps,
  finite-spin sign enclosures, real-norm estimates, infinite-spin/infinite-energy
  remainder records and unchanged original proof programs.
- `FULL_U/appendix.md`: the original mathematical description of the kernels and
  certification. It describes the broader original calculation; this smaller
  release targets only the three endpoints above.
- `FULL_U/DEPENDENCIES.json`: byte identities of retained original records.
- `FULL_U/bundle_run.py`: the original narrow path-relocation adapter, allowing
  Windows-authored record paths to run on other systems.
- `FULL_U/records/research/light_endpoint_replay.py`: a new adapter selecting
  C08's slope-zero finite merge and invoking C21's original sign checker.

Historical parent records and code needed by authenticated manifests remain
included. Their presence is not an assertion that all original boundary
certificates are distributed or checked here. In particular, the original
15-direction C08 producer remains as provenance, but this release runs its
endpoint projection instead of requiring the other fourteen targets.

Start the derivation with `DERIVATION_GUIDE.md`. The release manifest supplied
by the enclosing package covers this directory, including the new driver and
adapters.
