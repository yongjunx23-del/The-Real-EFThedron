# Three additional near-free-point certificates

The separate ledger is
`output/full_direction_20261007/certification/complete_free_point_certificate.json`.
The frozen fifteen-direction ledger and its packet manifest are unchanged.

The complete, outward-rounded inequalities are:

| Exact objective | Certified upper |
|---|---:|
| g0 - 32 g2 | 0.960024 |
| g0 - 64 g2 | 0.427128 |
| g0 - 128 g2 | 0.211575 |

All use the original fixed150-row witness, with the objective changed to
`(1,b)` and **no X4 or other multiplier change**. Thus the exact effective
X4 increment is zero in every case. There is no positive interpolation in
these three reported constants. The ledger records both exact rational
constants and approximate unit-normal plotting values. Their angles lie
just above270degrees; the angle labels are not the exact mathematical
definition of the normals.

The proof consists of the same complete decomposition as the previous
certificates:

1. Fresh interval integration of spins0,2,4,6 at280bits and340bits, including
   all source endpoint gaps and threshold caps. Exact partition coverage
   and totals are checked, and the larger of the two complete upper sums
   is used for each spin. These are same-implementation replays, not an
   independent implementation.
2. Exact overlap merge of the existing9790-leaf dense spectral partition
   with19790 real-kernel mass cells, separately for every even spin8..40.
   The negative shift only strengthens each saved V/p inequality.
3. All even spins42 and above, at every energy, use unchanged A40 and
   `V_new/p <= -(31/50-b)`. Consequently the remainder is bounded by
   `A40/[2(31/50-b)]`. This is a direct consequence of the already checked
   full-domain base-tail certificate; no new sampling or UV cutoff enters.

The old `b=-16` bound together with g2>=0 also bounds these three directions
by1.766916. Every new direct inequality is strictly stronger than that
consequence. This comparison does not assert superiority to any published
bound; the published comparison is a separate task.

## Reproduce

Use the same Python environment and frozen dependencies as the other2026-10-07
certificates. For each spin l in0,2,4,6, run:

```
python research/full_direction_20261007/head_family_extended.py --spin l --slopes=-32,-64,-128 --precision 280 --out output/full_direction_20261007/certification/free_point_head_p280_spinl.json
python research/full_direction_20261007/head_family_extended.py --spin l --slopes=-32,-64,-128 --precision 340 --out output/full_direction_20261007/certification/free_point_head_p340_spinl.json
```

Replace the literal l in the command and filename with the selected spin.
Then run:

```
python research/full_direction_20261007/combine_free_point_finite.py
python research/full_direction_20261007/assemble_free_point.py
```

Observed per-spin head runtimes were18–60seconds; the finite merge took
11seconds and final exact assembly about4seconds. Each output retains its
own timing and provenance hashes. No old source or certificate was edited.
