# Reproduce the repaired83 peak certificate campaign

This campaign targets the upper graphical-comparison peak near75.887% of the boundary's normalized arclength, whose outward normal is83.1759644845 degrees. Boundary arclength and normal angle are different coordinates. A full directional bound is claimed only when same-witness head, compact, ultraviolet and real-norm components have been assembled.

Witness: `output/weak_regions_20261007/upper/angle83_repaired.npz`.
SHA256: `68e39591f3f9a9d3b5cdfe3ade3924fecb227a329de4a9543729f0607e8f663a`.
Exact objective: `(535119985632953/4503599627370496,8943390007067311/9007199254740992)`.
All220 exact dyadic coefficients are preserved in `output/weak_regions_20261007/upper/angle83_repaired_exact_coefficients.json` and copied to `upper83/head/exact_witness_definition.json`. The full numerical-search constraints and repair provenance remain in the previously sealed upper numerical report. This campaign changes no coefficient or physical row.

Use the pinned Python and python-flint runtime from the reproduction archive. Run in an empty disposable reconstruction with the original relative paths; do not overwrite frozen proof records. For every process set:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
```

The head uses all even spins0..40 at every energy, with independent source table rebuilds at256bits, `N=128`, `Nuv=32`, `jmax=74`, and regular source/threshold endpoint caps. Producer integration uses256bits; every final partition is replayed at384bits. Centered models are used for spins0/2 and Fast models for4..40. The ordinary frozen Centered driver restricts its normalized branch to its valid domain. No widened Roy-source window is used.

```text
PYTHON research/peak_certification_20261007/upper83/prepare_head.py
PYTHON research/peak_certification_20261007/upper83/run_head.py --workers 3
PYTHON research/peak_certification_20261007/upper83/refine_spin2_threshold.py
PYTHON research/boundary_length_20261007/head/certify_head_constant.py --witness output/weak_regions_20261007/upper/angle83_repaired.npz --moments output/peak_certification_20261007/upper83/head/moments128 --out output/peak_certification_20261007/upper83/head/spin2_refined/replay_l2.json --replay output/peak_certification_20261007/upper83/head/spin2_refined/head_l2.json --bits 384
```

The optional but selected spin2 refinement replaces the quadratic threshold cap by the unchanged negative-constant cap. Its smaller gap is beta0..1/1024;64 new cells cover1/1024..1/256. Every original interior cell is preserved. The entire new partition, including all original cells, is independently replayed. The original and refined results remain separate; final assembly must explicitly select and bind the refined pair.

Prepare the shared tail/source packet using the `tail_upper83` instructions before compact work. The compact proof covers all even spins42 and above for4<mu≤49. Paired analytic400/512-bit envelopes choose17 energy rectangles; finite complements reach spin160. The finite normalized formula is regular atmu4, so the proof directly covers the threshold without using the failed low42 analytic shortcut.

```text
PYTHON research/peak_certification_20261007/compact_upper83/analytic_angleupper83.py --bits 400
PYTHON research/peak_certification_20261007/compact_upper83/analytic_angleupper83.py --bits 512
PYTHON research/peak_certification_20261007/compact_upper83/build_moments.py --kind roy
PYTHON research/peak_certification_20261007/compact_upper83/build_moments.py --kind strips
PYTHON research/peak_certification_20261007/compact_upper83/run_piecewiseupper83.py
PYTHON research/boundary_length_20261007/audit_compact_moments.py --witness output/weak_regions_20261007/upper/angle83_repaired.npz --compact-moments output/peak_certification_20261007/upper83/compact/moments --head-moments output/peak_certification_20261007/upper83/head/moments128
PYTHON research/peak_certification_20261007/compact_upper83/assemble_angleupper83.py
```

The literal adapter-generation script records the original relocation. The final saved adapter files additionally include the reviewed threshold-domain and spin-range changes and should be used as supplied; regeneration alone does not apply those changes. Complete dependency hashes select the final code. All commands and timings for finite compact producers/replays are recorded in `compact/commands.jsonl`. Every512-bit finite cell is independently replayed at640bits. Fresh compact source tables are degree-exact polynomial integrations with explicit truncation errors; the independent audit checks factorial/beta formulas and overlaps with the independently built head moments.

These are directed-interval computations plus analytic infinite-domain estimates, with shared implementations and internal review. They are not proof-assistant verification or an independent external derivation. Published graphical supports are comparison estimates extracted from the PDF, not exact author-supplied certificates.

## Assemble and recheck the selected complete bound

The completed UV packet is `upper83/tail/uv_complete.json`; its producer, replay and exact algebra checks are described in `research/peak_certification_20261007/UV_TAIL_PROOF_AND_REVIEW.md`, with executable command details in the neighboring `TAIL_REPRODUCTION_COMMANDS.md`. In the large-tau region the signed source2 moment estimate replaces a much looser absolute norm. The separate formula review is `REVIEW_SIGNED_SOURCE2_UPPER83.md`. Its geometric remainder remains explicit; no source term is dropped. Finite UV spins42..448 have paired400/480-bit cell enclosures, the analytic infinite-spin region begins at450, and the large-tau closure uses paired320/400-bit proofs.

First select the refined spin2 producer and its full replay. Assemble a complete bound with the global compact and UV coercivity constants. Preserve that result, then replace its finite42..448 tail component by the tighter local-margin accounting result:

```text
PYTHON research/peak_certification_20261007/upper83/assemble_selected_head.py --witness output/weak_regions_20261007/upper/angle83_repaired.npz --base output/peak_certification_20261007/upper83/head --spin2-directory spin2_refined
PYTHON research/peak_certification_20261007/assemble_upper83.py --out output/peak_certification_20261007/complete_upper83_global_certificate.json
PYTHON research/peak_certification_20261007/banded_upper83.py --angle upper83
PYTHON research/peak_certification_20261007/review_upper_tail.py
PYTHON research/peak_certification_20261007/assemble_upper83_banded.py
```

The last command writes the canonical `complete_upper83_certificate.json`; both assemblers refuse to overwrite an existing result. The banded finite cost is formed from4096 real-norm cells. On each cell its coercivity constant is the minimum across all204 even spins42..448 and every overlapping compact, closed UV or large-tau interval. Both neighboring UV cells are included at exact seams. Each reciprocal weight is rounded upward to40 significant binary bits, and products and sums are exact rationals. A separately implemented interval-index reconstruction checks all112,962 paired finite-UV cells and the4096-by204 coverage ledger; it reproduces the same exact finite cost, approximately0.00004384769205095601. This substitution changes no head or infinite-spin component.

Both complete packets retain the exact objective, exact support, exact component sums, all witness coefficients through the bound definition, and recursive dependency hashes. The displayed unit-normal support divides an outward9-decimal support by a256-bit exact rational lower bound on the direction norm and rounds upward again. Rechecking the bound is distinct from rerunning the optimizations: the saved coefficient vector is a fixed witness, and only the certificate computation is required to verify its half-plane.

## Frozen result and timing

The canonical complete packet is `output/peak_certification_20261007/complete_upper83_certificate.json`, SHA256 `d24ba5d73fb886571bc360fcdb802f6b708e2b14b12d5d5c0e6f3061d8a37cd9`. Its status is `PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND`. It recursively binds454 files. The exact rational values in the JSON are authoritative; the table below displays rounded values for readability.

| Quantity | Value |
|---|---:|
| Outward unit-normal complete support | 0.402794181 |
| Retained head, even spins0..40 | 0.402750321591258 |
| Finite even-spin42..448 all-energy tail | 0.00004384769205095601 |
| Infinite-spin compact contribution | 5.2069035557861494e-20 |
| Infinite-spin UV and cutoff-sliver contribution | 0.000000011374230839440732 |
| Entire omitted-spin tail | 0.0000438590662817955 |
| Preserved global-margin complete support | 0.402833739 |
| Extracted published dual support | 0.4056828806060012 |
| Extracted published primal support | 0.39228187330893666 |

The canonical support is0.00288869960600124, or approximately0.7121%, below the extracted dual support in this direction. This compares numerical locations of bounds, not the rigor of the published proofs. The exact multiplier witness is unchanged from the previously audited repaired83 numerical candidate. The certificate is an outer half-plane; it does not establish attainability, optimality, or agreement with the primal curve.

The selected head uses35089 ordinary integration cells, plus the explicitly recorded endpoint/source gaps. All21 spins have a full higher-precision replay. The compact assembly rechecks34 analytic producer/replay records and71096 finite producer/replay spin-energy cells. The finite UV ledger checks112962 producer/replay cells, and the local real-cost calculation uses4096 cells across204 spins.

Measured head/source timings:

| Stage | Elapsed seconds |
|---|---:|
| Two source builds, source matching and head preflight | 54.65064334869385 |
| Original21-spin producer/replay campaign, three workers | 2553.3327860832214 |
| Sum of the42 original per-process job durations | 7433.22034406662 |
| Refined spin2 producer, reusing original interior cells | 9.060791254043579 |
| Full refined spin2 replay at384 bits | 381.39204955101013 |

The summed process durations are not wall-clock time. The spin2 refinement and tail/compact work overlapped the original head campaign, so adding these durations would not give the total campaign time. Exact subprocess commands, exit codes and timings are preserved in the head preparation/campaign/job records and compact/tail command logs. The two final assemblies passed with the commands shown above; their mathematical input files are frozen.
