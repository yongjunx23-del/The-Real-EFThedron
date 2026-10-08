# Final consistency review of the boundary-length campaign

Reviewed 2026-10-07, read-only except for this new note. No optimization, interval-kernel rerun, or frozen-input edit was performed. This is an internal independent structural and saved-record review, not an independent implementation of the entire proof.

## Result

No new mathematical or bookkeeping defect was identified in the completed 85°, 255°, and 265° ledgers, `assemble_complete_v2.py`, `plot_lower_progress.py`, or the current `build_handoff.py`. The 259.294118509649° final packet was not yet present during this review; its extension was reviewed structurally, so this note does not certify its still-pending final head/assembly result.

The exact rational component sum, head/tail decomposition, stored normal squared, and upward six-decimal rounding agree for all three completed ledgers. Every listed dependency hash matched its current file: 230 manifest entries for 85°, 235 for 255°, 227 for 265°, and 529 distinct files in their union. Witness hashes also matched the NPZ bytes.

| Nominal angle | Complete outward bound | Graphical dual support | Graphical primal support | Relative difference from dual |
|---|---:|---:|---:|---:|
| 85° | 0.453860 | 0.457271030269 | 0.417705556599 | -0.745953722% |
| 255° | 1.005497 | 1.039322823180 | 0.928910672718 | -3.254602172% |
| 265° | 0.122076 | 0.124251322407 | 0.096464443809 | -1.750743867% |

These ratios were independently recomputed with the saved exact direction divided by its norm and the benchmark support API. All three bounds remain above the graphical primal support. The table's six-decimal certified numbers govern the exact stored rational objectives; the normals are numerically unit length, but are not defined as mathematically exact trigonometric values. `support_normalized_float` and the displayed percentages are comparison diagnostics, not additional outward interval certificates.

## Proof assembly and signs

The head covers precisely the even spins 0 through 40. The omitted-spin pieces cover finite even spins 42 through 448 and all even spins from 450 upward. Compact energy is 4 < mu <= 49; UV energy is mu >= 49. Common boundaries have zero integration measure and there is no spin gap.

The support convention is H(U,V) = sqrt(U²+V²)+V. For p > 0 and V/p <= -kappa < 0, rationalization gives H <= U²/(2 kappa p). The assembler uses that factor of two consistently. Its real-norm split occurs at a dyadic beta below sqrt(45/49); the remaining small compact-energy sliver is conservatively charged using min(compact kappa, UV kappa), rather than assuming it belongs wholly to the UV domain. The objective signs are taken from the exact NPZ normal, including the negative g0 and g2 components of the lower-sector witnesses.

Head checks bind witness, exact objective, source manifest and code dependencies; check the exact partition of beta in [0,1], source degree-exactness conditions, producer/replay sums, higher replay precision, and retain the larger per-spin upper bound. Where parallel replay chunks occur, their index coverage and single inclusion of endpoint gaps are checked. Compact and UV domain coverage is supplied by their separately assembled, hash-bound component packets. This is a layered proof review, not a fresh replay of those kernels here.

The v2 diff only adds the accepted 259.294118509649° witness alias and the existing constant-negative spin-2 cap with its exact directed endpoint-sign record. It retains the actual exact normal; it does not substitute the direction at 259°. The beta cap [0,1/256] has strictly negative upper endpoints at both recorded precisions. No negative sign is inferred from rounded Arb display text.

## Constraint inventory and plots

The actual NPZ names match the handoff's 220-row inventory: 39 forward H2–H40 rows; 7 crossing rows X4, X5, X6, X7_1, X7_2, X8_1, X8_2; 60 Roy rows; 64 tapered-strip rows; 36 forward-transfer derivative rows; and 14 rational cosine/sine rows. The Roy parameters and strip windows agree with the reported canonical geometry. Guards on these coefficients do not create extra physical assumptions or substitute finite sampling for a continuum proof.

The plot constructs the certified polygon from exact rational halfplanes and accepts only completed ledgers. Numerical candidate lines are placed in a separate, explicitly uncertified polygon. Dense angles evaluate support functions of these polygons; they are not new optimizations or radial distances. The normalized-arclength plot uses local segment normals and reports absent line intersections as gaps, without sweeping normals through corners. The published dual and primal curves are both plotted.

The benchmark support is the maximum over the extracted upper/lower vertices, with the documented isolated right-tip marker included for dual comparisons. This is a graphical estimate. The missing right cap receives no invented boundary segment or arclength. It is not active in the three checked comparisons; this does not turn the working extraction tolerance into a rigorous error bound on the authors' optimization. No claim of universal superiority over all published duals is justified without matching analyticity and amplitude assumptions.

The handoff now includes the previously missing 259 compact cutoff, constant-cap explanation, and v2 reproduction command. Final plots and tables should be regenerated only after its completed packet exists. No remaining essential edit was found in the reviewed handoff version.

## Limits of this review

The saved proof is an internal computer-assisted interval certificate conditional on the dispersion/source identities, amplitude class, exact canonical parameters, and reviewed analytic tail lemmas. Higher-precision replays mostly share their mathematical implementation; separate source identities, algebraic specialization, selected field comparisons and structural reviews provide additional checks but are not formal proof-assistant verification or external peer review. The graphical dual/primal comparisons retain their separate, explicitly weaker provenance.

## Addendum: completed 259.294118509649° packet

The final packet became available after the initial review above. Its saved-record consistency review **passes**, superseding the earlier pending-result qualification for this direction. The reviewed file is `output/boundary_length_20261007/complete_angle259_certificate.json`, SHA256 `cdb4e2c5d553d2b6b2158e1752b6683dec47c0579585e2c7ad0e9224d2338368`.

All 229 dependency hashes match the current files. The witness hash is `321695c10783e2e0f6e6e28589f132fdb7a68497314b1cf92eb71a07612c704f`. The head, compact and UV packets use the same witness and exact normal

`(-6692978864176391/36028797018963968, -4425208767245823/4503599627370496)`.

The exact rational component sum is approximately 0.4630915997746932 and its exact upward six-decimal ceiling is **0.463092**. Its components are head 0–40: 0.4607996033793012; finite compact tail: 0.00020585754028615728; finite UV plus sliver: 0.0020857329851238287; infinite compact tail: 2.1479917269285844e-18; infinite UV plus sliver: 4.0586998203212903e-7. The tail sum, normal squared and normalized diagnostic all agree with the stored ledger.

The domain checks confirm even spins 0–40 in the head, 42–448 in the finite tail, and every even spin from 450 onward in the infinite tail, with compact and UV meeting at mu=49. Compact kappa=1, finite UV kappa>0.02, and infinite UV kappa=0.02 have the required positive sign. Both exact directed spin-2 cap upper endpoints on beta in [0,1/256] are strictly negative; no sign is inferred from rounded display strings.

The graphical dual support is 0.47818495173466086, and the graphical primal support is 0.4195951391987455. The certified bound is **3.156300022% smaller than the graphical dual estimate**, while remaining above the graphical primal. The missing-right-cap sensitivity flag is false here. These benchmark numbers and their precision retain all graphical-provenance and shared-assumption caveats stated above. This addendum checked saved exact records and hashes; it did not rerun interval kernels or optimization.
