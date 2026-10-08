# Final independent review of comparison, handoff and corrected certificates

The four current complete packets passed the final registry-driven exact-accounting and SHA-256 review. The review hashed 764 unique bound dependencies. The authoritative selection is `output/weak_regions_20261007/CURRENT_CERTIFICATE_SELECTION.json`; the exact reviewed source and packet hashes are in `output/weak_regions_20261007/final_handoff_review.json`. Read `output/weak_regions_20261007/CORRECTION_NOTICE.md` before using the historical outputs.

## Applicability error found and corrected

The original arithmetic-only review did not detect that the new shifted and mixture wrappers sometimes selected the Centered normalized negative-V formula beyond its finite Roy-source window. Its upper cell endpoint must obey beta squared <=79/119, equivalent to mu<=119/10. Finite interval arithmetic and negative-V checks alone do not justify the formula there. The earlier shifted review `REVIEW_SHIFTED_HEAD.md` is explicitly superseded.

At each of 256 and 384 bits, the two shifted objectives selected unsafe choices in 16 distinct spin-2 cells (30 objective-cell choices), and the mixture selected nine unsafe spin-2 choices. The raw shifted complete packets and both sets of raw spin-2 ledgers are retained unchanged as invalid historical evidence. Saved PASS strings in those files have no current evidentiary force.

The shifted correction deterministically replaces every affected choice by the corresponding original certified base-cell upper, using exact objective-shift monotonicity. The corrected assembler reconstructs the entire correction and compares it exactly. The mixture correction recomputes actual interval fields and Taylor bounds, with a positive combination of base-cell uppers as fallback, and removes every unsafe selected normalized choice. Its only completed packet already requires the repaired records. The original baseline certificates and the independently optimized fresh 269.5-degree certificate used the original guarded models and are unaffected.

## Current complete packets

| Complete packet | Outward unit-normal upper |
|---|---:|
| `complete_shift269p0404_corrected.json` | 0.007050981 |
| `complete_shift269p5_corrected.json` | 0.002133785 |
| `complete_fresh269p5.json` | 0.001905964 |
| `complete_mixture86_certificate.json` | 0.592722686 |

The final audit checks exact component sums, squared direction norms, positive lower norm bounds and outward rounding by rational inequalities. It verifies all 220 multiplier names and exact values for each NPZ witness, or reconstructs every rational positive-mixture multiplier and objective from the weights and authenticated base NPZs. It checks every declared file binding and the hashes of both current and explicitly invalidated registry entries.

For both corrected shifted packets it independently reconstructs the repaired spin-2 ledger from the raw and original base records. For the fresh head it reconstructs all 21 exact partitions from 24,246 ordinary cells and 105 analytic gaps, checks complete coverage and higher-precision replay correspondence, and sums the larger producer/replay dyadic upper for each spin. For the mixture it verifies the exact convex base total, removed original components, retained remainder and replacement components; both spin-0/2 partitions and gap sums; the exact nine corrected indices; and the identity of every unaffected cell. The correction does not silently discard any cell or remainder.

## Comparison mathematics, coefficient export and final visual check

The handoff writer now distinguishes exact-rational JSON mixture witnesses from NPZ witnesses. It verifies file hashes, objective equality and coefficient counts, then exports the JSON rational multipliers directly or uses `Fraction(float(value))` for NPZ coefficients. This preserves exact binary64 values without substituting rounded decimals. The five numerical witnesses and four completed certificate witnesses are labeled separately.

- Certified polygon cuts use the exact unnormalized objective and `support_upper_exact`; normalization never changes only one side of an inequality.
- Unit-normal displays scale the normal and support consistently. Outward rounded normalized values are display bounds, while the exact saved inequality is authoritative.
- Numerical candidates do not enter the certified polygon. The numerical envelope retains all 29 previously qualified directions, prior certified cuts and the new audited candidates.
- Branch distance uses physical axes divided by their observed dual ranges. Upper and lower branches are measured separately; corners receive no extra length from a normal-angle sweep.
- Positive normal-line gap denotes a locally weaker outer intersection; negative denotes a locally tighter intersection. Missing line intersections remain missing. Neither sign is a global support-function comparison.
- Published primal and dual curves are digitized references. The primal tolerance is a graphical check, not an author-supplied rigorous error bar. The missing right cap is excluded from length statistics.
- The repaired upper 83-degree candidate is distinguished from its rejected raw refinement. A deterministic repair or mixture certificate is not counted as another independent optimization.

The regenerated full-boundary figure was visually inspected: both primal and dual references are shown, certification and numerical candidates have distinct styles, both boundary-length difference panels have readable signs and units, and the captions state the unproved numerical remainders and polygon-evaluation distinction. The handoff lists all four current filenames, hashes, exact objective/upper decompositions and exact multiplier tables; it includes the mandatory correction explanation.

The 801-station curves remain visualization samples. A separately generated exhaustive geometry record, `comparison/geometry_continuous_audit.json`, checks every affine-envelope breakpoint on every observed digitized reference segment, with exact rational geometry and directed length enclosures. Its current reported excluded fraction is 0.5319608584883521 of observed normalized boundary length, and maximum remaining positive normal-line gap is 0.041585538004571036. These results concern the digitized coordinates and do not establish exact dominance over the published optimum. The final handoff states this distinction. This review inspected the report's geometry conventions and presentation; it does not claim a separate rederivation of that entire geometry auditor.

## Reproducible packaging

The sealed base ZIP omits `research/markdown_export_20261007/numerical_inventory.json` and the two later root Markdown exports. The addendum must therefore include the complete Markdown-export research directory and both `ALL_CONSTRAINTS_HUMAN_READABLE_2026-10-07.md` and `ALL_BOUNDS_RECHECK_AND_MULTIPLIERS_2026-10-07.md`, as well as both new campaign trees and the current master handoff. The packager documents this requirement. The canonical registry prevents historical raw PASS records from being treated as current packets.

All four packets, the corrected ledgers and correction scripts were unchanged during the final review. The final independent executable review can be rerun directly against the saved files. It checks internal exact arithmetic, applicability corrections and file bindings; it does not independently derive the physical dispersion kernels, rerun every Arb interval evaluation, establish optimizer optimality, or provide proof-assistant verification.
