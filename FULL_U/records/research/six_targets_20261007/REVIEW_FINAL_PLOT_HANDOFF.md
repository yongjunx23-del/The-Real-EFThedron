# Final six-target plot and handoff review

PASS for the final artifacts and source hashes listed below. This is an internal independent saved-data and source review, not external peer review or a new physical certificate.

- `compare_results.py`: `8f9d944a3ca89872b0961e5fc6ce00272e9dc01729787570adfc7c13d11b89d8`
- `build_report.py`: `51c3bd58d2560ed7f9fa88db57ac0897c306edc805407ac39ef021e8745ed4ea`
- Root handoff `SIX_BOUNDARY_TARGETS_CERTIFICATION_2026-10-07.md`: `99e4afc2aafd9dc523ae0c0a362cbbcf0436a20894bcb16e29d17fa231273e54`
- Plot PNG: `93e8cf4972527bc796fb39cd872b0717f53227176de8072f5f27c84c188503ed`
- Plot PDF: `f59f1f43dc091eff4cc92e3ea7ce083d90ced45a9fbd78030d0f528531c7b611`

## Exact geometry and data checks

The separate script `lower/review_final_plot_handoff.py` reconstructs the intersection by exact rational edge clipping of the authenticated previous polygon with the six selected cuts. This is a different construction from the plot generator's pairwise line intersections. Its resulting set of 39 vertices agrees exactly with the saved new polygon. Every vertex satisfies every retained halfplane. The cut list preserves all previous cuts and appends exactly the six canonical completed packets, with exact objective and support values matching their authenticated certificates.

All six station supports and signed normal-line intersections were independently recomputed. Each station support strictly improves. Each requested station has a negative new local normal-line gap. The station normal for upper10% differs from the nearby newly certified witness normal; the polygon consequence, not that witness's reported normalized support, is correctly used in the station table. Likewise, station support and local normal-line distance are kept distinct throughout.

The saved continuous geometry audit reports 87.89300787056121% exclusion of the observed digitized dual length, compared with 73.8155042969071% previously. The report does not mistake the sampled 1201 stations per branch for separate optimization runs or a continuum proof of the physical kernels. The continuous geometry result concerns the supplied digitized line segments and does not establish rigorous digitization uncertainties.

## Exact witness export

All six exact objective inequalities, all component values, and every one of the 1320 retained multiplier rows (six tables of 220 rows, including zeros) were checked against the selected packet and witness bytes. Direct NPZ coefficients are exported as exact binary64 rationals. The lower/upper55 positive mixtures retain their exact rational coefficients and weights. The unchanged human-readable constraint catalogue is reproduced verbatim. The selected lower60 packet is the additional spin4 refinement, not the preserved weaker precursor.

The handoff clearly identifies the lower gap union, the lower15/45 midpoint refinement, the direct95 threshold cap and 64-cell bridge, the finite constraint-bank relaxation, all-spin/all-energy proof scope, and shared implementation limits. It does not claim independent unconstrained optimization of every plotted station or greater rigor than published dual papers. The graphical station table is labeled as floating evaluation of exact halfplanes; separately listed packet supports carry their outward rounding.

## Visual inspection

The rendered PNG was inspected directly. Both digitized published dual and primal branches appear on the full boundary plot. Previous and new certified polygons use distinct line styles/colors with a clear legend. Both difference panels span the full observed branch length from left to right, and all six requested length stations are marked. Axis labels, percentage labels and three footer notes are readable. The zero line and footer correctly identify negative gaps as locally tighter and positive gaps as locally weaker. The unknown right cap is explicitly excluded from the length comparison, while the certified polygon itself is closed.

No blocking scientific-labeling, geometry, multiplier-export, or rendering issue was found. Remaining positive-gap regions are visible; neither plot nor handoff claims uniform dominance over the published dual curve.

Machine-readable evidence and all bound source/artifact hashes are in `output/six_targets_20261007/final_plot_handoff_peer_review.json`. No physical interval kernels were rerun in this review.
