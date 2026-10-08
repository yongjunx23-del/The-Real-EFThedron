# Independent internal review of the final selection and banded-tail assembly

Reviewed `research/boundary_length_20261007/assemble_complete_v3.py` and
`research/boundary_length_20261007/banded_finite_cost.py` on 2026-10-07,
focusing on their changes from the previously reviewed complete assembler.
No mathematical defect was found in the new composition.

The selected-head logic checks the exact witness/objective and final
record paths through the existing head validator. The model-selection
metadata binds the selected budget, cap proof, and exact endpoint evidence.
The constant-cap endpoint domain is checked against the actual selected
spin-2 threshold gap, so the 268° interval `[0,1/1024]` is not silently
treated as the earlier `[0,1/256]` interval. The stopped quadratic-cap
checkpoint is not used when the caller supplies the alternate final head
summary. Exact negative endpoints, not coarse Arb display strings, are
tested.

For the banded finite tail, the saved real-norm density is nonnegative and
already bounds the aggregate over all 204 even spins 42–448 on each beta
cell. Using the minimum certified coercivity across every such spin and
every touching energy interval gives a valid common weight. Therefore
`A/(2*kappa)` bounds the support cost on the whole beta cell.

The exact decreasing map
`tau = ell*(ell+1)*(1-beta^2)/(4*beta^2)` is used. Cells crossing the
compact/UV seam include both bounds. Finite tau cells touching endpoints
are included on both sides where present; tau greater than256 uses the
separate large-tau proof. Its minimum spin108 is sufficient because a
physical UV interval with tau greater than256 requires ell at least108.
The finite tau cover ends at `min(256,ell*(ell+1)/45)`, exactly matching
the physical energy domain beginning at49.

The cell kappa uses the larger production/replay upper for V/p, negated;
the real integral uses the larger production/replay upper for A. The
reciprocal weight is rounded upward by exact rational comparison. All
cell products and sums are exact. The complete assembler independently
recomputes the saved interval-choice ledger and replaces the two old
finite-tail terms with the single banded term, avoiding double counting.
The infinite-spin contribution remains separate.

This review concerns formula composition, coverage, and bookkeeping. It
does not replace the underlying interval proofs or independently establish
the physical assumptions of the dispersion and unitarity framework.
