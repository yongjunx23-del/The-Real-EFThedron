# Internal review of energy-dependent finite-tail costs

Reviewed `banded_finite_cost.py` on 2026-10-07 against the earlier reviewed
`real220_banded_coercivity.py` construction. No new integration is used.

For every saved beta cell, A bounds the sum over even spins 42 through 448
of the integrated nonnegative real density U squared divided by p. A common
positive kappa valid for every such spin and every point of the cell implies
an integrated support cost at most A/(2 kappa). Choosing the minimum of the
actual local certified margins is conservative; it does not require those
margins to be equal or use interpolation between sampled values.

The map is exactly tau=ell(ell+1)(1-beta squared)/(4 beta squared), decreasing
in beta. The compact/UV seam is compared through beta squared=45/49. A cell
straddling the seam includes both compact and UV bounds, so the small real-norm
cutoff sliver is covered. UV intervals touching either endpoint are included;
large-tau pieces above 256 use their separate analytic proof. For spins whose
physical UV domain ends below 256, the saved finite cover ends at ell(ell+1)/45.
The complete even-spin list is 42,44,...,448 (204 spins).

The driver takes the larger directed upper from each original/replay UV cell,
checks that it is negative, and negates it for a safe margin. Each real-norm
cell retains max(original,replay). The reciprocal weight is rounded UP to
40 significant binary bits, with exact rational inequalities checked. Products
and sums are then exact dyadics. The final complete-bound assembler must replace
the old two finite-tail terms with this single finite cost, not add both.
The infinite-spin tail remains a separate contribution.

Reviewed checks include exact paired partitions, spin and energy coverage,
witness/objective identity, source hashes, increasing precision, finite-width
seams, and exact cell-cost sums. The saved ledger makes the interval choices
auditable. This is an internal formula/code review. It does not replace the
underlying interval proofs, establish their physical assumptions independently,
or constitute formal verification or external review.
