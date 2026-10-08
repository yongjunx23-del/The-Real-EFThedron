# Internal review of the 259.2941-degree compact adaptation

Reviewed on 2026-10-07 for optimized witness SHA256
`321695c10783e2e0f6e6e28589f132fdb7a68497314b1cf92eb71a07612c704f`.
The folder alias `angle259` does not round the objective, which remains the
exact rational pair saved in `angle259p2941_from255.npz`.

The analytic and finite models use the same reviewed formulas as the
265-degree campaign: all nonpositive shifted crossing powers, exact
Bernstein coefficient bounds, Jensen only for powers j>=1, piecewise energy
restrictions, source polynomial division and geometric remainders, the
exact objective contribution, and all derivative rows. Only the frozen
witness, normal, source contractions and manifest paths change. Both the
analytic and finite producer/replay inequalities require V/p<-1, matching
the reported kappa=1.

The largest finite compact spin is 74, so source tables with jmax/skmax74
and N256 suffice. Their independent exact-beta/shared-moment audit passes
against the independently rebuilt same-witness head tables. The separate
66-point real-kernel comparison also passes; numerical point comparisons
do not substitute for interval domain coverage.

The same reviewed assembler checks exact energy adjacency and complementary
finite/analytic spin coverage, paired 512/640-bit finite cells and paired
400/512-bit analytic records. No new mathematical assumption or formula
change was identified. This review is internal, not external or formal,
and is not by itself an all-energy support bound.
