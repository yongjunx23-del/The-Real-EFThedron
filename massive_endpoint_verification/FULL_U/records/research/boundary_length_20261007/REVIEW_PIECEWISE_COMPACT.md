# Internal review of compact-energy refinements

Reviewed on 2026-10-07. This file supplements the already frozen
`REVIEW_NEW_ENVELOPES.md`; neither file is an external or formal review.

## Piecewise 85-degree envelope

In `compact/analytic_piecewise.py`, restricting to an exact energy interval
`[a,b]` replaces the former largest numerator factor by `(b-4/3)^3`, the
smallest positive denominator by `a-4/3`, and the worst Legendre mean and
positive-monomial denominator endpoint by `b-4`. These substitutions retain
the directions of all previously reviewed inequalities. The physical-source
branch is used only for `a>=8`; the constants `1/8`, `1/4`, and `83/60` remain
conservative on the smaller interval. The source norm denominator remains
positive. On the low-energy branch, the proof bounds the larger domain
`4<mu<=b`, so using it only on `[a,b]` is valid.

All numerator terms in the monomial-ratio envelope are nonnegative. For
larger spin the positive Legendre monomial coefficient increases, while the
positive-coefficient mean increases. Thus the selected minimum spin gives
the worst upper bound throughout the physical-source branch. The low-energy
branch checks the derivative of the linear-source-minus-crossing envelope
at its first spin. Its negative crossing derivative only becomes more
negative at larger spin, so the all-spin extension is justified.

Choosing a cutoff by scanning finite candidate integers does not reduce
rigor: each accepted cutoff carries an interval-validated inequality that
covers every larger spin. The untested smaller spins remain explicitly in
the finite bridge. The paired 400/512-bit outputs may be assembled only with
matching intervals and complementary finite-spin coverage.

## 255-degree envelope with all negative shifted crossing powers

`compact/analytic_angle255.py` uses the same restrictions and the pointwise
crossing/Jensen lemma already reviewed in `REVIEW_NEW_ENVELOPES.md`.
For `j>=1`, its retained negative coefficients contribute
`-a_j (m-10)^j`, with `m>10`. No constant-power Jensen step is taken.
On the low-energy branch the derivative is exactly

    slope - sum_{j>=1} j*a_j*r*(r*ell-10)^(j-1).

Each subtracted term is nonnegative and nondecreasing for `ell>=L` where
`r*L>10`. A negative derivative at `L` therefore stays negative. This is a
valid strengthening of the old fifth-power-only envelope. The other
forward, source, and derivative-row estimates are unchanged.

The actual exact objective is retained: `a` is included in the forward
coefficient preparation and `b` is added once to the final constant.
The `4<mu<=9/2`, all-even-spins-at-least-42 low-energy certificate overlaps
the remaining rectangles at `9/2`; there must be no uncovered interval
between this endpoint and the first finite bridge.

No formula defect was identified. Completion still depends on checking
the exact rectangle union and all finite-spin replay partitions; a positive
analytic pilot or a successful formula review alone is not a full bound.
