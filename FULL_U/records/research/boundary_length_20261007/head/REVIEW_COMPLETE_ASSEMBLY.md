# Review of the complete 85°/255° certificate assembler

Reviewed `research/boundary_length_20261007/assemble_complete.py` and the
head, UV, compact, and split real-norm packet schemas on 7 October 2026.

No defect was found in the domain coverage or upper-bound combination.
The finite heads cover every even spin 0–40 at every energy; the two tail
spin ranges 42–448 and 450–infinity meet without a missing even spin.
The compact coercivity domain is `4 < mu <= 49`, and the UV coercivity
domain is `49 <= mu < infinity`. The exact beta cutoff for the real norm
lies below `sqrt(45/49)`. Thus its remainder contains a small compact
sliver as well as the UV domain. Using the smaller compact/UV coercivity
constant on that entire remainder is conservative and covers the sliver.

The combination uses `H(U,V) <= U^2/(-2V)` for negative V. Accordingly each
real norm is divided by twice a positive coercivity constant. All
components are nonnegative. Original and replay upper bounds are retained
by maxima before exact rational summation. The final six-place decimal
ceiling is checked against the exact sum; numerical normalization and the
graphical published comparison do not enter the exact inequality.

The assembler rechecks head record hashes, exact objectives and partitions,
all even-spin coverage, linked higher-precision replays, cell/gap sums,
source rebuild hashes, and source Gauss-degree sufficiency. Split real-norm
packets are checked for complete ordered partitions, region labels,
precision ordering, nonnegative bounds, per-cell maxima, and exact sums.
Compact and UV packets are accepted through their completed schemas and
recursive dependency manifests; their own assemblers perform the detailed
energy/spin rectangle or tau/eta coverage checks. These checks do not
replace a mathematical review of the underlying field and remainder code.

The review identified three imported real-norm helpers missing from the
outer manifest. The parent added bindings for the correct paths:

- `research/all_direction_certification_20261007/alternative/rational_head.py`
- `research/all_direction_certification_20261007/alternative/exact_bernstein.py`
- `research/lower_g0_upper_g2_20261006/matrix_real_tail_gram.py`

The correction was verified in the source before any completed certificate
was assembled. It does not alter the computed bounds or require changing
the frozen real-norm packets. `certify_real_tail.py` is already covered by
the finite-head dependency manifest.

This is an internal code and formula review, not an independent kernel
implementation or formal verification. The scientific comparison still
depends on matching the PDF's analytic assumptions, and its reference
curve is graphically extracted rather than exact author optimizer output.
