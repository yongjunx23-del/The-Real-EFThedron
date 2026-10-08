# Independent structural review of the positive-UV extension

Review date: 2026-10-07. This is a review of the mathematical construction and saved exact assembly arithmetic. It is not a fresh evaluation of every kernel, an independent implementation of all source kernels, or a formal proof-assistant verification.

## Files reviewed

- `research/full_direction_20261007/certify_positive_uv.py`
- `research/full_direction_20261007/assemble_positive_wide.py`
- `research/full_direction_20261007/consolidate_certificates.py`
- `research/full_direction_20261007/assemble_extended.py`
- `research/full_direction_20261007/combine_dense_family.py`
- `research/full_direction_20261007/refine_dense_source.py`
- `research/full_direction_20261007/CERTIFICATION_EXTENSION_NOTES.md`
- Relevant inherited estimates in `research/tail_certification/certify_uv.py`.
- `output/full_direction_20261007/certification/complete_fifteen_direction_certificate.json`.

No missing sign, factor, or spectral-domain coverage condition was found in these new constructions. This conclusion remains conditional on the previously reviewed source representations, moment-table enclosures, and physical assumptions recorded in the frozen certificates.

## Positive-support inequality

For the disk support function

\[
h(U,V)=\sqrt{U^2+V^2}+V,
\]

with \(p>0\) and \(\eta>0\), monotonicity and the Lipschitz constant two in the second argument give

\[
h(U,V)\le h(U,-\eta p)+2(V+\eta p)_+
\le \frac{U^2}{2\eta p}+2p(V/p+\eta)_+.
\]

The factor two in the positive contribution is present in the code. The all-energy real norm is used once, as \(A_{40}/(2\eta)\); it is not inadvertently restricted to or duplicated over the UV cells.

## Integrated all-spin majorant

Put \(T=\ell(\ell+1)\), \(\tau=T/(\mu-4)\), and \(Q=\mu-4/3=T/\tau+8/3\). For \(\mu\ge49\), the phase-space factor is bounded below by \(\rho_{49}=\sqrt{45/49}\). The positive Legendre Taylor expansion supplies the upper bound

\[
P_\ell\!\left(1+\frac{8\tau}{3T}\right)
\le I_0\!\left(2\sqrt{4\tau/3}\right).
\]

With the stored measure and definition of \(p\), the transformed mass is therefore bounded by

\[
p\,d\nu_\ell\le
\frac{(2\ell+1)\tau I_0(2\sqrt{4\tau/3})}
{2\pi\rho_{49}\ell^2(\ell+1)^2}\,d\tau.
\]

The telescoping identity

\[
\frac{2\ell+1}{\ell^2(\ell+1)^2}
=\frac1{\ell^2}-\frac1{(\ell+1)^2}
\]

bounds the even-spin sum starting at \(L=42\) by the sum over every integer spin, namely \(1/L^2\). Thus the cell mass used in the script,

\[
\frac{(b-a)b I_0(2\sqrt{4b/3})}{2\pi\rho_{49}L^2},
\]

is an upper bound on a \(\tau\)-cell \([a,b]\). Each particular spin's physical restriction \(\tau\le T/45\) may remove part of that cell. Integrating the larger full cell is safe because the majorant is nonnegative. The endpoint \(\tau=0\) is the \(\mu=\infty\) limit and creates no omitted energy interval.

## Envelope signs and coverage

- The compact region \(4<\mu\le49\) and the large joint variable region \(\tau\ge256\) must separately satisfy \(V/p\le-\eta\). The candidate-selection code checks both before accepting any value of \(\eta\).
- The bounded UV cells cover \(0\le\tau\le256\). Together these regions cover every even \(\ell\ge42\), \(\mu>4\).
- The signed denominator choice in the UV ratio is preserved: a nonnegative numerator is divided by the lower bound on the denominator, and a negative numerator by its upper bound.
- The modified crossing multiplier remains negative, as explicitly asserted. Multiplication by a lower bound on \(X_4/p\) therefore gives an upper bound on its contribution. Dropping its nonpositive large-\(\tau\) contribution remains safe.
- The source functions, their absolute norms, and the Roy sign certificate remain those of the bound frozen witness; only the objective and crossing multiplier change.
- A failed wider case is not accepted through the overall `PARTIAL` status. The assembled positive cases are individually required to pass. In particular, the unsuccessful default \(b=4\) trial does not enter the final fifteen inequalities.
- Spins \(0,2,4,6\), then \(8,10,\ldots,40\), then every even spin from \(42\) onward form a gapless partition.

## Exact ledger review

An independent lightweight check was saved as:

`output/full_direction_20261007/left/fifteen_arithmetic_review.json`

It is bound to the fifteen-ledger SHA-256:

`95f0fc8368dd75dcf5da2ec213d1576f38745b6d9e1cde40e453ea4d272975f5`

Using exact fractions, that check verified all direct component sums, selection of the stronger inherited/direct result, nonnegative interpolation weights, matching objective normals, equality of the stored source constants to the eligible endpoint constants, effective crossing increments, final sums, and outward six-decimal rounding. No kernel values were recomputed by that check.

The direct wider positive family uses \(d=18b/5\). Consequently the final positive interpolated family also has that effective increment, including at slopes previously certified using \(d=4b\). The consolidation records direct and effective increments distinctly; interpreting the final interpolated result as the old direct witness would be incorrect.

The fifteen results cover right-facing normals \((1,b)\), with \(-16\le b\le2\), plus valid nonnegative combinations of certified inequalities. They do not constitute certificates for the new left-sector 220-row numerical candidates or competitive certified supports around the full circle. A separate coarse global constraint is a different statement and must remain labeled separately.

## Independent-check scope

This review adds independent algebraic and exact-assembly scrutiny. Existing higher-precision replays are valuable numerical verification but generally share the source implementations. The correct claim is an interval-certified result under the stated assumptions, with the recorded replay and review scope. It is not a claim that every representation was independently rederived in a separate numerical implementation or formally verified.
