# Left180 compact coercivity

The fixed functional has objective exactly `(-1,0)` and witness
`output/full_direction_20261007/left/angle180_final.npz`, SHA256
`e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
All source moments are newly contracted with this witness. The robust-g2
moment values are not reused.

The target component is `V_l(mu)/p_l(mu)<-10` for every even `l>=42` and
every `4<mu<=49`. It has three overlapping pieces:

* Every even `l>=42`, `4<mu<=5`: the sharp analytic low-energy envelope.
* Even `l=42,44,...,98`, `5<=mu<=49`: fresh interval evaluation over exact
  gapless rational partitions.
* Every even `l>=100`, `4<mu<=49`: the sharp analytic high-spin envelope.

The analytic pieces are independently evaluated at 400 and 512 bits. The
finite piece uses a 512-bit producer and a fresh 640-bit replay of every
accepted spin/energy cell. The integer-spin seam between 98 and 100 and
the energy seam at 5 are covered. The physical threshold is open.

## Identical kernel algebra, new functional

The finite evaluator is a separate copy of the frozen, reviewed robust-g2
evaluator. Its only changes are the witness and cache paths, default last
spin 98, and an assertion of the exact objective `(-1,0)`. The algebra of
the scaled Legendre polynomials, exact removable-pole division, Roy PV
orthogonality, tapered-strip PV endpoint expressions, derivative rows,
and combined crossing numerator is unchanged. See
`research/competitive_20261007/robust/COMPACT_BRIDGE_PROOF.md` for those
formulas. The source norm inputs and every source moment are replaced by
the newly bound left-witness data.

In particular, `Rstar_cc` retains its full coefficient `k(k-1)`. The X4
and higher crossing terms are combined over `Q^6 Rstar` before interval
evaluation. The constant objective term in this numerator is now zero;
the `-g0` objective is included in the forward subtraction coefficient
`alpha`. All 220 row coefficients are retained and interpreted as exact
binary64 rationals.

The new source tables have angular order through 98, 256 geometric
moments, and Roy spectral polynomial order through 98. Certified Gaussian
nodes and weights enclose degree-exact polynomial integrals at 512 bits.
The explicit order-256 geometric residual is still included in every
field evaluation. The separately certified combined strip norm is used
for each individual strip, a conservative choice. A positive interval
lower bound for `Rstar` is required before division.

An independent floating implementation was checked at spins 42 and 98
over energies 4.5 through 49 before launching interval production. This
is a diagnostic implementation check, not the justification for a
continuum claim. The separate moment audit checks all 297 transfer
moments against exact beta integrals and 56,034 shared moment enclosures
against independently constructed head tables. The moment builder's
polynomial exactness theorem is the primary integration justification;
overlap alone is not treated as a certificate.

## Analytic spin closure

The same exact crossing polynomial has a fifth-degree negative tail
after shifting the Taylor index by ten. Exact rational Bernstein bounds
control its coefficients throughout `0<=1/Q<=3/8`. The finite low-index
part is bounded separately. Jensen's inequality bounds the fifth power
using the positive Legendre Taylor mean.

The source positivity certificate is the left witness's exact Roy
certificate. Its absolute source norms, signed K2 moments, and endpoint
derivative bounds are also left-witness quantities. The complex-transfer
bound `|omega|<4/3` and source geometry are unchanged. The positive
Legendre monomial bound controls forward terms uniformly in spin and
energy. In the low-energy region, the resulting spin envelope has a
negative derivative at the first retained spin and remains decreasing.

For the left functional there is no added `+1` upper-g2 objective term.
At the declared endpoints, the saved bounds are approximately
`-394627.76` for the all-spin42 low-energy component and `-501.05` for the
all-spin100 middle-energy component. Acceptance uses exact directed
bounds below `-10`, not these printed approximations.

## Interval geometry and replay

The prior g2 partitions supply initial interval endpoints only. No g2
source moment, multiplier, kernel value, or inequality is reused as a
left-witness bound. Every seed cell is evaluated anew, and any failure is
bisected. Each spin's final cells must form an exact gapless partition of
the claimed energy interval.

For parallel production the energy interval is divided into the six
half-unit pieces of `[5,8]`, followed by `[8,12]`, `[12,24]`, `[24,36]`,
`[36,48]`, and `[48,49]`. The fresh replay groups identical intervals
across spins only to share computation. It verifies every original
spin/interval pair, exact partition equality, finite values, and strict
upper bounds below `-10`.

The final assembler checks all exact seams, spin lists, witness and
dependency hashes, and paired-precision component claims. The result is
only a compact-domain coercivity component. A lower-g0 inequality also
requires the finite-head integration, the real-part norm, and all UV
regions. Shared-code replay is not an independent implementation or a
formal proof-assistant certificate.

## Reproduction from frozen source inputs

```text
python research/competitive_20261007/left/build_compact_moments.py --kind roy
python research/competitive_20261007/left/build_compact_moments.py --kind strips
python research/competitive_20261007/left/certify_analytic_compact.py --bits 400
python research/competitive_20261007/left/certify_analytic_compact.py --bits 512
python research/competitive_20261007/left/run_compact_slices.py
python research/competitive_20261007/left/replay_compact_parallel.py --input output/competitive_20261007/left/compact/slice_8_12_512.json --out output/competitive_20261007/left/compact/slice_8_12_640_parallel.json --jobs 6
python research/competitive_20261007/left/replay_compact_parallel.py --input output/competitive_20261007/left/compact/slice_12_24_512.json --out output/competitive_20261007/left/compact/slice_12_24_640_parallel.json --jobs 6
python research/competitive_20261007/left/merge_compact_slices.py
python research/competitive_20261007/left/assemble_compact_coercivity.py
```

The two moment builders may run concurrently. The slice runner waits for
their manifests, performs its independent floating preflight, and then
runs fresh 512/640-bit computations. The final two commands require all
pieces to have passed. The selected 640-bit replays of `[8,12]` and
`[12,24]` use the explicit six-process commands, which can overlap the
slice runner; their redundant serial replays are not selected. The
independent moment-audit packet and the
left/refined source cache are prerequisites checked by the final hash
manifest. Existing frozen results should be replayed in a disposable
copy if bit-for-bit preservation is required.
