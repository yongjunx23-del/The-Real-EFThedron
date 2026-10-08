# Directional certification extension, 2026-10-07

All new code is in `research/full_direction_20261007/`; all certificate outputs
are in `output/full_direction_20261007/certification/`. Frozen originals are
read-only dependencies. The base witness is the150-row `bbbc5990...` NPZ.

## Final consolidated result

`complete_fifteen_direction_certificate.json` is the authoritative new ledger.
It improves the upper endpoint to `g0<=2.705955` (exact upper total
approximately2.705954615160051) and contains fifteen complete rational
directions from slope-16 through+2. Examples are
`g0-16g2<=1.766916`, `g0-2g2<=2.585048`, `g0+(1/2)g2<=2.737537`,
`g0+g2<=2.769119`, `g0+2g2<=2.833035`. All displayed decimals are outward.

The ledger distinguishes the direct component bound from an improved exact
positive combination. In particular, interpolation between the preceding4b
family and new3.6b family can change the effective X4 increment at the same
slope. `crossing_increment_exact` and `direct_crossing_increment_exact`
describe the direct component record; use
`certified_effective_crossing_increment_exact` to reconstruct the possibly
interpolated certified functional. All other multipliers stay fixed.

## Direct directional family

Use the exact normal `(1,b)`, keep all original multipliers, and change only
the X4 multiplier by `d=0` for negative slopes and `d=(18/5)b` for positive
slopes. Thus `U` and its all-spin real norm remain unchanged, whereas
`V_new/p=V_old/p+b+d chi_l`. The threshold spin2 change is `-b/80`, since
`chi_2(4+)=-9/32`. This is a smaller positive repair than the preceding4b
family, while retaining a strictly negative threshold shift.

The first extension consists of slopes `-2,-4,-8,-16,3/10`. Complete results
and exact component totals are in `complete_extended_certificate.json`.
Heads0,2,4,6 are replayed at280/340bits, finite8..40 uses the previously
reviewed per-spin real-norm and spectral vector tables, and the all-spin
tail is replayed with the actual repair. Negative slopes strengthen the old
global negative V/p bound automatically. These are direct fixed-witness
bounds, not mere interpolation of the older seven slopes.

## A tail lemma allowing positive V

Uniform negativity of V is sufficient but unnecessary. For any exact
positive eta and positive p, the disk implies

```
h(U,V)=sqrt(U^2+V^2)+V
      <= U^2/(2 eta p) + 2 p (V/p+eta)_+.
```

One proof uses `|Ur|<=U^2/(2 eta p)+eta p r^2/2`, `r^2<=2b`, and
`0<=b<=2`. Equivalently h is increasing and2-Lipschitz in V and can be
compared with `h(U,-eta p)<=U^2/(2 eta p)`. This is an inequality for the
same exact disk support, not a new physical assumption.

Put `k=l(l+1)`, `tau=k/(mu-4)`, and restrict to `mu>=M=49`, even `l>=L=42`.
Then `rho>=sqrt(45/49)`, `Q>=k/tau`, and positive Legendre Taylor
coefficients give `P_l(1+8tau/(3k))<=I0(2sqrt(4tau/3))`. Since
`|dmu|=k/tau^2 dtau`,

```
p_l dnu_l <= (2l+1)/k^2 * tau I0(2sqrt(4tau/3))
             /(2 pi sqrt(45/49)) dtau.
```

Summing the larger set of all integer spins telescopes exactly:

```
sum_(l>=L) (2l+1)/(l^2(l+1)^2)
 = sum_(l>=L) [1/l^2-1/(l+1)^2] = 1/L^2.
```

Thus on a joint-tau interval `[a,b]`, the complete even-spin p-mass is at most

```
(b-a)b I0(2sqrt(4b/3))/(2 pi sqrt(45/49) L^2).
```

The upper Riemann estimate is valid because `tau I0(2sqrt(4tau/3))` is
nonnegative and increasing. It may integrate over extra mu values below49
after the change of variables, which only enlarges this nonnegative bound.

`certify_positive_uv.py` obtains a directed upper `v_c` on V/p in every
old validated UV tau cell, using the exact modified X4 coefficient. It then
sums exactly

```
A40/(2 eta) + 2 sum_c max(v_c+eta,0) * p_mass_upper_c.
```

The compact region4<mu<=49 and the unbounded tau>=256 region must obey
`V/p<=-eta`; those sign conditions are explicitly checked. Hence the sum
covers the entire spectral domain and all omitted spins, not just sampled
positive islands. Eta is selected from an explicit rational finite list;
each passing choice independently gives a valid bound. No external g2 cap
or numerical estimate of an omitted spin sum enters this argument.

This extends the passing tail cases to b=1/2,1,2 at d=9/5,18/5,36/5.
Their tail upper costs are approximately0.00017887431,0.00066719612,
0.00243267111. The b4 sufficient compact envelope fails; this says nothing
about physical infeasibility or the actual fixed-witness support.

## Finite-block refinement for the wider positive family

The old per-spin spectral partition was sufficient at b1/2. At b1 and2,
only14 and39 old coarse cells, respectively, failed the affine sign test.
`refine_positive_finite.py` retains8765 already valid cells and replaces the
39 failures using117 fresh320-bit interval evaluations. The result has
8843 cells, covers mu4..101 and compactified w0..1/97 without gaps, and
checks every spin8,10,...,40. It stores raw old-witness V and X4 bounds;
`combine_positive_finite.py` explicitly applies the actual d=(18/5)b.

No source/kernel/moment approximation is changed. The refinement only
reduces interval overestimation. `assemble_positive_wide.py` verifies the
exact affine signs again, checks the head partitions and exact totals,
checks all positive-UV correction sums, and binds hashes. Heads are run at
320/360bits; positive-UV certificates at280/360bits. These are fresh
invocations of the same reviewed implementation, not independent code.

A further bounded refinement, `refine_dense_source.py`, forces spectral
interval widths at most1/32 in the physical source region12<=mu<=48. It
retains8564 prior cells and replaces279 coarse cells using2173 fresh
320-bit interval evaluations, producing9790 spectral leaves. Its metadata
field `initial_failed_cells` counts all cells selected for refinement,
including cells selected solely by width; it does not indicate a physical
sign failure. Runtime was71seconds. `combine_dense_family.py` merges that
partition with19790 existing real-mass cells for all15 exact family cases
(70seconds). `consolidate_certificates.py` substitutes any improved finite
contribution into the prior complete certificates and then checks exact
positive combinations. This improves the finite g0 block from
0.00263719134027 to approximately0.00197579924052.

## Weak whole-plane baseline

`certified_baseline.py` verifies the existing Gram packet's hashes and
combines its335 rational cuts, g2 positivity, and the seven earlier complete
directional inequalities. It gives72 rational-normal support certificates
by exact nonnegative combinations. Floating-point LP is used only to select
the pair of cuts; every combination is checked with Fraction arithmetic.
The global box `-81087<=g0<=2.706617`, `0<=g2<=16556` independently proves a
finite bound for every real direction. This is a weak derived baseline,
not competitive certification around the whole published island.

The Gram closure uses the first transfer derivative and consequently the
stated local uniform convergence needed for differentiating its null.
None of the new competitive right-side fixed150 certificates requires
this separate Gram closure. Plot labels for the72 baseline points are
approximate angles; their actual normals are exact stored rational pairs.
`whole_plane_baseline_updated.json` uses the new fifteen-direction ledger
instead of the older seven points; the other sides remain deliberately weak.

## Reproduction

Use the workspace Python with python-flint0.9.0 available from
`tmp/certification_deps`; the scripts insert its path. NumPy and SciPy are
needed for the baseline and frozen numerical helper imports.

```
python research/full_direction_20261007/certified_baseline.py
python research/full_direction_20261007/extend_fixed_family.py tail --name tail_extended_18_5.json
python research/full_direction_20261007/extend_fixed_family.py finite --name finite_extended_18_5.json
```

For each spin l in0,2,4,6, run `head_family_extended.py` with
`--spin l --slopes=-2,-4,-8,-16 --precision280` and the matching
`negative_head_spinl.json` output, then slope3/10 to
`positive_head_spinl.json`. Replay all five slopes at340bits to
`replay_extended_head_spinl.json`. Exact command options use spaces, e.g.
`--precision 280`. All `--out` paths must include the new output directory.
Then run `assemble_extended.py`.

For the wider family:

```
python research/full_direction_20261007/certify_positive_uv.py
python research/full_direction_20261007/certify_positive_uv.py --slopes=1/2,1,2 --precision 360 --output positive_uv_tail_replay.json
python research/full_direction_20261007/refine_positive_finite.py
python research/full_direction_20261007/combine_positive_finite.py
```

Run head slope1/2 at320bits to `half_head_spinl.json`, slopes1,2 at320bits
to `wide_head_spinl.json`, and all three at360bits to
`wide_replay_head_spinl.json`, for each l. Then run
`assemble_positive_wide.py`. Every final assembler checks its input hashes;
do not edit source after creating certificates without regenerating the
affected records. Existing components are frozen dependencies, not files
to overwrite.

To reproduce the final additional tightening and fifteen-point ledger, run:

```
python research/full_direction_20261007/refine_dense_source.py
python research/full_direction_20261007/combine_dense_family.py
python research/full_direction_20261007/consolidate_certificates.py
python research/full_direction_20261007/certified_baseline_updated.py
```

Reproduction timings are retained in each JSON. The main head batches took
approximately17–128seconds per spin depending on case count and precision;
the positive-UV proof/replay took30–39seconds. These elapsed times reflect
concurrent processes and should not be added as a serial wall-clock estimate.
