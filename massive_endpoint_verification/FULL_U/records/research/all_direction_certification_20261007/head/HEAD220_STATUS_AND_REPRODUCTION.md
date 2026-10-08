# Rigorous finite-spin head work for the 220-row witnesses

This document concerns **components of fixed-witness bounds**, not a completed
competitive bound. A finite-spin head result cannot be added to a numerical
spin-tail estimate and called certified. The tail branch must cover every
remaining even spin and every energy, using the same witness hash.

## Fixed inputs

* Upper g2: `output/full_direction_20261007/right/g2_guarded192_refined.npz`,
  SHA256 `825dcff5caca223ae757bc027b2e1438d80e9878f4c2801db80e2a5790452ca7`.
  Its numerical head0..40 is1.045147553850433. Direction is exactly(0,1).
* Left: `output/full_direction_20261007/left/angle180_final.npz`, SHA256
  `e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
  Direction is exactly(-1,0).

The220 rows are all retained:39H2..40,7X4..8,60Roy(B=11.9,J=0,2,4,6,n=0..14),
64 tapered strip rows(two windows, j=0..3,n=0..7),36 transfer derivatives
(orders1..3,H_j2..12 plusG_regular), and14 cos/sin rational generators.
See `research/full_direction_20261007/right/RIGHT220_REPRODUCTION.md` for the
full numerical constraint inventory and prior rejected quadrature candidates.

All stored multipliers and direction components are exact binary64 rationals.
Labeled rational parameters are interpreted as exact decimal rationals;
canonical strip T is exactly4/3, despite the shortened stored name.

## Completed records (additional records may appear in the directories)

* `head/g2/head_l0.json`: rigorous all-energy spin0 upper0.89996298971257,
  1295 cells,60.23s. `head/g2/head_l0_replay.json` separately recomputes every
  cell and all5 gaps, verifies exact coverage[0,1], and passes ceiling0.9.
  Replay time37.44s.
* `head/g2/taylor_tight/head_l0.json`: improved rigorous spin0 upper
  **0.880589838585546**, total102.66s. Numerical spin0 is0.876109808358857.
* `head/g2/remaining/head_l2.json`: rigorous all-energy spin2 upper
  0.148013813184739,4000 cells,342.41s; desired0.1111314 was not reached.
* `head/g2/taylor/head_l2.json`: improved all-energy spin2 upper
  0.126607141294514,4000 cells,253.05s; desired0.1111314 was not reached.
* `head/{g2,left}/threshold_caps.json`: rigorously negative V/p on a finite
  threshold neighborhood for every even spin2..40, including explicit
  support-mass upper bounds. These are not merely threshold limits.

All `head/...` paths in this document are relative to
`output/all_direction_certification_20261007/`.

The ongoing strongest integration branch is `head/g2/centered/`. Its outputs
must be read, not inferred from an optimizer or an intermediate console line.
The exact dyadic endpoints in each JSON are authoritative; displayed decimals
are convenience approximations. A record whose `target_met` isfalse can still
provide the displayed upper bound if every cell has a finite proved enclosure.

## Mathematical enclosure

Use beta=sqrt(1-4/mu), e=1-beta^2, and scaled fields U/e^2,V/e^2. The
single-spin integral is4(2ell+1)/pi times the integral from0 to1 of
H(U/e^2,V/e^2), where H(U,V)=sqrt(U^2+V^2)+V.

`finite220.py` adapts the exact source-polynomial constructor and enclosure
methods of `research/head_certification_lowspin.py`. The original frozen code
is neither changed nor used with its hardcoded old moment-cache constructor.
The new constructor explicitly loads witness-matched tables from its argument.

For finite energy, principal-value source integrals are exact polynomial/log
expressions. The removable crossed pole is canceled by exact polynomial
division. Other source denominators use finite geometric expansions with
explicit absolute remainders. No source quadrature error is silently dropped.

For mu>=128, the source kernel uses a32-term series in x=1/mu with remainder
bounded by4*norm*Tmax^j*B^2*(B*x)^32/(1-B*x), with the corresponding factor
for the other source term. This includes beta=1 and thus energy infinity.

Outer cells use both interval Riemann support bounds and the convex
interpolation inequality. With cell widthh and analytic second-derivative
Taylor coefficientsU[2],V[2], the interpolation remainder is
h^3/8*(2|U[2]|+4|V[2]|)+4h*source_error. All operations use256-bit Arb.

Small neighborhoods of source endpoints mu=11.9,12,24,48 are omitted from
bulk cells and bounded separately by exact analytic log-integral bounds.
Their rational beta endpoints enclose the irrational endpoint and are checked.
The partition includes every beta from0 to1 with no omitted sliver.

## Exact forward and derivative rows

`forward220.py` handles harmonics, both rational generator types, higher
crossings and all36 derivative rows. The inner transfer derivatives use exact
Taylor algebra through order3, not finite differences. For ell>=2, substitute
t=beta^2*tau before extracting coefficients. Factor
(1-omega)^(2p-2)/beta^(2p-2), leaving an explicit beta^-2. A complex value is
represented as(real,imag/beta), so beta*U is finite at threshold.

The G_regular subtraction is retained separately and normalized by p using
the exact Legendre numerator. Its threshold terms do not get discarded.
`forward_combined.py` combines all coefficients of each derivative order
before interval evaluation, expressing the polynomial in omega as an exact
Taylor polynomial in omega-1. This changes the enclosure width, not the
functional or row inventory. It is implemented for the known j<=12 bank.

`forward_validation.json` compares independent numerical kernels for
ell0,2,4,14,40 at9 energies. Additional full-field comparisons for ell0 and2
use21 energies from4.0001 to1e6. These diagnostics help catch implementation
errors but do not themselves constitute a proof. Combined/uncombined Arb
evaluations overlap at those independent point checks.

## Threshold certificates

For the basic cap, rigorously prove V/p<0 on beta in[0,b]. If Ub bounds
beta*U/e^2, then H<=U^2/(-2V), and the remaining inverse-p factor has an
explicit beta^(2ell-2). Integrate that power exactly. This handles the genuine
U~1/beta behavior of derivative rows, which the older bounded-U cap cannot.

The improved `finite_taylor.py` cap writes V/p as an even analytic partA(beta)
plus explicit beta^(2ell)*log(beta) source terms and a uniform geometric
error. Arb encloses A''/2 on the whole cap. Because A is even, A'(0)=0
identically. Bound beta^(2ell-2)|log beta| by its endpoint maximum on this
small interval. For ell2 and4 this proves
V/p<=constant_upper+quadratic_upper*beta^2, with both coefficients negative.
Then use the quadratic term in the support denominator and integrate
beta^(2ell-4) exactly. For the repaired g2 witness atb=1/256:

* ell2: constant<=-5.25617621776e-7, quadratic<=-6.55791638184;
  cap mass<5.08744768390e-5.
* ell4: constant<=-2.23039153886e-8, quadratic<=-38.5036193832;
  cap mass<7.99360876737e-14.

Higher finite spins already have comfortable constant threshold margins and
use the simpler cap. The source moment remainders are included in these signs.

`finite_centered.py` further evaluates normalized V/p and the real field on
bulk cells as the midpoint value plus midpoint first derivative times the
radius, plus an interval second-derivative remainder. This preserves the
cancellations at the midpoint instead of subtracting large interval ranges.
It is used only in the low-beta region where the Roy source is active and
the physical strips are inactive. These are enclosure improvements, not new
physical assumptions or sampled sign checks.

## Source tables

`build_moments.py` calls the existing pure rigorous builders. All source
integrands are polynomials; the Gauss rule order is checked against their
degree, and nodes/weights are enclosed with Arb. Thus this use of a Gauss
rule is exact polynomial integration enclosed by intervals, not empirical
quadrature convergence.

For the N128 tables used for threshold and later head work:

* Moment jmax74, kmax128, Roy skmax128; UV jmax40,N32.
* Roy order109: degree bounds216 in energy and207 in angle.
* Strip order103: degree bounds136 in energy and204 in angle.
* g2 build20.33s; left build20.90s.
* Earlier spin0 record used N64 tables, with their own hashed manifest.

Every moment and UV table carries the witness hash. Manifests record table
hashes, precision, polynomial degrees and absolute source norms. The simple
coefficient-sum norms here are conservative; they need not equal optimized
source-norm estimates used by the tail agent.

## Reproduction

Run from the workspace root using the bundled Python executable. The scripts
insert `tmp/certification_deps` to load python-flint0.9.0.

```powershell
$py='C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$base='research/all_direction_certification_20261007/head'
$out='output/all_direction_certification_20261007/head/g2'
$w='output/full_direction_20261007/right/g2_guarded192_refined.npz'
& $py "$base/build_moments.py" --witness $w --out "$out/moments128" --N 128
& $py "$base/threshold_batch.py" --witness $w --moments "$out/moments128" --out "$out/threshold_caps.json"
& $py "$base/head_centered_batch.py" --witness $w --moments "$out/moments128" --out "$out/centered" --numerical-costs output/full_direction_20261007/right/head384_candidates/g2_guarded192_refined.json --start 2 --stop 6 --relative-slack .01 --absolute-slack .0001 --max-cells 6000
```

Use a **fresh output directory** for new trials; batch drivers skip existing
per-spin records. `numerical-costs` supplies only requested stopping targets;
it is not used to prove any cell or sign. For the original head0 replay:

```powershell
& $py "$base/replay_head.py" --witness $w --moments "$out/moments" --record "$out/head_l0.json" --ceiling .9 --out "$out/head_l0_replay.json"
```

The replay checks witness and dependency hashes, exact rational partition
coverage, all singular gaps, and every cell. It does not rerun optimization
and it does not accept the stored cell upper bounds on trust. A separate
adapter is required for centered/Taylor records so it instantiates the same
enclosure class and checks the same gap definition.

## Work still required

The competitive220 claim remains incomplete until rigorous all-energy head
components cover the chosen finite-spin cutoff and a matching analytic tail
covers all larger spins. Neither the numerical384-node head nor sampled
high-spin signs fill that gap. Weak unrelated all-angle bounds do not certify
the numerical1.0452 candidate. New robust220 optimization is being pursued in
parallel to trade a small objective increase for easier threshold/tail proofs.

Earlier stopped head4 trials did not produce completed component records and
must not be quoted as certificates. Their console upper bounds only described
an unfinished adaptive run. All frozen prior numerical and certificate files
remain unchanged.
