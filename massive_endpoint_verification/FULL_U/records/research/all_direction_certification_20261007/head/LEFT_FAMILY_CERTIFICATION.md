# Exact lower-left objective-shift family

Base witness is the fixed left180 NPZ, SHA256
`e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b`.
For each exact rational b in{1/4,1/2,1,2,4,8}, retain every multiplier and
source row but replace the objective(-1,0) by(-1,-b). This is a derived
witness definition, not a claim that a different NPZ hash supplied the
original source tables.

The real fieldU, corrected coefficientalpha, derivative coefficients,
Roy/strip polynomials and all source moments are unchanged. The imaginary
field satisfiesV_b=V_0-b*p, where
p=P_ell(1+2C/(mu-4))/(mu-C)^3 andC=4/3. Sincep>0 for everymu>4,
H(U,V)=sqrt(U^2+V^2)+V decreases pointwise withb. Thus every existing
threshold/source-gap bound and every completed base-tail cost is also a
valid bound for the shifted functional. If local tail coercivity is used,
kappa_b=kappa_0+b provides an optional further improvement.

For beta^2=1-4/mu, e=1-beta^2, dbar=4beta^2 andQbar=4-Ce, the elementary
shift kernel in the scaled head fields is

```
p/e^2 = e R_ell(dbar,Ce) / (dbar^ell Qbar^3),
R_ell(dbar,Ce)=sum_k binom(ell,k)binom(ell+k,k)dbar^(ell-k)(Ce)^k.
```

`replay_left_family.py` first replays the original saved rational partition
and bounds the six shifted heads on that same full partition. It keeps the
base direction immutable and caches only underlying field/jet enclosures.
The original support endpoint cache is never reused as a shifted support
value; new endpoint supports are computed from reusable endpoint fields.

Each shifted cell takes the best of independently valid bounds:

* The original cell upper bound, by monotonicity inV.
* A shifted interval Riemann support bound.
* Convex endpoint interpolation with the original real-field jet, the
  shifted imaginary jet, and the unchanged uniform source remainder.
* The normalized negative-V estimate usingV_0/p-b where that original
  normalized branch is valid: spin0 only belowbeta.5, higher spins below
  beta.99. No old spin0 source-activity formula is used abovebeta.5.
* A bound from a smaller positiveb, again by monotonicity.

The interval kernels and derivatives are enclosed by Arb; source geometric
remainders are unchanged. All source/threshold gaps are recomputed for the
base objective and their proved upper bounds are reused for everyb. No gap
or infinite-energy endpoint is removed from the partition.

Original replay records retain their existing schema and remain separate.
For spins whose original replays were already frozen, additional original
replays go into`family/base_replays/`; old replay files are not overwritten.
The shifted records are`left/family/head_l{ell}.json`, with exact objectives,
base witness SHA, code/source hashes, gap bounds, per-cell dyadic bounds and
per-spin sums. A full family direction requires every even spin0..40 and a
matching base-tail certificate, or a separately improved shifted-tail bound.
Until those are assembled, these are finite-spin components only.

The point-identity diagnostic`left/family_shift_diagnostic.json` passed18
checks for spin28 atbeta.5,.9,.99 and all six shifts against direct field
evaluation with the objective changed. This is an implementation diagnostic;
the interval bounds and the analytic shift identity provide the proof.

Unit-direction support for a plot is the final unnormalized cost divided by
sqrt(1+b^2), using outward interval rounding. Plot angles are approximately
194.036,206.565,225,243.435,255.964 and262.875 degrees; the exact objects are
the rational directions(-1,-b), not rounded decimal angles.

A frozen ten-slope variant, `replay_left_family10.py`, extends the list with
16,32,64,128. Only the exact slope list and descriptive count differ from
the six-slope implementation. New remaining replays use this variant and
write shifted records in`left/family10/`. Earlier six-slope files remain
untouched and may be supplemented by separate ten-slope runs. No direction
is marked complete until all21 spin components contain that exact slope.
The additional plot angles are approximately266.424,268.210,269.105 and
269.552degrees. Their mathematical directions remain exact(-1,-b).

The resulting fixed family is not automatically competitive. In particular,
the b=4 spin0 contribution alone is6.751779199397851/sqrt(17), approximately
1.6376 at255.964degrees. An optimized numerical witness at255degrees has
support approximately0.99367. These directions are not identical, but the
size of the difference rules out describing the fixed family as a uniform
match to the optimized lower arc. The user requested competitive bounds,
so supplement-only passes to fill the four largest slopes were stopped
before launch. Already-running computations and all existing records are
preserved. Remaining required base replays use the original six-slope
driver, to finish the potentially useful directions closer to180degrees.

`assemble_left_family.py` creates`left/family_packets.json`. Each exact
slope receives a separate status, exact dyadic head sum when complete,
all21 spin checks, original/head/family paths and hashes, source rebuild
checks, and dependency hashes. It verifies per-cell dyadic sums against
the recorded per-spin upper bound and verifies the exact rational union
of gaps and ordinary cells. Missing slopes/spins are explicit and no full
bound is emitted for an incomplete packet. Completed packets still need
the same-base all-spin tail and outward unit-direction normalization.

For some later higher-spin replays,`replay_left_family_fast.py` replaces
the old source-polynomial evaluation with the algebraically reviewed
`FiniteFast` contraction. It otherwise copies the six-slope driver exactly,
requires spin>6, and records both its own source hash and`finite_fast.py`.
Every component identifies the selected implementation through those
hashes. The252 interval/jet overlap diagnostics preceding adoption are in
`left/fast_interval_diagnostic.json`; they are implementation checks only.
