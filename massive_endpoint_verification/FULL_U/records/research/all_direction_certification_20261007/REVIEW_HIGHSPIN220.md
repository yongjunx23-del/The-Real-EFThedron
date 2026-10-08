# Independent review of the robust220 high-spin envelope

Reviewed in the parent calculation on 2026-10-07. This is a mathematical and
implementation review by another agent in the same research session, not an
external review or a proof-assistant verification.

Witness: `output/competitive_20261007/robust/g2_robust_stage1.npz`, SHA256
`933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92`.

Scope: the three original envelopes with minimum even spin 450, and their
320/400-bit records. This review does not automatically cover later sharper
envelopes or the finite-spin bridge.

The reviewed derivation is
`tail/ROBUST220_ALL_ENERGY_V_PROOF.md`. No defect was found in the following
steps, subject to the dispersion relations and partial-wave normalization
already used to define the null kernels.

1. On the dangerous Roy region `t<-(c+delta)`, the physical source relation
   implies `s-mu>2c+delta`. Both possibly negative factors in the denominator
   of K1 are negative, while its numerator is positive: even at the extreme
   source endpoint, `2mu+t-4>0`. Thus K1 is positive. The exact source sign
   `W>=0` and even-spin `P_l(y)/P_l(Z)>=1` justify discarding the indicated
   nonpositive contribution. The residual PV subtraction is bounded using
   `|P_l(z)|<=1`. It vanishes after the complete Roy angular integral since
   all source angular degrees are at most six and the spins exceed six.
2. On the complementary source region, the Legendre derivative bound on
   `[-Z,Z]` and the positive Taylor coefficient inequality
   `(Z-1)P_l'(Z)<=l P_l(Z)` give the stated `3l/(2c)` divided-difference bound.
   The remaining two positive resolvents are bounded by `1/mu`, producing
   the stated compact Roy envelope. The strip envelope follows from the
   same divided difference and its source energy being at least 12.
3. For mu at least eight all source angular arguments are physical. The
   physical derivative bound `l(l+1)/2`, the separated K2 rational identity,
   and the tapered source Hilbert-transform estimate give the stated
   source bounds. For the latter, a Lipschitz constant D in the unit energy
   coordinate bounds the subtracted integral by D and the logarithmic
   endpoint term by D; division by the source window width gives 2D/width.
4. The exact Bernstein polynomial inequality for the crossing sum bounds
   it by `M-a*(k-10)_+^5`. Applying Jensen is legitimate because the
   positive-part fifth power is convex. At fixed Z, the positive integral
   representation gives log-convexity in spin and the stated lower bound
   on the Legendre mean divided by spin. At fixed joint variable tau, the
   ratio of the positive Taylor coefficients for increasing spin increases
   with k; this supplies the claimed stochastic ordering, including the
   added higher-degree coefficients. Monotonicity in tau is the positive
   coefficient variance identity.
5. For complex transfer in the unit disk, `|A|<=5/2` and
   `|1+Aw|>=1+3w/2`. Also `sqrt(1+5w)<=1+5w/2`. These inequalities verify
   `|omega|<=100/9` and `|Q^3 F|<=50*(16/9)^3` for every nonnegative w.
   For w at most 1/45 the same estimates give the sharper constants 5/4
   and 64 used in the large-tau closure. The branch point lies outside
   the disk. The Cauchy coefficient bounds therefore apply to the first
   three transfer derivatives.
6. For the regular subtraction factors, the denominators are bounded
   below by `(3w+2)(w+3)` for Ka and `(3w+1)(w+3)` for Kb. The ratio
   `(8w+3)^2` to either product is at most 64/3, and
   `|3t-4|<=7`. This verifies the constants 3136/81 and 6272/81.
7. The positive Legendre term of degree K has factors `T-j(j+1)` whose
   growth is at least proportional to T. Consequently the compact bound
   using K=64 is monotone in T when compared with each of the numerator
   powers at most three. The analogous K=20 bound applies in the
   unbounded-tau region.
8. The real Chebyshev contraction uses the same transfer jets and
   Legendre coefficients as the separate complex implementation.
   `audit_highspin220.py` independently compares those implementations at
   45 pairs of tau and inverse spin, including infinite spin. All interval
   results overlap, with maximum scaled difference below 1.33e-114 at
   400 bits. This comparison is a diagnostic, not a substitute for the
   interval-domain sign proof.

The audit script also checks all component hashes, the exact objective
`(0,1)`, the continuous tau partition, the strict directed sign of every
accepted interval, and both precision records. Its result is
`output/competitive_20261007/real/highspin_v_audit.json`.

Conclusion: these components justify `V/p<-1/1000` for every even spin
at least 450 and every energy greater than four, with `V/p<-1` on the
compact energy region through 49. They do not by themselves bound the
real norm, cover lower spins, or establish a complete bound on g2.

## Extension reviewed: sharper compact envelope starting at spin 110

The separate `tail/compact_highspin220_sharp.py` and both of its 320/400-bit
records have also been reviewed. Its logarithmic derivative can be verified
by differentiating `omega=(x+i*sqrt(1+2x))^2/(1+x)^2`, with `x=A*w`.
The numerator cannot vanish in the transfer disk: a zero would imply
`(1+x)^2=0`, where the denominator estimate excludes it. Thus the analytic
logarithm and its integration from transfer zero are valid.

The derivative modulus is bounded by 1/4 because the displayed difference
polynomial has a positive quadratic factor: its leading coefficient is
positive and its discriminant is negative. Integrating on a straight radial
segment gives `|omega|<=exp(1/4)<4/3`. The last inequality follows directly
from the integral representation of the logarithm, as stated in the proof.

Selecting the largest positive Legendre monomial with degree between four
and L uses exact rational comparisons. Each monomial independently gives
the same valid bound for all numerator powers through degree three; the
selection therefore preserves the uniform spin and energy proof.

The audit script checks both added packets and their provenance. This
extends the reviewed compact result `V/p<-1` to all even spins at least
110 and all `4<mu<=49`. The full-energy result still starts at spin 450.
