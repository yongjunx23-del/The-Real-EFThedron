# Deriving the three endpoint inequalities

Use the coupling normalization and physical assumptions stated in
`FULL_U/appendix.md`. The external-energy variables in the null equations
are distinct from the spectral integration variable, which covers all energies
above threshold. No spectral-energy or spin cutoff is silently dropped.

For each direction `(n0,n2)`, add the specified null equations with their frozen
rational multipliers to the dispersion representation of `n0*g0+n2*g2`.
The result has a real kernel U and an absorptive kernel V. Full unitarity is
`r²+(b-1)² <= 1`, with `r=Re f` and `b=Im f`. Therefore, pointwise,

    U*r + V*b <= V + sqrt(U²+V²).

Integrate this support function over energy and sum all even spins. Any upper
bound on that integral is a valid inequality on the selected coupling
combination. The multipliers need not be optimal for the inequality to be valid.

## Constraints

C08 uses the 150-row set: 39 forward harmonics, 14 resummed harmonics, one
crossing null, 60 folded Roy averages, and 36 fixed-transfer strip averages.
C20 and C21 use the 220-row set: 39 forward harmonics, 14 resummed harmonics,
seven crossing nulls, 60 folded Roy averages, 64 strip averages and 36 transfer
 derivatives. The precise kernels, index ranges and resummation parameters
are given in `FULL_U/appendix.md`, and each retained coefficient and row name
is in `FULL_U/data/multipliers_exact.json`.

The Roy source energies range from 4 to 11.9. The strip source energies cover
[12,24] and [24,48], with transfer in [-4/3,0]. They are homogeneous analytic
nulls after the subtraction coupling has been eliminated. The symmetric-point
coupling map retains its crossed subtraction term. These are the actual
selected constraint sets, not every possible allowed analytic null.

## All-spin control

For low spins, the code integrates the disk support using directed Arb interval
arithmetic. Finite-width caps cover threshold and source endpoints, and
rational cell partitions cover the whole compactified energy interval.

For higher spins, a proved negative absorptive margin `V <= -kappa*p`, with
`p>0`, gives

    V + sqrt(U²+V²) <= U²/(2*kappa*p).

The code integrates an upper enclosure for that real-kernel norm and sums its
analytic spin remainder. This controls a continuum, rather than checking the
sign only on sampled energies or spins.

C08 integrates spins 0,2,4,6 directly, controls spins 8–40 with separate
cellwise margins, and uses an infinite-spin proof from spin42 onward.
C20 and C21 integrate spins 0–40 directly, use finite-spin estimates for
42–448 and an analytic remainder from450 onward. Compact energy ends at49;
the UV proofs then cover all larger energies. The compactified UV proof
includes zero inverse spin, so its tail is not an extrapolation.

## From components to the displayed numbers

C08 is the sum of the four refined head costs, the original slope-zero dense
finite-block cost, and the fixed all-spin tail cost. No positive mixture is
used for this endpoint.

C20 uses the refined scalar head, the other heads through spin40, the compact
real norm divided by twice its unit coercivity margin, and UV contributions
with finite-spin margin1/100 and infinite-spin margin1/1000.

C21 uses the refined heads through40, the locally weighted finite42–448 norm,
and the compact and UV infinite-spin tails. The selected finite component is
`local_cost2048.json`, not the producer's older default.

The assembler sums exact rational upper endpoints. Rounding the final constants
outwards gives `g0 <= 2.705955`, `g2 <= 1.052351`, and `g0 >= -7.751502`.
See the resulting `C08_assembly.json`, `C20_assembly.json`, `C21_assembly.json`
reports for exact values and component accounting. Verification does not rerun
an optimizer. The README distinguishes saved-certificate checks, the optional
fresh finite-spin calculation, and a complete fresh interval proof.
