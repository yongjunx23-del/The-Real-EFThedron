# Fresh optimized 269.5-degree witness: UV proof

The witness is `output/weak_regions_20261007/lower/search/angle269p5.npz`,
SHA256 `f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae`.
The exact objective is
`(-2515252609197481/288230376151711744,
-9006856288192973/9007199254740992)`.
All 220 stored binary64 coefficients are interpreted as exact dyadics.
This witness is different from the monotone shifted 268-degree witness.

These new wrappers preserve the mathematical kernels of the completed
268-degree campaign. The changes are the input/output paths, witness
SHA256, exact objective, and a stricter proof target kappa=1/10 instead
of 1/500. The finite scheduler runs four independent groups in parallel.
The adaptation generator records the original source scripts; the final
component packets bind the actual new script contents.

The source cache uses exact rationals with N=12 and K=80. It is rebuilt
independently with the same builder, entry by entry. The algebra audit
also specializes the uniform spin polynomial to spins 450, 452, and 600,
compares exact source/crossing polynomials, and checks 18 field enclosures.
These tests do not substitute for the full interval cover.

Every even spin 42 through 448 has an adaptive rational tau partition,
where tau=ell(ell+1)/(mu-4), through min(256,ell(ell+1)/45).
Production uses 400-bit Arb arithmetic. Replay recomputes exactly the
same partition at 480 bits. Every cell must have a finite directed
upper for V/p strictly below -1/10.

For every even spin at least 450, eta=1/[ell(ell+1)] is covered by the
closed interval [0,1/(450*451)]. The uniform polynomial envelope is
certified on the whole eta interval and an adaptive tau partition of
[0,256], again at 400/480 bits with kappa=1/10. Eta=0 is included.

The separate analytic large-tau bound covers tau>=256 for every even
spin at least 108, with mu>=49. Exact crossing Bernstein coefficients,
all negative shifted powers, and Jensen's inequality prove kappa=1.
The analytic bound is evaluated at 320 and 400 bits. Spins below 108
cannot have tau>=256 when mu>=49. Thus there is no missing region.

The assembler checks all exact partitions, paired precision records,
directions, source hashes, spin/energy coverage, and the large-tau seam.
The weaker of the two precision margins is retained locally. The
global finite margin is the minimum local margin and 1; the infinite
global margin is the explicit 1/10.

`run_fresh269p5_components.py`, one directory above, preserves the exact
commands, logs, thread environment, and elapsed times for the source,
UV, independent field, and real-norm jobs. The UV assembly command is:

```text
python research/weak_regions_20261007/tail269p5/assemble_uv269p5.py --angle 269p5
```

A subsequent `--verify-only` invocation checks the frozen result without rewriting it.
The standard all-energy real-norm calculations and the separate compact
coercivity proof must also pass for this same witness before a complete
coupling bound can be claimed. Finite heads are separate as well.

This is an internal computer-assisted continuum proof, conditional on
the physical assumptions and reviewed algebra. It is not a proof-assistant
formalization or an independent external implementation.
