# Same-witness ultraviolet certificate for the upper86 candidate

This campaign uses the reviewed upper83 UV mathematics without changing an inequality, truncation order, energy cutoff, spin cutoff, eta-panel count or acceptance margin. Its exact multiplier vector and objective are different. The new input and all actual interval evaluations are rebuilt for that vector; no upper83 numerical interval bound is reused as a bound for upper86.

The fixed input is `output/weak_regions_20261007/upper/angle86_refined.npz`, SHA256 `14e2b30df27e95ede876a8ad75ca190904547fb4ace535c14c6a11725d2df97f`. The objective is exactly `(8816243287893001/144115188075855872,8990329285288823/9007199254740992)`. Every prepared component checks these identities and binds the source kernels and caches by file hashes.

The domains and formulas are those described in the sealed `research/peak_certification_20261007/UV_TAIL_PROOF_AND_REVIEW.md` and `REVIEW_SIGNED_SOURCE2_UPPER83.md`. Finite even spins42..448 cover mu≥49 with tau=l(l+1)/(mu−4), until min(256,l(l+1)/45). The uniform representation covers every even spin≥450 using eta=1/[l(l+1)] in the closed interval0..1/(450*451), partitioned into16 exact equal intervals. Paired400/480bit computations enclose all tau0..256. The eta0 endpoint is included, so taking a large-spin limit alone is never used to substitute for the continuum.

For tau≥256 the paired320/400bit analytic proof applies for every even spin≥108. It keeps every negative shifted crossing power and uses the positive-Legendre mean bound and Jensen's inequality. The exact overlap condition106*107/45<256 covers all smaller finite spins completely with the finite UV partition. The signed source2 expansion is

`−S2(z) ≤ −m0 + sum(k=1..11) |mk| zmax^(2k) + norm2*r^24/(1−r²)`,

where `zmax=3/143`, `r=(317/60)zmax<1`, and `m0=Hs` is checked against the independently contracted joint moments. The constant convention is unchanged: the standalone forward envelope excludesHs, so the signed−S2 term is inserted directly. Every required moment and the absolute-norm remainder remain in the saved packet.

The source cache is regenerated exactly, comparing972 rational moments and9 geometric-remainder coefficients. The finite and uniform source/crossing polynomials are checked for exact identity at spins450,452,600, with18 independently evaluated interval overlaps at512bits. These are algebra checks in addition to the domain certificate, not substitutes for it.

`assemble_uv_peak.py` checks the exact partition endpoints, all same-witness identities, negativity of every retained cell, matched producer/replay partitions, analytic overlap, source dependencies and local margin accounting. The final `--verify-only` pass reads the frozen records without replacing them. UV coercivity alone is not a support bound; the complete directional assembler also requires the compact4<mu≤49 proof, retained-spin integrals and real-source-square norms.

The local witness adapter and signed source review are recorded in `REVIEW_SIGNED_SOURCE2_UPPER86.md`; complete-assembly reviews are separate. These internal reviews and higher-precision replays share kernels. No proof-assistant or external independent verification is asserted.
