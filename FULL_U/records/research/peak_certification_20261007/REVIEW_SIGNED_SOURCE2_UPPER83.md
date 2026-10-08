# Independent formula review of the upper83 signed-source ultraviolet estimate

Reviewed by agent `weak_upper_search` against the frozen `uniform_source_uv.py`, `finite_uv_bridge220.py`, `contracted_forward.py`, and `compact_highspin220.forward_bounds` formulas.

Producer: `research/peak_certification_20261007/tail_upper83/certify_angleupper83_large_tau_moments.py`.
SHA256: `7a0a03110b849011d5c694d4aa0839da56e0f6a8679e135918f73e96c4e58862`.
Witness SHA256: `68e39591f3f9a9d3b5cdfe3ade3924fecb227a329de4a9543729f0607e8f663a`.

For mu≥49, let `w=1/(mu-4)` and `z=w/(1+(8/3)w)=1/(mu-4/3)`. Then `0≤z≤3/143`. The exact signed source2 contraction has the expansion

`S2(z) = sum(k=0..11, m_k*z^(2k)) + remainder`.

The already reviewed source-moment machinery bounds the remainder uniformly by

`abs(remainder) <= norm2*r^24/(1-r^2)`, where `r=(317/60)*(3/143)<1`.

Consequently a valid uniform upper enclosure of the term that enters V/p is

`-S2(z) <= -m_0 + sum(k=1..11, abs(m_k)*(3/143)^(2k)) + norm2*r^24/(1-r^2)`.

Retaining the sign of the exact constant moment improves the previous absolute-norm estimate. All moments and the norm are exact rationals in the hash-bound same-witness prepared source cache. The producer checks the witness identity, `m_0=Hs`, at least12 retained moments, and `r<1`; its prepared source manifest binds `compact_norms.pkl` and the norm/moment constructors.

The cancellation convention was checked explicitly. `ContractedRegion.numerator` includes `-Hs` in its outside term, and the finite/uniform source helpers add `+Hs-S2`. In contrast, `forward_bounds` returns a forward-only absolute outside bound without an Hs contribution. The unbounded-tau producer therefore inserts the signed `-S2` enclosure directly, with no extra Hs term. This matches the complete source term and avoids either double subtraction or an omitted constant.

The remaining proof retains the previously reviewed positive-Legendre ordering and all-powers Jensen/Bernstein crossing coercivity for every even spin at least108 and tau≥256. Its exact objective is included separately. The improved estimate must be combined with the finite-tau proof and other energy/spin components; this review alone is not a complete support bound.

Review result: the new signed-moment enclosure and its source convention are mathematically sound under the same kernel/source assumptions as the frozen proof helpers. Directed interval runs at the two stated precisions and complete coverage/hash checks remain required. This is an internal independent formula review, not external or proof-assistant verification.
