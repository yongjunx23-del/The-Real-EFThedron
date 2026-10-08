# Independent source review of the upper83 banded assembly

Reviewed on 2026-10-07 by the lower-bound sub-agent, independently of the upper-bound assembler author. This is an internal source and exact-accounting review, not an external mathematical audit or a proof-assistant verification.

Reviewed source: `research/peak_certification_20261007/assemble_upper83_banded.py`.

SHA-256: `38bc5bd8db046bf8188d63626c27669f0ce9b7f53e307f30a657d3d2be2cf0f4`.

Result: **PASS for the source-level component substitution and its checks.** This review does not assert a completed directional bound before the global certificate and this assembler both finish successfully.

The wrapper preserves the complete same-witness global certificate and substitutes only its two finite-spin tail contributions for spins 42 through 448. It verifies the old exact five-component sum, nonnegativity, witness identity, tail sum, and current hashes before substitution. The head and both infinite-spin contributions remain the exact values of the original complete packet.

The replacement requires the independently reconstructed coverage ledger in `upper83/tail/independent_coverage_accounting.json`, whose source `review_upper_tail.py` was also inspected. That ledger independently indexes every touching closed UV interval, rather than importing the banded producer's coverage helper. It checks 112,962 producer/replay UV cells and 4,096 real-norm cells for all 204 finite spins. It binds the exact replacement file and real-norm replay file. The reviewed ledger's replacement cost is approximately `0.00004384769205095601`.

The wrapper verifies the full real-cell partition from zero to one; matching endpoints, regions, and spin counts; the maximum of producer and replay real-norm bounds; positive local coercivity constants; outward reciprocal weights at least `1/(2*kappa)`; and every exact weighted contribution. It checks regional and total sums and requires equality with the independent ledger. It then removes precisely the old finite compact and finite UV/sliver contributions, verifies the new finite cost is smaller, and recomputes the complete support by exact rational arithmetic.

The direction normalization uses an exact rational lower enclosure of the Euclidean norm, checked by integer-square-root inequalities. Because the total support is nonnegative, division by this lower norm and upward decimal rounding is conservative. The final packet binds this review, the wrapper source, the global assembler, the original complete packet, and the union of their proof dependencies. The original complete packet is preserved and the new output refuses overwrite.

No mathematical or schema defect was found in this source revision. The replacement's physical-kernel correctness continues to depend on the established compact/UV coercivity and real-norm implementations and the assumptions stated in the complete global certificate; the independent ledger checks coverage and accounting, not independent rederivation of those physical kernels.
