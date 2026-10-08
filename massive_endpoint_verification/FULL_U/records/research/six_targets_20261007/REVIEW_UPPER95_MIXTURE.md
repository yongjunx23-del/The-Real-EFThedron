# Internal peer review of the upper95 exact mixture

Reviewed on 2026-10-07 by the lower-target agent, independently of the upper-target author. Result: **PASS for the reviewed source construction, domain handling, and exact assembly logic.** A completed directional bound still requires every selected producer/replay record and successful final assembly. This is not external mathematical review or formal verification.

## Reviewed source hashes

| Source | SHA-256 |
|---|---|
| `upper/upper95_common.py` | `1ad3f3962a7282f57e8aa435147c9395f31d082b4fda3119937eadfd57c2f787` |
| `upper/upper95_head.py` | `a9fe177217c984f4fa4d43c6dde0ee95b70c8b07ecae2bfad3b4cc16c11555b7` |
| `upper/assemble_upper95_mixture.py` | `8a49fa59ab7673fab78c377d510d5137ebc9fb1f31439a59bbd5aef5e987ef61` |
| `upper/assemble_upper55.py`, imported accounting helpers | `1bd478bb0c45fc3ea4d8f4a4c29f56bbea7922f5d2584398087372e6076c0c5e` |

Paths in this table are relative to `research/six_targets_20261007`. The reviewer also read the original `DirectionalZero`, `DirectionalLow`, `SpinZero`, `LowSpin`, `Finite220`, and `FiniteCentered` implementations used by these sources.

## Field conventions and exact witness

The mixture combines the completed 220-row upper83 witness with the completed legacy b6 witness. The second contribution consists of the exact saved 150-row witness, objective `(1,6)`, and the additional crossing coefficient `107/5` multiplying X4. The original Directional classes add exactly `6*p+(107/5)*X4` to V, with no change in U; their derivative and normalized versions use the corresponding same addition.

Both field implementations use the same beta-scaled U,V conventions and the outer factor `4*(2*ell+1)/pi`. The normalized representation returns the same scaled U, `Vn=e^2*V_scaled/p`, `pinv=1/p`, and `e=1-beta^2`; hence `V_scaled*pinv*e^2` agrees with Vn. The positive rational mixture weights are solved and checked exactly against the objective vectors. The exported hierarchical coefficient definition preserves both bases, every exact base coefficient, both exact weights, and the additional X4 coefficient. It must not be interpreted as two identically named coefficient arrays without retaining their distinct source-basis meanings.

The source cache check reconstructs the legacy N64 cache and checks all nested exact fields without replacing the original source bytes. The completed parent certificates bind the underlying source and moment files.

## Unequal gaps and interval validity

The common gap domain is the union of both bases' analytic endpoint/source gaps. Each base contributes its complete old gap bound plus the whole saved nonnegative ordinary-cell bound for every cell whose interior intersects an enlarged gap. This can overcount, which is conservative. The complement is obtained by refining both endpoint lists and checked as an exact partition of beta in `[0,1]` together with the merged gaps. No width-proportional scaling of a cell-integral bound is used.

Ordinary cells use the minimum of valid interval Riemann, convex-endpoint/Taylor, saved whole-cell convex fallback, and applicable negative-V candidates. The spin2 normalized shortcut requires `beta_hi^2 <= 79/119`, equivalent to `mu <= 119/10`. Both Roy sources have exact domain `4..119/10`, so their normalized real-source formula is active throughout this permitted range. Spin0 does not use the normalized shortcut. Other cells use source-aware actual fields and series.

## Assembly accounting

The assembler reconstructs all retained upper83 head components as maxima of their producer/replay values. For b6 it reconstructs the paired 280/340-bit contributions at spins 0,2,4,6 and authenticates the complete 8..40 block, verifying equality with the exact finite-head component of its completed parent certificate. Thus removing the retained spin0/2 costs from the positive weighted complete bound leaves a nonnegative rigorous remainder containing every other spin and both infinite domains.

The new spin0/2 records are checked at 256 and 384 bits, with matching exact objectives, weights, coefficient definitions, complete cell/gap geometry, enlarged-gap ledgers, all sums, and allowed normalized domains. The current assembler also independently reconstructs every selected `convex_base_cells` bound. It retains the larger direct paired bound, capped by the original weighted whole-spin bound. Exact rational arithmetic then substitutes only the two selected spin contributions. Normalization uses an integer-square-root-verified lower norm enclosure and upward decimal rounding. All selected commands must complete successfully, and dependency hashes are bound recursively.

## Independent supplemental checks

Record: `output/six_targets_20261007/upper95/peer_field_and_gap_checks.json`.

SHA-256: `198a0f50784f8f7eccaf5a878b2254008563901ae072e928af8aaa812751af76`.

The separate peer script independently reconstructs the ordinary complement by midpoint classification and checks every enlarged-gap whole-cell charge by strict interval intersection. It also performs 18 interval-overlap diagnostics on the legacy directional additions and normalized field units, including spin2 points through beta=0.81, and checks both exact Roy source endpoints. These are supplemental diagnostics supporting the source review, not sampled substitutes for the continuum producer/replay proof. Paired calculations share mathematical kernels and therefore do not constitute independent implementations.
