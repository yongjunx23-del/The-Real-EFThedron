# Independent review of the85.8424675-degree positive-mixture adaptation

Agent `weak_upper_search` compared `extra_upper/mixture_head.py` and `assemble_extra_upper.py` against the previously completed, corrected86-degree mixture implementation.

The producer changes exactly the target-direction file and the applicability guard of its normalized estimate. The new guard permits Centered spin0 only below beta1/2 and Centered spin2 only when `beta_upper^2 <= 79/119`, equivalent to `mu_upper <= 119/10`. This is the Roy real-source activity domain; the strip real sources are inactive there. Other cells retain the actual-field, convex-endpoint, or positive-combination fallback enclosures. No physical coefficient or source kernel is modified.

The assembler changes the target/output locations, reads the directly valid producer/replay records, and removes the obsolete historical-repair ledger. In its place it independently checks the safe applicability guard on every selected normalized cell. The old correction machinery is therefore unnecessary for this new calculation, because the unsupported branch is never selected.

The exact positive weights, all220 exact rational mixed multipliers and objective identity, complete base-certificate verification, and explicit source-table hash bindings remain unchanged. The two base head decompositions are reconstructed exactly before subtracting retained spin0/2 costs. The new complete bound retains the larger of the two direct whole-spin enclosures, optionally capped by the old valid convex contribution. Source gaps use the identical exact geometry and outward base costs. The entire ordinary partition plus gaps is checked exactly. The direction norm uses an exact integer-square-root lower bound and the displayed support is rounded outward.

Review result: the inspected adaptation preserves the mathematical argument and closes the previously identified Centered-source applicability problem at the producer and assembly stages. This review concerns source logic; the final packet additionally requires all four complete head runs, successful exact assembly, and dependency hash verification. It is internal review with shared mathematical implementations, not external or proof-assistant verification.

The completed assembler binds this review, its own final source and all actual head-producer sources in its dependency manifest. Earlier sealed proof files remain unchanged.
