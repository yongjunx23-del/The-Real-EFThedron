# Internal review: 265-degree spin-2 threshold cap

Reviewed `certify_head_constant.py` and the paired finite-width preflight
record for witness SHA256
`5aae4902c2db47f5e3a5f5a37d87c136cd7ec3cb191752d8f3b2f7ee04d6d8c3`.

The wrapper overrides only `threshold_cap` and delegates it to the unchanged,
previously reviewed `FiniteSpin.threshold_cap`. The interior is still the
unchanged `FiniteCentered` evaluator. No multiplier, constraint, or kernel is
changed. The earlier quadratic-cap failure is a limitation of that enclosure
method, not evidence that the witness has a positive threshold singularity.

The selected cap encloses `V/p` by a strictly negative interval over the
entire finite interval `0<=beta<=1/256`. It also encloses `beta*U`, retaining
the possible derivative-row `1/beta` behavior. The unchanged formula uses
`sqrt(U^2+V^2)+V <= U^2/(-2V)` and integrates the remaining explicit positive
power `beta^(2*ell-2)` analytically. For ell=2 this is integrable. All source
geometric and logarithmic remainders remain included.

At 256 and 384 bits, the cap is strictly negative and the resulting support
cost is about `1.10e-7`. The full partition must use exactly this endpoint,
include the other source endpoint caps, and be replayed through the wrapper
so the `constant_centered` model is registered. The wrapper adds its own hash
to the unchanged driver's dependency manifest.

No mathematical defect was identified in this composition of existing
methods. This preflight is not a full finite-spin integral or an all-spin
certificate. It is an internal review, not external or formal verification.
