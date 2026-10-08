# Internal review: 259.2941-degree spin-2 constant cap

The same unchanged FiniteSpin.threshold_cap composition reviewed in
REVIEW_CONSTANT_CAP.md is used for witness
`321695c10783e2e0f6e6e28589f132fdb7a68497314b1cf92eb71a07612c704f`.
The finite-width cap is beta in [0,1/256]. At both 256 and 384 bits the actual
directed V/p endpoints are approximately -0.0503143215053 and
-0.0167198889856, strictly negative throughout the cap. The exact dyadic
endpoints are saved in constant_cap_exact_sign.json.

Arb's coarse string [+/-0.0504] discards useful sign information in its
display; it is not the interval endpoint used by the assertion or bound.
The threshold cost uses the actual directed negative upper endpoint and
includes the source remainders and the beta*U bound. Its integral cost is
less than 7.592e-7. The endpoint observer records the unchanged method's
local interval; it does not alter the arithmetic or model.

The earlier quadratic-cap method failed its applicability condition. This
successful constant-negative method does not require a negative quadratic
coefficient. The full ordinary-cell partition still needs its production
and higher-precision replay. No new formula defect was identified; this is
an internal review, not formal or external verification.
