# Lower-target mixture source review

Reviewed internally on 2026-10-07 by the parent agent, independently of the adapting agent. This is source review, not external mathematical review or formal verification. A complete result still requires successful recorded computations and the fail-closed assembler.

## Construction

For positive exact rational weights w1,w2, the objective and each of the 220 multipliers are combined with exactly those weights. The full disk support is H(U,V)=sqrt(U²+V²)+V. It is nonnegative, positively homogeneous and subadditive. Consequently every unchanged spin/tail contribution can retain the positive weighted sum of the two complete base costs. For each reintegrated spin, either its independently enclosed combined integral or that original convex bound can be retained, whichever is smaller.

The parent read `lower/mixture_head.py`, `lower/mixture_head_spin4_v2.py`, and `lower/assemble.py`. The physical field and Taylor formulas are the same as the sealed mixture implementation. Spins 0 and 2 use FiniteCentered; spin 4 uses FiniteFast. These classes implement the same physical fields with their respective source-aware interval formulas. CellCache delegates to the selected class. Numerical target values do not enter proof costs.

## Domain and interval checks

Both bases must have identical analytic gap intervals for each selected spin. All ordinary-cell endpoints are merged exactly. The resulting ordinary cells plus gaps must cover beta in [0,1] without holes or overlap. On a subdivided cell the fallback uses the entire containing old cell integral, never a width-proportional fraction. Its possible overcount is conservative because H is nonnegative. Gap contributions are combined with positive weights and upper endpoints.

The Riemann, convex-endpoint/Taylor and negative-V alternatives are upper bounds, so their minimum is valid. The normalized negative-V shortcut is restricted to hi<1/2. This avoids extension of the Centered source-activity formula beyond its supported range. Spin4 v2 catches ValueError/ZeroDivisionError only around an optional Taylor alternative; the already established whole-cell fallback survives. A failed Taylor alternative contributes no claimed bound. Failed preliminary runs are diagnostics, not selected proof records.

## Exact accounting

The assembler authenticates base witness identity, all 220 coefficients, exact objective and positive weights. It reconstructs all 21 retained base head components as max(producer,replay), verifies exact equality with each complete base certificate, and subtracts only selected components from that exact decomposition. It independently checks producer partitions and gap totals, the new 256/384-bit paired partitions, legal methods, negative-V domain restrictions, exact fallback costs, source identities, and component sums. Refinement is explicitly spins [0,2,4] for lower15/lower45 and [0,2] for lower60. Remaining spins and both infinite domains stay in the nonnegative inherited remainder. The unit-normal conversion divides an outward bound by a proved rational lower bound for the norm.

The adaptation is mathematically acceptable under the assumptions and implementation limitations of the complete bases. Final acceptance additionally requires the recorded jobs, moment/source identities and exact complete assembly to pass. Paired runs share an evaluator and are not independent implementations.
