# Internal review: 261.3814-degree compact adaptation

Reviewed on 2026-10-07 against the frozen 259.294118509649-degree versions.
The finite evaluator changes only the witness path and hash, exact objective,
source-cache path, and maximum finite spin from 74 to 70. The objective agrees
with the exact binary64 rationals in the NPZ. Source division, geometric
remainders, removable poles, derivative terms, and the combined forward and
crossing polynomial are unchanged. The analytic driver changes only witness
imports and output/dependency paths; its interval list and formulas are identical.

The freshly contracted source tables must have jmax=skmax=70 and N=256,
covering all finite spins used. Their exact-beta/shared-table audit and the
separate 512/640-bit finite covers and 400/512-bit analytic covers remain required.
A smaller finite cutoff is valid only where the saved analytic envelope starts
at the immediately following even spin. The final assembler must check that
union on every energy rectangle. This review does not assert completion of
those jobs or certify a complete coupling bound. No formula defect was found
in the reviewed adaptation. This is an internal review, not formal verification.
