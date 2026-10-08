# Internal review of the 268-degree spin-2 cap refinement

The original quadratic-negative cap passed its sign checks but contributed
more cost than the requested numerical quality target. It remains saved as
an unsuccessful refinement attempt, not a false proof.

The replacement uses the unchanged `FiniteConstantCentered` wrapper and its
existing negative-constant cap. At both 256 and 384 bits the saved direct
dyadic endpoints are strictly negative on beta in [0,1/1024]. This narrower
threshold interval is allowed because the ordinary integration partition
starts at the same endpoint. The final head and its replay must check that
there is neither a gap nor overlap in the full beta coverage.

The new model-selection record, final budget, paired cap screen, and exact
endpoint record are saved separately. They do not replace the original
preflight, selected budget, or stopped spin-2 checkpoint. The final assembler
must use the constant-model spin-2 record and require that its threshold gap
matches the cap record exactly, rather than assume the previous [0,1/256]
threshold size. The mathematical cap and ordinary-cell evaluators are unchanged.

Reviewed on 2026-10-07. This cap review is not completion of the ordinary
integration or its higher-precision replay and is not formal verification.
