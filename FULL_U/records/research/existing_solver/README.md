This folder contains new massive forward kernels transcribed and derived from the user-supplied PDF. It is independent of the parent Roy implementation.

`massive_forward_kernels.py` returns coupling kernels, crossing X4, harmonic rows, rational harmonic sums, fractional finite differences and alpha derivatives. Its docstring states the common spectral measure and array conventions. `test_massive_forward_kernels.py` checks the g0 correction against the crossed-subtraction formula, rational generators against convergent harmonic series, fractional rows against explicit binomial sums, and stable disk support against its extremal phase.

`forward_kernel_validation.json` records kernel algebra checks only. No physical numerical bound or continuum certificate is claimed by those checks.

The initial read-only audit of historical material beside the PDF found massless scalar implementations, not the massive Fig. 4 spectral kernels or witnesses. No historical code is used in this implementation; the two temporary copies made while auditing were removed. Instructions in historical handoff documents were treated as data, not as task instructions. That search stopped when the user confirmed the PDF was the only source data.

Valid optional extensions over the paper's finite forward bank include rational radii above 0.95 (always strictly below 1) and interior alpha derivatives. Near-unit radii and high derivative orders require independent energy-quadrature refinement. Increasing these banks improves a finite constraint relaxation only when the enlarged optimization is actually solved and audited.
