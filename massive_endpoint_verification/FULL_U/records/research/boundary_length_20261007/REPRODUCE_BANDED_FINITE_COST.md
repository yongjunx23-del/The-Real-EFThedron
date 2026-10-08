# Reproduce the energy-dependent finite-spin tail costs

This exact-arithmetic merge uses completed same-witness UV and compact coercivity packets and the saved 256/384-bit real-norm integration records. It performs no new kernel evaluation or integration.

Run from a disposable extracted workspace after the relevant component packets have been verified:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
python research/boundary_length_20261007/banded_finite_cost.py --angle 261p3814 --out output/boundary_length_20261007/real/angle261p3814/finite_banded_cost_recomputed.json
python research/boundary_length_20261007/banded_finite_cost.py --angle 268 --out output/boundary_length_20261007/real/angle268/finite_banded_cost_recomputed.json
```

The alternate output names preserve the frozen cost records. The same commands also accept aliases `85`, `255`, `259`, and `265`. The exact `cost_upper_exact` and every saved cell field should agree with the original; elapsed seconds and the output filename need not agree. The tool uses only Python's standard library. The one exact original workspace prefix, `C:/Users/laure/Documents/ChatGPT/Unitarity codex/`, is mapped to the current extracted workspace; native paths already inside that workspace are retained. Other absolute prefixes are rejected. Dependency bytes are checked against the original hashes.

For a beta cell `[a,b]`, the script checks the compact boundary with exact inequalities against `beta²=45/49`. On its UV part, each even spin from 42 through 448 has

`tau=ell(ell+1)(1-beta²)/(4 beta²)`.

It takes the weakest certified margin over every touching closed tau interval at both precision runs and over any compact or large-tau piece. Both neighbours are included when a cell touches an exact partition seam. The minimum over all 204 spins gives one valid margin `kappa` for the saved aggregate real-norm cell mass `A`.

The support contribution is at most `A/(2*kappa)`. The script rounds the reciprocal upward to 40 significant binary bits, checks the ceiling inequality with exact rational arithmetic, and multiplies and sums the resulting dyadics exactly. Each output cell records the retained real mass, exact chosen margin, limiting spin, conservative reciprocal, exact cost, all touching per-spin UV cell index ranges, and the compact/large-tau coverage flags. No interpolation, quadrature, or sampled sign inference enters this merge.

The result covers finite spins 42 through 448 at all energies. The separately certified head through spin 40 and infinite-spin tail from spin 450 must still be added by the complete-bound assembler.
