# Original 150-row search settings

These settings explain how the base witness was found. They are not prerequisites for verifying the saved witness, and none of the sampled guards below replaces the later interval proof.

`output/angles_resummed/angle_000p000.json` records the configuration: `n=64`, `lmax=24`, `hmax=40`, `roy_order=14`, `roy_spin=6`, `transfer_order=5`, `transfer_taper=1`, `cross_degree=4`, radii `[.8,.95]`, `fractional=false`, objective `(1,0)`. The saved status is `iteration_limit` after60 accepted iterations. Its finite-grid value is2.7057636978806574, not the final certified upper2.706617.

## All original search restrictions

The original `angle_search.py`/`search_problem.py` path uses only the common high-p sign guard in `dual_newton.py`: at every integrated or zero-weight probe column with `p>1e4`, require `V/p<=-2e-6`. The fixed record reports1891 such columns. The sign guard is an inequality in the dual multipliers, not an extra condition on physical amplitudes.

There is one zero-weight threshold point `mu=4+1e-10` at every even head spin0..24. Its weight is zero in the objective but it can enter the high-p sign mask. No additional Roy-source sample columns, independent X4-sign column, or spin26..40 columns are appended by this original driver. Those414 extra guards belong to the later enriched numerical search. Exact Roy-source positivity and all-spin crossing/coercivity properties were verified afterward for the frozen150-row coefficients.

There are no explicit amplitude-derivative tests, imposed elastic equalities, global g2 priors or primal-amplitude truncation constraints in this dual search.

## Grid and source integration

The spectral grid uses beta=sqrt(1-4/mu). Its sorted edge set is

```
0, .12, .35, sqrt(1-4/11.9), sqrt(1-4/12),
sqrt(1-4/24), .95, sqrt(1-4/48), .995, 1
```

The expression is sorted numerically; this gives9 panels. Use64 ordinary linearly mapped Gauss-Legendre nodes in every panel. This original grid does not use the later sine-squared clustering or the56-panel rational-peak grid. The physical measure is `4*(2ell+1)/pi * d beta/(1-beta^2)^2`. Each of13 spins has576 positive-weight nodes and one zero-weight threshold probe, for7501 total columns.

The retained source defaults are Roy48 energy by40 angular Gauss nodes; untapered strip40 energy by32 transfer nodes, with taper power1 assembled by exact finite Legendre-basis algebra (evaluated numerically in the search). Real PV residues and removable crossed poles are handled analytically by the kernel routines. The interval proof later replaces these source quadratures with exact polynomial moments and explicit series remainders.

## Solver settings and replay

`angle_search.run()` calls `dual_newton.solve_constrained_dual()` with eps2e-6, ftol3e-11 andgtol1e-6. The retained solver defaults add threshold1e4, trust_radius10, eigen_floor1e-8. At each step the equilibrated-Hessian eigenvalue floor is `max(max_eigenvalue*1e-8,1e-12)`. Its CLARABEL QP has max_iter250, gap tolerances1e-10, feasibility tolerance1e-11, and no explicit max_threads argument. It uses a .995 feasible step fraction, Armijo1e-4, at most70 halvings, acceptance residual<=max(eps*.01,1e-13)=2e-8, and requires both function decrease and Newton descent tolerances for convergence. A violated start is projected to an interior feasible point; projection failure is reported. These are computational choices, not physical constraints.

The exact archived initial warm file for this particular150-row run is not identified by its compact JSON configuration. Do not invent it or claim a bit-identical optimization history from zero. The exact final coefficient table is sufficient to replicate its certificate. The following starts a new search with the same finite problem, using the frozen witness as its warm start, and writes elsewhere:

```powershell
& $Python research/angle_search.py --angles 0 --n 64 --lmax 24 --hmax 40 --roy-order 14 --roy-spin 6 --transfer-order 5 --cross-degree 4 --radii .8 .95 --max-iter 60 --warm output/angles_resummed/angle_000p000.npz --out output/replication_original_search
```

It may improve or slightly change the coefficients; such a new witness requires its own full certification. To reproduce the reported bound exactly, keep the supplied150-row witness unchanged and use the certification guide.

Sources: `research/angle_search.py`, `research/search_problem.py`, `research/massive_search.py:spectral_grid`, `research/dual_newton.py`, `research/roy_kernels.py:build_roy_bank`, `research/transfer_kernels.py:build_transfer_bank`, and the saved angle JSON. Current code hashes are included in the portable bundle manifest.
