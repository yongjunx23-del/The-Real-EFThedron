# Further numerical endpoint refinement, 2026-10-06

The measured partial cost decreases from 2.69894367611644 to **2.69888154040664**, a numerical improvement of **6.21357098e-05**. This is **not a certified upper bound**. The independently certified directional bounds in the parent campaign use a different, older witness.

## Saved witness and measured components

- Witness: `output/directional_certification/numerical/bridge/bp0_stopped.npz`.
- SHA-256: `fab65bffd85740b1d6043124bd1f2a0324547e1c033f98274340ab0b92beb2e4`.
- Direction: `(1,0)`.
- Every even spin 0..40, all energies, fine numerical cost: 2.698881478706302.
- Every even spin 42..100, 4.02<mu<49: 3.447242743556079e-08.
- Every even spin 42..100, all mu>49: 2.722791006942304e-08.
- Sum: 2.69888154040664.

The unmeasured domains are all energies for spins>100 and 4<mu<4.02 for spins42..100. None of these three numerical integrals has a proved quadrature-error bound. The sum therefore remains a partial numerical cost.

## What changed

The physical inventory remains the previous 266 rows: H2..40; 33 rational cosine and18 sine rows; X4;95 Roy rows with J<=8,n<=18,B=11.9;80 tapered-strip rows on [12,24] and[24,48] with j<=4,n<=7,T=4/3. No physical assumption or new null relation was introduced.

First, the old 266-row witness was refined for16 accepted steps on n24 peak panels, source and angular quadrature80. Its spin<=24 training objective decreased from2.6988150684049574 to2.698812924859148. That intermediate is `b0_refined.npz`.

Second, spins26..40 were added to the actual objective, rather than merely checked afterwards. They use the existing overflow-safe20-node first beta panel and16-node later peak panels. Every even spin0..40 is now trained. The same-grid warm cost was2.6989400302853745; the frozen checkpoint after6 accepted steps has transformed training cost2.69888107370747. The run was stopped at this saved accepted iterate and does not claim optimizer convergence.

The source guard and threshold scheme remains: sampled Roy-source positivity, lambda_X4<=0,88 extra spin26..40 guard points, and V/p<=-2e-6 wherever p>1e4 in the training/guard inventory. These are restrictions on the dual search, not extra restrictions on physical amplitudes. Harmonic subtraction is an invertible numerical reparametrization; saved multipliers are in the physical row basis.

## Independent replay

The fine head uses32 nodes per rational-peak panel,160 on the first head panel,20 on the bridge threshold panel, and96x96 source/angular quadrature. Its difference from the coarser training objective is 4.04998832e-07. This is a numerical convergence diagnostic, not an error certificate. All three existing checks passed: threshold through spin20 at delta1e-12, through spin40 at delta1e-6, and the integrated large-p sign guard.

The source-tail audit uses48 nodes and80x80 source/angular order. The UV audit uses128 nodes in each of seven inverse-mu panels, covering the entire mu>49 ray, with80x80 source/angular order. Both have zero sampled positive-V contribution; their maximum sampled V/p values are -410.802652 and -15.5995592, respectively. These sampled signs are not all-spin coercivity certificates.

## Reproduce

From the workspace, using the bundled Python:

```powershell
$py='C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $py research/directional_certification/numerical_direction_search.py --warm output/directional_certification/numerical/b0_refined.npz --b 0 --peak-panels --train-bridge --n 16 --kernel-order 64 --angular-order 64 --strip-rational-harmonics --eigen-floor 1e-10 --ftol 1e-8 --gtol 1e-4 --max-iter 6 --out output/directional_certification/numerical/replay_endpoint
& $py research/directional_certification/numerical_campaign_audit.py output/directional_certification/numerical/bridge/bp0_stopped.npz --outdir output/directional_certification/numerical/replay_fine --n 32 --threshold-order 160 --kernel-order 96 --angular-order 96 --tails
```

The generic direction driver supports normal `(1,b)` and the campaign driver can reuse its kernel cache across directions. No nonzero-angle optimization was completed in this refinement; the parent campaign's new certified angular points are separate results.
