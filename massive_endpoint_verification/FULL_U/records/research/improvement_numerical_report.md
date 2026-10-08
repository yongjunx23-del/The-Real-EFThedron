# Upper-g0 numerical campaign, 2026-10-06

The frozen new candidate is `output/improvement_20261006/right_roy18j8_strip7j4.npz`, SHA-256 `032ed315425850abd5ac75557795d959141672ba26e2712f2b46bb42a61c03c8`. It contains only names, multipliers, direction=(1,0). It is a numerically improved witness, not a full exclusion certificate.

It has 266 physical rows: H2..40 (39);33 rational cosine rows and18 rational sine rows;X4;95 Roy rows with B=11.9,J=0,2,4,6,8,n=0..18;80 tapered strips over [12,24] and[24,48],T=4/3,j=0..4,n=0..7,taperpower1. Roy/rational rows were inherited from warm starts, so the final JSON config describes requested additions, not the complete inventory. Row names are authoritative. No fractional or odd-Roy rows occur.

## Preferred measured values

- Fine head through every even spin40: **2.6989434549793883**, `right_final_peak_n32t160_audit.json`.
- Coarser head:2.698942932248419, `right_final_peak_n24t80_audit.json`; difference5.22730969e-7.
- Added every even spin42..100,4.02<mu<49:1.4296527419287956e-7, `right_final_source_tail_n48.json`. Coarser result1.4296534642440466e-7. All sampled V/p are negative; maximum−399.93531541548236.
- Added every even spin42..100,entire mu>49:7.817178136767922e-8, `right_final_uv_tail_n128.json`. Coarser result7.817178136751019e-8. All sampled V/p are negative; maximum−15.756258108617779.
- Combined measured support: **2.698943676116444**. Improvement relative to previous measured2.705761189389571 is approximately0.00681751327.

The earlier restricted tau=ell²/mu∈[1,30] result,7.704989053299316e-8 (`right_final_tail.json`), is contained in the full UV integral above. **Do not add it again.**

Unmeasured domains: mu∈(4,4.02) at spins42..100 and all energies at spins>100. No rigorous quadrature error is supplied for any numerical head/tail integral. The spin26..40 bridge first threshold panel uses only20 nodes, and this is still an unproved endpoint approximation. Consequently this total is a partial numerical dual cost, not a physical upper bound.

## Exact grids and quadrature

`improvement_peak_audit.py` splits velocity beta=sqrt(1−4/mu) at the old source and fixed edges0,.12,.35,sqrt(1−4/11.9),sqrt(1−4/12),sqrt(1−4/24),sqrt(1−4/48),.95,.995,1, plus rational peak centers tan(acos(c)/4) and center±2width, where width=(1+center²)(−log r)/4, for r>=.94 and |c|<1. New edges are used only in(.12,.995). The final inventory has56 panels; exact edge values are stored in both final peak-audit JSON files.

Fine replay:32 Gauss nodes per later panel, first linear beta panel160 for spins0..24; first panel20 for spins26..40. Later panels use beta=a+(b−a)sin²(pi(x+1)/4). This is1920 head energy nodes and1780 bridge nodes. Source-energy and angular quadrature orders are both96.

Coarser replay:24 later-panel nodes, first head80, first bridge20; source80×80. This is1400 head nodes and1340 bridge nodes. All threshold checks match angle_audit.py: nonzero spins<=20 atdelta1e-12 with zero sign tolerance, and all nonzero spins<=40 atdelta1e-6; all pass, as does the integrated large-p sign flag.

Ordinary9-panel audits were also run: n80/source64×56 gives2.698822801827036 and n160/source80×80 gives2.698938021967555. These retain rational-peak quadrature bias and are diagnostics, not the preferred estimates.

Source-window tail:8 velocity panels with physical edges4.02,4.1,6,8,11.9,12,24,48,49,using sine-squared maps. Coarse24 nodes/source64×64; fine48 nodes/source80×80. UV tail:7 linear inverse-mu panels with x edges0,1e-6,1e-4,1e-3,1/200,1/100,1/60,1/49. Coarse64 nodes/source64×64; fine128 nodes/source80×80. These map the entire mu>49 ray without a finite UV cutoff.

## Search, conditioning, and guards

The final search used a peak-aware16-node grid, first head panel80,13 even spins0..24, and a zero-cost delta1e-10 threshold probe:12493 physical cells. Source quadrature64×56.414 artificial zero-cost guards restrict325 Roy source samples, one X4 sign, and88 added-spin points (26..40 at11 energies4+1e-6,4.01,4.1,5,8,12.1,20,50,100,1000,1e5). They are numerical search restrictions, not physical null rows.

`improvement_dual_newton.py` is an isolated copy of the previous optimizer with an accepted-iterate checkpoint callback. Core solver files were not edited. `improvement_right_search.py` subtracts each rational row's known H2..H40 expansion during optimization using an invertible triangular row transform, then saves the coefficients in the original physical basis. Final transformed and physical objectives differed by5.4623e-14. This cured the strong overlap of rational rows with the harmonic bank.

Final threshold margin eps2e-6,eigenvaluefloor1e-10,ftol1e-7,gtol.0005,maxiter20. It finished at the iteration limit, not optimizer convergence. Training cost2.698815096042; the difference to the independently measured head is largely support from spins26..40 absent in training. Training improvements alone were not used as bounds.

## Checks and partial certificates

`right_final_roy_sign_certificate.json` exactly certifies the dangerous Roy-source polynomial W>0 with binary64 multipliers interpreted as rationals: lower bound3.8707281815178683e-5; degree26×8; denominator-clearing power8;8 Bernstein leaves,depth3. This is only one source-sector sign certificate, not the full all-spin V bound. Necessary crossing sign:lambdaX4=−32.12360385689058. Sampled W minimum.08057567465468707 on151×151 source mesh.

`extended_transfer_validation.json` extends the independent original-kernel and crossing-pole-amplitude checks to transfer j<=4 and untapered energy n<=9 (covering the n7 taper expansion). Maximum normalized kernel discrepancy2.7674e-11; maximum crossing-pole residual3.7085e-15. These are numerical algebra checks.

Conditional real-tail estimates in `right_final_conditional_real_tail.json`:A40≈.00555478892,A100≈5.04262009e-7. With G=4.3695,the Cauchy estimates are.220325 and.00209923 respectively. At L40 the absolute-angular envelope loses substantial cancellation. If a global V<=−10p coercivity proof were supplied,the same approximateA40 would instead giveA40/20≈.000277739. These are **not certified integrals**; the adaptive integration also issued a roundoff warning. Root's separate Arb work supersedes any numerical envelope estimate once proved.

## Rejected attempt and reproducibility

`right_enriched.npz` is rejected: apparent ordinary-grid training2.6702534 rises to2.7693339 on peak-aware replay. Sign checks alone missed narrow rational peaks. See `rejected_naive_rational_enrichment.json`. The unconditioned peak search was stopped; the useful staged warm starts are `right_peak_conditioned.npz` (187 rows) and `right_roy18j8_stopped.npz` (222 rows,17 accepted iterates,stopped after step size became tiny). The latter file is a saved checkpoint, not an optimizer convergence result.

Run the commands below from the workspace with the bundled Python. New output paths prevent overwriting frozen results.

```powershell
$py='C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$w='output/improvement_20261006/right_roy18j8_strip7j4.npz'
& $py research/improvement_peak_audit.py $w --n 32 --threshold-order 160 --kernel-order 96 --angular-order 96 --out output/improvement_replay/head.json
& $py research/improvement_source_tail_audit.py $w --n 48 --kernel-order 80 --angular-order 80 --out output/improvement_replay/source_tail.json
& $py research/improvement_uv_tail_audit.py $w --n 128 --kernel-order 80 --angular-order 80 --out output/improvement_replay/uv_tail.json
& $py research/certify_roy_dangerous_sign.py $w --output output/improvement_replay/roy_sign.json
& $py research/improvement_transfer_audit.py
& $py research/improvement_right_search.py --warm $w --peak-panels --n 32 --kernel-order 80 --angular-order 80 --strip-rational-harmonics --eigen-floor 1e-10 --max-iter 20 --out output/improvement_replay/refined
```

The audit scripts create their output directories; the exact Roy certificate script requires its output parent to exist (the preceding head replay creates it). The extended-transfer test deliberately writes its research-campaign validation JSON.
