# Replication audit: seven massive 4D directional certificates

Audit date: 6 October 2026. This is a new read-only audit of the preserved workflow. No frozen source, witness, cache, or result was edited to prepare it. The focus is the finite block \(\ell=8,10,\ldots,40\), the entire \(\ell\ge42\) tail, and their dependencies and assembly. The separate head guide should supply the detailed derivation and fresh generation of spins \(0,2,4,6\).

## 1. What is being reproduced

The authoritative current result is `output/directional_certification/complete_directional_certificate.json`, SHA-256:

```text
ad71650935d05b15962016bfa85be028856a5b873b299518bd92201ff993497b
```

It contains seven exact rational inequalities, printed below with their outward-rounded constants. The direction is \(g_0+b g_2\), not a unit-normal direction.

| \(b\) | Final outward \(C\) | Crossing multiplier increment \(d\) |
|---:|---:|---:|
| \(-1\) | 2.646050 | 0 |
| \(-1/2\) | 2.676329 | 0 |
| \(-1/4\) | 2.691473 | 0 |
| \(-1/10\) | 2.700559 | 0 |
| 0 | **2.706617** | 0 |
| \(1/10\) | 2.713168 | \(2/5\) |
| \(1/4\) | 2.722994 | 1 |

The trial \(b=1/2,d=2\) is rejected by this certification attempt. A failed sufficient envelope does not establish a physical violation or infeasibility of the true amplitude problem.

The underlying fixed witness is `output/angles_resummed/angle_000p000.npz`, SHA-256:

```text
bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c
```

Its 150 rows are 39 forward harmonics H2–H40; 14 rational cosine rows; X4; 60 Roy rows with \(B=119/10\), \(J=0,2,4,6\), \(n=0,\ldots,14\); and 36 power-one tapered strip rows on \([12,24]\) and \([24,48]\), with \(j=0,1,2\), \(n=0,\ldots,5\). Every saved binary64 multiplier is treated as an **exact rational number**. Decimal row parameters are exact decimals; the legacy label `T1.33333` means exactly \(T=4/3\). Never substitute printed rounded multipliers.

For each direction, \(U\) is unchanged and \(V_{b,d}=V_0+b p+dX_4\), with \(d=0\) for \(b\le0\), \(d=4b\) for \(b>0\). There is no requirement to rerun the optimizer. The optimizer's cost is not the certificate.

The intended amplitude assumptions and the validity of the null identities remain hypotheses. This is an interval-arithmetic, computer-assisted certificate with reviewed analytic estimates, not a proof-assistant formalization or a proof of those physical hypotheses from weaker axioms.

## 2. Environment and preservation rules

Run from the workspace or extracted bundle root. The tested platform is Windows PowerShell, Python 3.12.14. The original interpreter was:

```powershell
$Python = 'C:/Users/laure/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $Python --version
```

On another computer set `$Python` to its Python executable. The core proof uses NumPy, SymPy and python-flint/Arb. The complete archived research environment is pinned by:

```powershell
& $Python -m pip install -r requirements-handoff.txt
& $Python -m pip install --target tmp/certification_deps -r requirements-certification.txt
```

The latter pins `python-flint==0.9.0`; scripts prepend `tmp/certification_deps` to the module search path. The search requirements also include plotting, PDF and optimization packages that are unnecessary for some certificate-only steps. Keep the project-relative directory layout intact. Several saved manifest keys contain Windows backslashes; some readers normalize them, others do not. The full archive has been exercised on Windows, not established as turnkey cross-platform software.

Run Python normally, **never with `-O` or `PYTHONOPTIMIZE` enabled**: assertions implement necessary certificate checks.

**Three workflows must not be confused.**

1. **Frozen verification:** keep all canonical files byte-for-byte unchanged and verify hashes, partitions, exact sums and recorded proof claims. This is fast, but does not rerun every interval kernel.
2. **Fresh mathematical replay:** recompute components with the same code in a disposable copy or write explicitly supported outputs to new paths. Compare their verified bounds and partitions. This need not reproduce JSON bytes: timestamps, timings, absolute paths, library versions and outward ball formatting can differ.
3. **Fresh production of a new complete provenance chain:** regenerate every dependent record after any upstream file changes. This is more than invoking the final assembler.

Do not run fixed-output scripts in the only preserved copy. Use a second extracted/copy of the whole project tree. Do not bypass hash assertions, alter stored `PASS` flags or relabel a fresh result with an old digest. A hash proves identity, not correctness.

Important traps:

- `verify_fixed_tail.py --assemble-only` **writes** `fixed_tail_certificate.json` with `replayed_components_this_run=false`. The older complete assembler requires this flag to be true. It is not a harmless read-only check of the canonical tail ledger.
- `verify_complete_upper_bound.py --assemble-only` also writes an output unless `--out` points elsewhere.
- `assemble_family.py`, finite-vector producers, finite-family merger and most audit scripts have hard-coded canonical output names.
- Re-running `verify_moment_tables.py` preserves the moment tables but overwrites its verification record, changing its hash because runtime metadata changes.
- Re-running the final assembler is not an independent re-evaluation of its input interval enclosures.

## 3. Dependency graph

The mathematical partition is disjoint in spin: head \(0,2,4,6\); finite block \(8,10,\ldots,40\); tail \(42,44,\ldots\). Every block covers all \(\mu>4\).

```text
fixed150-row witness + exact parameter conventions
  |
  +-- exact Roy source sign and Bernstein source caches
  |     +-- compact V tail: 4<mu<=49, all even ell>=42
  |     +-- joint UV V tail: mu>=49, all even ell>=42
  |     +-- real A40 tail: beta in[0,1], all even ell>40
  |            --> fixed_tail_certificate.json
  |
  +-- certified polynomial source moment tables (N32, N64, UV)
  |       --> all_moment_tables_verified.json
  |       +-- original head0/2/4/6 certificates
  |       +-- corrected finite V kappa14 partition and full replays
  |
  +-- finite real A8..40<.235 partition
         +-- original local finite block and explicit1e-6 replay slack
                --> complete_upper_certificate.json (old2.710793)

complete_upper_certificate.json + original head partitions
  --> head_family.py at256/320bits, all8 trial slopes
  --> refined b=0 heads, subdivision4, at256/320bits
  --> preserved head audit records

corrected original finite V/A partitions + moment/cache dependencies
  --> finite_v_vectors.json at280bits (8804 cells)
  --> finite_real_vectors.json at280bits (19790 cells)
  --> combine_finite_family.py at280bits
  --> finite_vectors_independent_audit.json

fixed_tail_certificate.json + unchanged exact caches
  --> certify_tail_4b.py (actual200-bit envelopes)
  --> review_tail_family.py (actual320-bit positive-case replay)

head families + refined heads + finite family + tail4b + replay/audit records
  --> assemble_family.py
  --> complete_directional_certificate.json (seven final inequalities)
```

The old complete certificate is a provenance and input gate, even though the latest numerical total uses a sharper finite block and refined heads. `head_family.verify_dependencies()` validates every hash in the old complete ledger. Both finite-vector producers validate the old `V8_40_kappa14_local_replay.json` manifest. Thus replacing a base cache or moment table invalidates more than the final tail alone.

Key prerequisites include:

- `output/certification_head/complete_upper_certificate.json`;
- `output/certification_tail/fixed_tail_certificate.json` and its complete manifest;
- `output/certification_head/all_moment_tables_verified.json` and its nine tables;
- original head files `spin0_component_tight.json`, `spin2_tight_component.json`, `spin4_tight_component.json`, `spin6_tight_component.json`;
- `output/certification_head/V8_40_kappa14_exact.json`, its independent and local replay records;
- `output/certification_tail/finite_A8_40_tighter.json`;
- `output/certification_tail/source_cache_v2_depth10.pkl`, `compact_source_cache.pkl` and `expanded_roy_sign.json`.

The latest ledger contains 102 manifest entries in the audited snapshot, corresponding to84 distinct normalized file paths; some appear with both slash conventions. The earlier `CERTIFIED_UPDATE_README.md` describes the older2.710793 result; it is useful for prerequisites but is not the latest numerical claim.

## 4. Infinite-spin source data and sign proof

### 4.1 Exact source positivity

`research/tail_certification/certify_expanded_roy_sign.py` constructs the combined Roy polynomial from the exact multipliers. It proves positivity on

\[
26/5\le s\le119/10,\qquad 0\le z\le1-\frac{6/5}{s-4},
\]

equivalently the relevant source region \(t\le-3/5\). It clears \((s-4)^6\), maps the region to a square, and checks positive rational Bernstein coefficients after alternating dyadic subdivisions. Degree is \((20,6)\); the saved witness passes with2 leaves at depth1. The maximum permitted depth is20. The positive lower bound, before the irrelevant positive folded factor2, is about \(2.4279240847\times10^{-6}\), saved as an exact fraction as well as a display float.

This is one source-sign lemma, not a complete \(V/p\) certificate.

```powershell
& $Python research/tail_certification/certify_expanded_roy_sign.py output/angles_resummed/angle_000p000.npz --output output/certification_tail/expanded_roy_sign.json
```

Pass the witness explicitly: this script's historical default witness and output paths are not the final ones.

### 4.2 Exact Bernstein caches

`certify_uv.py` constructs \(K_{1,\infty}=2(s-c)(c-u)\), \(K_{2,\infty}=2(t-c)^2\), \(c=4/3\), contracts them with the full source polynomial, and computes exact rational signed integrals and Bernstein bounds on absolute/positive parts. The retained cache is **v2**, depth10, with1024 source boxes per polynomial. The positive-source correction integrates both its mass and its \(|t|\) moment, rather than using only a cell's worst \(|t|\) endpoint.

`compact_source_data.py` uses depth10 and40 exact source moments. It also bounds the positive/negative source masses and endpoint-vanishing strip derivative norms. Representative conservative masses are Roy absolute \(<5.618\), combined strip absolute \(<0.2401\). The strip estimate bounds each window separately as well, because each nonnegative contribution is at most their sum.

The v2 cache binds witness and builder hashes. The compact cache binds witness and its builder dependency hash. Do not substitute the historical v1 cache. `--rebuild` in the wrapper reruns the components, but the UV producer is allowed to reuse a matching valid exact cache; it does not necessarily recompute source Bernstein coefficients from zero.

```powershell
& $Python research/tail_certification/certify_uv.py --witness output/angles_resummed/angle_000p000.npz --depth 10 --L 42 --M 49 --dynamic --step .004 --output output/certification_tail/uv_L42_M49_v2.json
& $Python research/tail_certification/compact_source_data.py
```

For a genuinely uncached reconstruction, use a disposable tree lacking the corresponding cache, or call the documented `source_data(witness,10)` and `compute(witness,depth=10,N=40)` functions in a new driver. Do not remove the only preserved cache.

## 5. Infinite-spin absorptive sign: every even \(\ell\ge42\)

The proof is written in `research/tail_certification/PROOF.txt`. It uses \(Q=\mu-c\), \(\delta=\mu-4\), \(Z=Q/\delta\), \(p_\ell=P_\ell(Z)/Q^3>0\), and

\[
\frac{X_4}{p_\ell}=\frac{2Z^2P_\ell''-3ZP_\ell'}{Q^2P_\ell}
\ge\frac{\ell(2\ell-5)}{Q^2}>0.
\]

The fixed total X4 multiplier is negative. The inequality direction would be wrong if this lower bound were multiplied by a positive added coefficient alone; positive-direction tails must use the modified **total** negative coefficient.

### 5.1 Compact energy \(4<\mu\le49\)

`certify_compact.py --L42 --step .02` covers closed intervals of width at most1/50, with an explicit boundary at12. It uses200 cells on[4,8] and2050 on[8,49]. There is no omitted threshold slice. At the initial frozen run, Arb precision is160 bits.

For4<mu<=8, source positivity in the nonphysical angular sector, divided-difference bounds and the X4 term yield a polynomial-in-spin envelope. Both its value atL42 and its spin derivative must be negative; this covers every larger spin. For8<=mu<=49, a40-term expansion of the signed \(K_2\) moment is enclosed with geometric remainder

\[
\mathrm{norm}_2\frac{r^{80}}{1-r^2},\qquad r=\frac{317/60}{Q}<1.
\]

The Roy principal-value residue vanishes exactly at these spins. Strip endpoint tapering permits a uniform PV bound \(2D/(B-A)\). Physical Legendre estimates and monotonicity of \(\ell(\ell+1)/P_\ell(Z)\) justify replacing larger spins by L42 in the positive remainder terms.

```powershell
& $Python research/tail_certification/certify_compact.py --L 42 --step .02 --output output/certification_tail/compact_L42.json
```

The saved upper bounds are approximately−543.68381 on the low interval and−34.600213 on the middle interval. The low-energy spin derivative upper is about−79.808666. These are component diagnostics; the final proof uses the full outward intervals.

### 5.2 Unbounded energy \(\mu\ge49\)

`certify_uv.py` covers the joint variable

\[
\tau=\frac{\ell(\ell+1)}{\mu-4}.
\]

The actual settings are `--depth10 --L42 --M49 --dynamic --step .004`, Arb160. There are2000 tau intervals of width1/250 on[0,8],248 width-one intervals on[8,256], and an analytic closure for everytau>=256:2248 finite cells in total. Each cell[lo,hi] uses the stronger exact lower energy \(M_{cell}=\max(49,4+42\cdot43/hi)\). Thus this covers moving high-energy/high-spin limits, not only fixed-energy spin asymptotics.

The exact finite-energy corrections compare the kernels to their UV limits. In particular, one valid relative \(K_1\) error is

\[
\epsilon_1\le\frac{8/M+2153/M^2+14088/M^3}{(1-48/M)(1-20/(3M))}.
\]

The implementation takes the minimum with another independently valid ratio bound. \(K_2\) uses the exact rational envelope in `errors()`. Harmonic UV limits are combined with their signs before the prefactor error is bounded. Rational cosine terms have a uniform \(O(M^{-2})\) real-part error from the convergent harmonic-series derivative bound, with factor

\[
\frac{r(1+r)}{(1-r)^3}+\frac{4r}{(1-r)^2}+\frac8{1-r}.
\]

At fixedtau, positive Legendre coefficients bound all spins between \(P_{42}(1+2c\tau/[42\cdot43])\) and \(I_0(2\sqrt{c\tau})\). Sign-aware division uses the lower denominator for a positive numerator and upper denominator for a negative numerator. Derivative envelopes bound X4. Forlargetau the proof drops the favorable negative X4 term and uses monotonicity plus the finite numerator maximum.

The retained worst normalized V upper is about−0.622570850675398, at tau[1.512,1.516]; the unbounded-tau upper is below−4.77654. Combining compact and UV regions permits the simple exact global choice

\[
V_\ell/p_\ell\le-31/50\qquad(\ell\ge42\text{ even},\;\mu>4).
\]

## 6. Infinite-spin real norm and complete tail cost

`research/certify_real_tail.py` bounds

\[
A_{40}=\sum_{\ell=42,44,\ldots}\int_4^\infty U_\ell^2/p_\ell\,d\nu_\ell.
\]

Notice the parameter is **L40** for this script because its domain is evenell>L, whereas the sign scripts use **L42**, meaning ell>=L.

For each integer2<=k<=30, exact positive Legendre coefficients and a telescoping even-spin sum give

\[
\sum_{\ell>40,\,even}\frac{2\ell+1}{P_\ell(Z)}\le C_k(\mu-4)^k,
\quad C_k=\frac{(k!)^2(3/4)^k}{2(k-1)\prod_{r=43-k}^{40+k}r}.
\]

This eliminates the infinite spin sum analytically. Any candidate k gives a valid upper bound; choosing the numerically smallest valid one affects sharpness only. The energy variable beta=sqrt(1−4/mu) covers exactly[0,1]. The last cell uses onlyk2 or3 so the infinity endpoint has no divergent intermediate factor.

The forward imaginary combination is evaluated after dividing by its known endpoint zero. Its harmonic polynomial uses a degree-five midpoint Taylor expansion with a sixth-order derivative remainder; rational terms use a positive completed-square denominator. Strip projections are controlled by angular L2/Bernstein bounds. The final route is `--combined-strips`, combining the absolute real fields before squaring on each energy cell. The printed `strip_A_upper_separate_crude` is a large alternative bound and is **not added** to the combined result.

The exact command is:

```powershell
& $Python research/certify_real_tail.py output/angles_resummed/angle_000p000.npz --L 40 --target 2e-5 --max-cells 8000 --precision 160 --combined-strips --out output/certification_tail/real_A40.json
```

It begins with256 dyadic beta intervals; bisection is guided by cell mass, while each accepted mass is a rigorous interval upper density times exact width. Nonfinite enclosures subdivide, with minimum width2^−30 before an unresolved error. The saved run has3617 cells and

\[
A_{40}\le\frac{957792012329886398729379805777930133606695589535}{2^{175}}
\approx1.999963039190242\times10^{-5}<2\times10^{-5}.
\]

No quadrature-convergence estimate is used here. The older A100 certificate is not needed by the final L42 assembly.

For \(V\le-\kappa p\), the disk support satisfies

\[
H(U,V)=\sqrt{U^2+V^2}+V\le U^2/(2\kappa p).
\]

The fixed g0 tail cost is \(A_{40}/(2\cdot31/50)\), approximately0.0000161287341870181, strictly less than1/62000. This includes every evenell>=42 at every energy. It uses no external g2 cap.

### 6.1 Full frozen-tail wrapper

In a disposable working copy, this runs the source sign, UV envelope/cache, compact cache/envelope and real norm, then assembles them:

```powershell
& $Python research/tail_certification/verify_fixed_tail.py --rebuild
```

The default without flags is also a full replay. The wrapper explicitly supplies the final L/depth/M/grid/target settings. Standalone defaults generally do not match the final certificate. It checks the exact witness, domains, source sign, every relevant envelope against−31/50, A<1/50000, and the exact tail quotient. Its final manifest binds all component sources, proofs, caches and records.

## 7. Finite block \(8,10,\ldots,40\): exact source kernels

The mathematical derivation is in `research/tail_certification/FINITE_HEAD_COMPONENTS.txt`; later local and per-spin refinements improve the budget, not the kernel identities.

### 7.1 Source moments and their verification

The finite-V evaluator reads three N32 tables:

| Table | Original precision | Certified Gauss order | Polynomial bounds |
|---|---:|---:|---|
| `roy_moments.json` |220|48|jmax40, kmax32, skmax40; degrees94 energy,77 angle|
| `strip_12_24_moments.json` |240|40|jmax40,kmax32; degrees38 energy,73 angle|
| `strip_24_48_moments.json` |240|40|same|

Gauss nodes and weights are interval-certified; the integrands are polynomials of degree at most2n−1, so polynomial exactness removes the quadrature truncation error. Interval errors in roots, weights and arithmetic remain enclosed. The separate source geometric expansions below are still approximations and have explicit remainders.

Roy moments are \(\int W t^j(u+4)^k\), \(\int Wt^j(s-4)^a\), \(\int Wt^j\), \(\int W(t+2)^k\). Roy W includes the folded angular factor2. Strip transfer center/radius are both2/3; u centers/radii are40/3,20/3 for[12,24], and94/3,38/3 for[24,48]. Energy centers/radii are18,6 and36,12.

`verify_moment_tables.py` rebuilds all nine N32/N64/UV tables at340bits and requires each fresh ball to be wholly contained in the corresponding saved ball. Merely overlapping intervals do not pass. The saved verification covers24344 entries. This validates the actual stored tables rather than replacing them with new ones.

```powershell
# Working copy only: this overwrites its verification record, not the tables.
& $Python research/verify_moment_tables.py
```

For generating N32 tables in a new production tree, the retained interfaces are:

```powershell
& $Python research/certified_roy_moments.py output/angles_resummed/angle_000p000.npz --jmax 40 --kmax 32 --skmax 40 --out output/certification_head/roy_moments.json
& $Python research/certified_strip_moments.py output/angles_resummed/angle_000p000.npz --jmax 40 --kmax 32 --outdir output/certification_head
```

N64 and UV tables additionally serve the low-spin heads; their reconstruction belongs with the head workflow. Replacing any table requires regenerating its verifiers and downstream hash-bound records.

### 7.2 Finite V evaluation and remainder bounds

`certified_finite_v.py` normalizes by \(R_\ell=\delta^\ell P_\ell(1+2t/\delta)\), whose coefficients are exact integers. This avoids threshold overflow. Forell>=8 the Roy source's degreeJmax6 makes the PV residue vanish exactly, leaving a polynomial contraction with the saved moments. Removable crossing denominators are divided out algebraically, grouping equal moment coefficients before ball evaluation.

N32 geometric expansions handle \(1/(\mu-u)\) and \(1/(\mu-t)\). Their generic remainder is

\[
\frac{q^{32}}{D(1-q)},\qquad q=R/D<1,
\]

multiplied by the appropriate positive source norm and normalized Legendre envelope. Roy uses u center/radius4/4 and t center/radius2/2. For Roymu<159/20=7.95 the nonphysical angular ratio is explicitly enclosed; otherwise physical Legendre magnitude<=1 suffices. Strips always use physical angles.

For strips belowmu101, the energy PV is evaluated with exact Legendre-Q/polynomial formulas. Tapering removes endpoint divergences; intervals crossing an endpoint use the continuous xlog|x| form, with an enforced local |x|<1/4 condition. For the UV branch the strip energy denominator is also expanded toN32 with its explicit remainder.

After common-denominator assembly, the top two powers cancel analytically because the full kernel is O(mu^−3). Code checks the corresponding coefficient balls contain zero, then removes those coefficients **using the proved identity**, not by deciding a small float is zero. The reciprocal variable \(w=1/(\mu-4)\) evaluates the UV polynomials atw0, including infinity exactly.

### 7.3 Corrected original sign partition

The reusable seed is `V8_40_kappa14_exact.json`:8798 accepted leaves, all17 spins finite and below−14 on every leaf. The saved worst upper is approximately−14.000004870715432. The domains are exact closed covers ofmu[4,101] andw[0,1/97], the latter representingmu[101,infinity]. Initial mu cuts are4,5,8,119/10,12,24,48,49,100,101.

```powershell
# Fresh sign partition, in a disposable working tree.
& $Python research/tail_certification/certify_finite_v_domain.py --kappa 14 --maxcells 100000 --output output/certification_head/V8_40_kappa14_exact.json

# Or start from preserved cells and adaptively subdivide any failures.
& $Python research/tail_certification/certify_finite_v_domain.py --kappa 14 --maxcells 100000 --replay output/certification_head/V8_40_kappa14_exact.json --output output/certification_head/V8_40_kappa14_recomputed.json
```

Precision is hard-coded260 bits; maximum subdivision depth35; the maxcells flag bounds trials, not just accepted leaves. Exact acceptance is **every spin finite and its upper<-kappa**. The old2744-cell `V8_40.json` and the pre-correction kappa14 file are rejected historical results. They must not be restored as certificates. A float maximum once hid a nonfinite spin entry; the corrected evaluator tests all entries explicitly.

The old complete ledger additionally requires a full280-bit replay and a full local-upper replay with explicit1e−6 allowance:

```powershell
& $Python research/replay_finite_v_for_assembly.py output/certification_head/V8_40_kappa14_exact.json --out output/certification_head/V8_40_kappa14_independent_replay.json
& $Python research/replay_local_finite_v.py output/certification_head/V8_40_kappa14_exact.json --slack 0.000001 --out output/certification_head/V8_40_kappa14_local_replay.json
```

Each covers149566 spin-cell inequalities. The local replay checks freshupper<=storedupper+1e−6. This allowance is applied once in the older aggregate local-block assembly. It is not a universal fudge factor added to the latest freshly computed vector data.

### 7.4 Original real-norm partition

`certify_finite_real.py` computes \(\sum_{\ell=8}^{40}\int U_\ell^2/p_\ell\) onbeta[0,1] by180-bit upper Riemann sums. Initial256 dyadic intervals are refined by largest cell mass. Its exact physical strip projection is a terminating polynomial \(H_{\ell j}(c/(\mu-4))\); interval polynomial evaluation is shifted around the midpoint. Source-window boundary intervals include the union with zero. \(1/P_\ell(Z)\) is bounded at the beta cell's right endpoint by monotonicity.

```powershell
& $Python research/tail_certification/certify_finite_real.py --first 8 --last 40 --target .235 --max-cells 30000 --output output/certification_tail/finite_A8_40_tighter.json
```

The saved target is47/200=.235,19790 cells, bound about.2349979786611204. The archived record does not retain the original command's max-cells argument;30000 is an explicit sufficient budget for reproducing the target/partition, not a claim about that omitted metadata. Defaults `.25` and12000 would produce a different, weaker seed.

## 8. Latest per-spin finite block

The new producer takes the original partitions as seeds but recomputes the interval quantities per spin at280bits:

```powershell
& $Python research/directional_certification/certify_finite_vectors.py real
& $Python research/directional_certification/certify_finite_vectors.py v
& $Python research/directional_certification/combine_finite_family.py
```

These have no output or precision flags. Use a disposable copy. They read the old local replay manifest and reject stale prerequisites.

`real` produces19790 beta cells, with a directed real mass for each of17 spins. `v` starts from8798 V cells and subdivides six additional times to8804 leaves. It stores separately an upper bound forVold/p andX4/p for every spin. Maximum extra subdivision depth is25. Threshold and UV ratios are computed as exact rational functions ofdelta andw, respectively.

The `v` producer's trial metadata usesd=32b/9 forpositiveb because it imports the older base case list. **This does not define the final certificate.** The final merger owns an explicitd4b list and recomputes every affine sign

\[
k_{\ell,I}=-(V_{\ell,I}^{upper}+b+dX_{\ell,I}^{upper})>0.
\]

It checks all17 spins and every spectral cell. The conditiond>=0 is essential when usingXupper. The b1/2 case fails and is explicitly excluded.

For each beta cellJ, transform its endpoints exactly bymu=4/(1−beta²). Take the minimumkappa over **all** intersecting spectral cells separately for each spin, including the unbounded last interval. The block cost is the directed sum ofmassell,J/(2minI kappaell,I). Two closed intervals meeting at one endpoint allow the redundant neighbor to be omitted, because the retained interval covers that endpoint. Infinity is represented by an unbounded final interval, not a finite cutoff.

The retained g0 contribution is approximately.002637191340269, compared with the prior aggregate local budget.0059764297. Bounds for all seven cases are in `finite_block_family.json`. The file has overallstatusPASS even though its per-case list includesoneFAILED_COERCIVITY; inspect individual results.

`review_finite_vectors.py` checks all exact partitions, every saved affine sign, an independently coded interval-overlap merge and directed340-bit summation. It also freshly evaluates130 spectral cells and97 real cells at340bits. This last check is **sampled**, not a second full evaluation of the8804-by17 vector enclosures. The rigorous all-cell enclosures come from the reviewed original280-bit producer. Different valid interval calculations need not yield nested balls.

```powershell
& $Python research/directional_certification/review_finite_vectors.py
```

It writes `finite_vectors_independent_audit.json` and records that limited scope explicitly. It checks sampled excess tolerances only as diagnostics, not as the mathematical certificate's error budget. Its exact merge agrees with the producer to far better than1e−12.

## 9. Directional infinite tail

Use the4b wrapper, not the unwrapped base script:

```powershell
& $Python research/directional_certification/certify_tail_4b.py
& $Python research/directional_certification/review_tail_family.py
```

The wrapper substitutes the same exactb,d pairs as the final head/finite components. It runs the joint UV envelopes withL42,M49,step.004 and the compact envelopes withstep.02. The unchanged source data, moments and source positivity remain valid because onlyb and the X4 multiplier change.

Forb<=0,d0, the original exact margin shifts tokappa=31/50−b. Forpositiveb the total X4 coefficient remains negative; both compact and UV envelopes are reevaluated, thenb is added. Results:

| \(b\) | Producer global kappa, approximately | Assembled policy |
|---:|---:|---|
|−1|1.62|exact31/50−b|
|−1/2|1.12|exact31/50−b|
|−1/4|.87|exact31/50−b|
|−1/10|.72|exact31/50−b|
|0|.62|exact31/50|
|1/10|.4238827564153527|reduce by1e−10 after replay|
|1/4|.1257852055240021|reduce by1e−10 after replay|
|1/2|negative|FAILED_ENVELOPE; not certified|

The shared \(A_{40}\) norm is unchanged. Divide it by2kappa for each direction. The positive-case reductions are exact rational1/10^10, not a displayed-decimal subtraction.

**Actual precision versus stale labels:** `certify_directional_tail.main()` sets the shared Flintctx to200bits. Both its UV and compact calls therefore use200bits, even though the compact function's returned arithmetic string saysArb160. `review_tail_family.py` sets the same sharedctx to320bits; both calls use320bits. This was checked by object identity and readingctx.prec. The replay's wording “compact at frozen kernel precision” is stale. The formulas are the frozen formulas; the actual arithmetic precision is320.

The replay only reruns the two passing positive cases. It binds the input tail4b file and its source hash, proves its fresh upper<=old upper+1e−10<0, and records the weakened kappa. Negative cases need no fresh analytic envelope because the exact shift is monotone.

## 10. Assembly, exact arithmetic and current outputs

With valid preserved prerequisites, the following produce the final directional stages in a disposable tree:

```powershell
# Detailed head-method validation belongs to the head guide.
foreach ($Spin in 0,2,4,6) {
  & $Python research/directional_certification/head_family.py --spin $Spin --precision 256 --out "output/directional_certification/head_spin$Spin.json"
  & $Python research/directional_certification/head_family.py --spin $Spin --precision 320 --out "output/directional_certification/replay_head_spin$Spin.json"
  & $Python research/directional_certification/head_family.py --spin $Spin --slopes=0 --subdivide 4 --precision 256 --out "output/directional_certification/refined_spin$Spin.json"
  & $Python research/directional_certification/head_family.py --spin $Spin --slopes=0 --subdivide 4 --precision 320 --out "output/directional_certification/refined_replay_spin$Spin.json"
}
# New record-audit helper: all32 directional pairs and4 refined pairs.
# These are new records for the actual newly generated inputs, not copied PASS flags.
& $Python research/review_head_replay_records.py --out output/directional_certification/head_family_independent_audit.json
& $Python research/review_head_replay_records.py --out output/directional_certification/refined_heads_independent_audit.json
& $Python research/directional_certification/certify_finite_vectors.py real
& $Python research/directional_certification/certify_finite_vectors.py v
& $Python research/directional_certification/combine_finite_family.py
& $Python research/directional_certification/review_finite_vectors.py
& $Python research/directional_certification/certify_tail_4b.py
& $Python research/directional_certification/review_tail_family.py
# Only after every required audit record is valid for the regenerated inputs:
& $Python research/directional_certification/assemble_family.py
```

The final assembler:

1. Requires exact witness and source hashes and all nested prerequisite hashes.
2. Verifies each head's exact rational[0,1] partition, nonnegative masses, exact total and negative threshold cap.
3. Joins records by exact rational(b,d), not angle labels; accepts only the seven passing finite/tail cases.
4. Takes the larger head total from256/320bit runs. Forb<=0 it may also use the refined g0 cap by monotonicity ofH inV, taking the smaller of two separately valid bounds.
5. Adds the exact stored finite-block upper and exactA/(2kappa) tail using Python Fractions.
6. Searches positive combinations of pairs of certified slope inequalities. This lowers some rows; it is exact convex algebra, not a new quadrature computation. `direct_witness_upper_exact` preserves the pre-interpolation value. For example the finalb=.1 constant is about2.713167121015, while its direct certificate is about2.713644276966.
7. Rounds the final total **upward** at six decimals. Unit-angle supports and angle labels are floating displays only.

The independently computed common decimal caps in audit records are conservative comparison conveniences. The final assembler uses the exact larger verified head totals and exact finite block; it does not add the audits' display cap margins a second time. Likewise the old finite-local1e−6 slack belongs to the old aggregate certificate, not the new vector block. The positive-tail1e−10 kappa allowance is used by the latest assembly.

## 11. Preserved audit history and safe replay procedure

The preserved source is enough to reconstruct and review the mathematics. The new helpers below close a record-generation gap, but this is still **not a turnkey byte-identical full-from-scratch build** of every historical artifact.

- The original frozen `head_family_independent_audit.json` was generated by an inline audit. Its original standalone producer was not retained. The final assembler requires this record and validates its recorded input hashes. Re-running any head file changes the corresponding hash, so the old record cannot simply be reused.
- The old `review_refined_heads.py` explicitly covers spins2,4,6. The latest ledger also includes refinedspin0 and its fresh replay, which the assembler checks directly; do not describe that older audit record as covering refinedspin0.
- **New solution:** `research/review_head_replay_records.py --out PATH` audits all36 pairs (72 records):32 directional spin/slope pairs and4 refinedg0 pairs, including refinedspin0. It checks both precision labels256/320, exact partitions/sums/threshold signs, parameters and provenance, and compares totals within1e−12 as a diagnostic. It does not rerun kernels. The commands in Section10 create honest replacement audit records bound to the newly generated files; the existing assembler accepts their `PASS` plus complete manifest format. The original frozen records need not be modified.
- Several producers write only to fixed filenames. Their timing metadata changes on a replay even when exact numerical outputs agree.
- A higher-precision run is not guaranteed to reproduce narrower stored balls. The declared theorem must use explicit verified enclosures or documented slack; midpoint agreement is insufficient.

**Safe frozen verification:** the new standard-library-only verifier is:

```powershell
& $Python research/verify_replication_records.py
```

It takes no arguments and writes nothing. It checks the expected frozen ledger hash, all102 manifest entries (84 distinct normalized file paths), every saved directional/refined head partition, rational sums and threshold signs, then reconstructs all7 exact totals and positive interpolation. It does **not** reevaluate interval kernels, finite vectors or tail envelopes. A changed ledger intentionally fails the frozen-digest check; this is not a generic verifier for newly generated ledgers. Do not call a writing `--assemble-only` command over the canonical tail record. The final assembler may be run in a separate copied tree using unchanged preserved inputs; its output can be compared mathematically, not necessarily bytewise.

**Safe mathematical replay:** create a disposable copy with all frozen prerequisites. Recompute one branch there. Compare its exact domain cover, successful inequalities, source conventions, remainders and resulting outward bound with the frozen result. Retain the original branch and parent ledgers in the archive. If rebuilding the whole graph, regenerate parent manifests only after fresh child checks pass; do not edit their hashes manually to silence failures.

**Fresh full production:** proceed topologically: exact source tables/caches; their validations; original all-spin tail; original finite/head components and replays; original complete ledger; then directional heads/vectors/tails; new audit records; final directional assembly. The new head-audit helper can now produce the required replacement records. Broader rebuilds still require re-generating every upstream dependent record in the proper order; the Section10 commands deliberately assume the preserved valid base prerequisites. This audit does not claim that all base head components can be regenerated from an empty output directory by that short command block.

The exact head record-audit algorithm is: for all4spins×8slopes plus4refined pairs, verify source/witness/dependency hashes and exact(b,d); verify each cells+gaps partition covers[0,1] without gaps or overlaps; verify nonnegative rational cell bounds and their exact sum; verify the stored negative threshold cap where applicable; compare original/replay totals; save new record hashes. This audit records completed interval runs—it must not manufacture them from copied totals.

## 12. Runtime expectations and scope of review

Observed original wall times, not guarantees:

| Stage | Approximate retained runtime |
|---|---:|
| realA40 combined tail,3617cells |2.7s|
| original finiteA<.235,19790cells |43s|
| corrected kappa14 finiteV,8798leaves |147s|
| new finite real vectors |38s|
| new finite V/X4 vectors |230s|
| new seven-case finite merge |23s|
| directional tail4b |24s|
| independent finite merge plus sampled kernels |32s|
| refined head320bit runs, one spin |roughly200–250s in the parallel audit|

Fresh exact-cache construction and all head runs add time; do not promise these timings across machines. Parallelize independent processes only when they write distinct files. Do not let a component overwrite an input while another process is hashing or reading it. Audit scripts sometimes cache imported kernel data within a process; start a fresh process after changing any input.

The latest producer rigorously encloses all finite vector cells. Independent review included all exact partitions/merges and sampled fresh vector evaluations. Earlier aggregate kappa14 V had a full all-leaf replay. Head and positive-tail replays largely share the same kernel implementation. There is no claim that every formula has two independent implementations or formal verification. Source formulas and analytic inequalities were reviewed separately, and remain part of the trusted mathematical/code base.

Historical numerical curves, the enriched266-row candidate, old fractional-kernel trials, old rejected V partitions, and rounded published graphical outlines are not inputs to the seven certified inequalities. Reproducing the plot does not reproduce the proof. Reproducing the proof does not show the witness is optimal or that an amplitude saturates it.

For further mathematical details use `research/tail_certification/PROOF.txt`, `research/real_tail_certificate_derivation.md`, `research/tail_certification/FINITE_HEAD_COMPONENTS.txt`, `research/head_certification_spin0_PROOF.txt`, the frozen kernel sources, and `FURTHER_DIRECTIONAL_CERTIFICATION_2026-10-06.md`. Their old component-only warnings remain useful scope statements; the latest complete ledger is the final assembled claim.
