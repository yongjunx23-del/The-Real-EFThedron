"""Build point-scope mathematical appendix and immutable search provenance."""
from pathlib import Path
import re, json, hashlib

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p): return str(p.relative_to(ROOT)).replace('\\','/')
def path(s):
    p=Path(str(s).replace('\\','/'))
    return p if p.is_absolute() else ROOT/p

catalogue=(ROOT/'ALL_CONSTRAINTS_HUMAN_READABLE_2026-10-07.md').read_text(encoding='utf-8')
sections={int(m.group(1)):m.group(2) for m in re.finditer(r'^## (\d+)\. (.*?)(?=^## \d+\.|\Z)',catalogue,re.M|re.S)}
def section(n):
    s=sections[n].strip()
    return '## '+s

intro=r'''# Constraints and construction of the dual witnesses

This appendix defines the constraints used by the certified boundary in the final two-panel figure. The 39 active half-planes of the upper panel include the 22 distinct inequalities represented in the lower panel (the minimum-\g_0 inequality is displayed at both ends of the perimeter). Each half-plane is an exclusion inequality, not a constructed scattering amplitude. The accompanying registry ties every inequality and plotted position to its exact normal, multiplier vector, and selected proof packet. The smooth connecting segments in the lower panel are display guides and introduce no new certified inequalities.

There are two physical row banks, of sizes 150 and 220. The analytic inequality \g_2\ge0 is also included. The counts specify retained identities, not linear independence or a complete characterization of physical amplitudes. No independent Gram or SOS coupling constraint is used. The rational source-polynomial inequalities below constrain a dual witness during its search; they do not impose new assumptions on amplitudes.

The certificates use the frozen coefficients exactly. A binary64 coefficient denotes its exact dyadic rational, while a rational combination of witnesses retains its exact rational weights. Nominal angles identify records only: the rational normal supplied with the inequality is authoritative. Multipliers in the data file are expressed in the physical rows defined below, after undoing every numerical conditioning transformation.
'''.replace('minimum-\\g_0','minimum-\\(g_0\\)').replace('inequality \\g_2\\ge0','inequality \\(g_2\\ge0\\)')

bank=r'''## The two retained physical row banks

Every range includes all Cartesian-product combinations. The ordering is the saved row-name ordering, retained in the exact multiplier export.

| Family | Bank A: 150 rows | Bank B: 220 rows |
|---|---|---|
| Forward harmonics | \(H_j,\ j=2,\ldots,40\):39 | Same:39 |
| Rational cosine sums | \(r=4/5,19/20,\ \zeta\in C_7\):14 | \(r=4/5,\ \zeta\in C'_7\):7 |
| Rational sine sums | None | \(r=4/5,\ \zeta\in C'_7\):7 |
| Crossing | \(X_4\):1 | \(X_4,X_5,X_6,X_{7,1},X_{7,2},X_{8,1},X_{8,2}\):7 |
| Folded Roy | \(B=119/10,\ J=0,2,4,6,\ n=0,\ldots,14\):60 | Same:60 |
| Tapered strips | Two windows, \(j=0,1,2,\ n=0,\ldots,5\):36 | Two windows, \(j=0,\ldots,3,\ n=0,\ldots,7\):64 |
| Regular transfer derivatives | None | \(p=1,2,3\), each \(H_j,\ j=2,\ldots,12\), and one \(G_{\rm regular}\):36 |

The strip windows are \([12,24]\) and \([24,48]\), with transfer width \(T=4/3\) and taper power one. The parameter sets are exactly
\[
C_7=\{-1,-3/4,-1/2,0,1/2,3/4,1\},\qquad
C'_7=\{-9/10,-3/4,-1/2,0,1/2,3/4,9/10\}.
\]
Bank A supplies the right-facing family obtained from the original \((1,0)\) witness by exact objective shifts and, for one selected inequality, a positive combination of two certified costs. Every selected Bank A inequality retains the same original multipliers; none of the older positive-slope crossing repairs is used in this figure. Bank B supplies the independently searched left, upper, and lower witnesses and their rational mixtures. A mixture combines matching row definitions with exact coefficients; it does not introduce another constraint family.
'''

# Formula text is kept verbatim, except scope and cross-reference wording.
formula_parts=[section(2),section(3),bank]+[section(i) for i in range(5,12)]
formulas='\n\n'.join(formula_parts)
formulas=formulas.replace('Catalogue B','Bank B').replace('§7','the crossing subsection').replace('§11','the transfer-derivative subsection')
formulas=formulas.replace('These 36 rows occur only in Bank B.','These 36 rows occur only in Bank B; the derivative-order symbol \\(p\\) in this subsection is distinct from the positive kernel \\(p_\\ell\\).')

search=r'''## Numerical optimization and the complete search restrictions

The optimization minimizes a finite quadrature approximation to \(\langle H(U,V)\rangle\) over the free real multipliers. The full physical disk is eliminated by its exact support function rather than sampled in phase. The finite objective omits positive contributions from spins above its training range, so its minimum is not an exclusion bound. Certification subsequently encloses the complete fixed-witness integral. Neither a finite iteration budget nor a solver's successful quadratic subproblem proves global optimality of the reported witness.

For columns \(k\) with nonnegative quadrature weights \(W_k\), the numerical convex objective and its derivatives, away from \(R_k=\sqrt{U_k^2+V_k^2}=0\), are
\[
J(\lambda)=\sum_k W_k H(U_k,V_k),\qquad
\partial_iJ=\sum_kW_k\left[\frac{U_ku_{ik}+V_kv_{ik}}{R_k}+v_{ik}\right],
\]
\[
\partial_i\partial_jJ=\sum_k\frac{W_k}{R_k^3}
(V_ku_{ik}-U_kv_{ik})(V_ku_{jk}-U_kv_{jk}).
\]
The floating evaluator computes \(V/p\) before multiplying by \(p\), uses the cancellation-free support expression for negative \(V\), and floors its numerical radius at \(10^{-300}\). These devices stabilize discovery; the certificate encloses the unfloored exact support.

### Training quadrature

The original Bank A search includes even spins \(0,2,\ldots,24\), with 64 Gauss--Legendre nodes on each of nine velocity panels. The velocity is \(\beta=\sqrt{1-4/\mu}\), and the sorted endpoint set is
\[
\{0,.12,.35,\sqrt{1-4/(119/10)},\sqrt{1-4/12},
\sqrt{1-4/24},\sqrt{1-4/48},.95,.995,1\}.
\]
The exact measure in this chart is \(4(2\ell+1)d\beta/[\pi(1-\beta^2)^2]\). A zero-weight point \(\mu=4+10^{-10}\) is appended at each training spin. There are 7,501 columns, of which 1,891 satisfy the high-\(p\) guard mask. The source quadratures use Roy orders 48 in energy and 40 in angle and strip orders 40 in energy and 32 in transfer; the taper is assembled by its finite Legendre expansion.

The ordinary Bank B cache uses every even spin through 40, 48 spectral nodes per panel, and source/angular orders 64. At spins 26--40 the first threshold panel uses 20 nodes. The denser upper-sector cache uses 192 spectral nodes and source/angular orders 96, with the same 20-node exception in the first panel at spins 26--40. Its saved metadata records 1,729 energy columns for each spin 0--24 and 1,556 for each spin 26--40, including the zero-weight probe. Principal values and removable poles are treated by the formulas above in both searches. These quadratures are numerical approximations used for discovery; the proof does not assume their error vanishes.

### Common sign guard and its normalization

On each training or zero-weight probe column with \(p_\ell(\mu)>10^4\), impose \(V/p\le-\varepsilon\), with \(\varepsilon=2\times10^{-6}\) for Bank A and \(10^{-8}\) for Bank B. A spin-zero probe need not satisfy the mask or have negative \(V\). The original Bank A search has no additional Roy-source samples, independent crossing-sign column, or added-spin grid.

Additional Bank B multiplier guards are zero-weight artificial columns, not additional amplitude constraints. For a raw affine guard \(e+d\cdot\lambda\le-m\), put \(s_d=\max(1,\max_i|d_i|)\) and encode
\[
v_{\rm art}=\frac{e+m+d\cdot\lambda}{s_d},\qquad
U_{\rm art}=0,\quad p_{\rm art}=10^5,\quad W_{\rm art}=0.
\]
The common sign guard therefore enforces \(e+d\cdot\lambda\le-m-\varepsilon s_d\). The highest-crossing and source guard definitions below specify their own affine normalization; the crossing-series constant term includes an explicit \(-\varepsilon\) offset to cancel this cushion. The exact numerical matrices are preserved with the run metadata rather than inferred from rounded parameter labels.

### Sampled Bank B source and spectral guards

All selected ordinary Bank B searches retain the following:

1. A Roy-source screen with 41 equally spaced \(s\in[20/3,119/10]\) and 21 equally spaced points \(z\in[0,1-(8/3)/(s-4)]\) at each \(s\), including coincident endpoint samples. It requires a positive requested source margin for
\[
W_R(s,z)=\sum_{J,n}\lambda_{BJn}P_n\!\left(2(s-4)/(B-4)-1\right)P_J(z).
\]
2. A highest-order crossing coefficient guard \(\sum_i\lambda_i\beta_{i,6}\le0\), with the common strict cushion. This is not a constraint on the sign of \(\lambda_{X_4}\) alone.
3. The total-field guards in the table below, using the actual \(a v_0/p+b+\sum_i\lambda_i v_i/p\). The first group has requested margin zero; the remaining groups use the run's tail margin.

| Even spins | Spectral energies |
|---|---|
| 26--40 | \(4+10^{-6},4+10^{-3},4+10^{-2},4.1,4.5,5,6,8,12.1,20,30,50,75,100,150,250,500,1000,10^4,10^6\) |
| 42--72 | \(13,15,\ldots,103,150,250,500,1000,10^4,10^6\) |
| 80,100,140,200 | \(8,10,12.1,16,20,30,50,75,100,150,250,500,1000,10^4,10^6\) |

4. Joint spin/energy guards \(\mu=4+\ell(\ell+1)/\tau\), for even \(\ell=42,\ldots,100\) and \(\ell\in\{120,160,200,300,500,1000\}\). The \(\tau\) set is the union of 12 geometrically spaced values from \(10^{-3}\) to 1, the arithmetic lists \(1,1+1/8,\ldots,7\) and \(7,7+1/4,\ldots,40\), and 14 geometrically spaced values from 40 to 150; duplicate binary64 values are removed. The tail margin applies to every such guard.
5. All 9,408 Bernstein source inequalities and 42 crossing-series inequalities described next.

A dense sample is not a continuum sign proof. Source and spectral signs used by certification are recomputed for the frozen multiplier vector with exact or directed interval arithmetic.
'''

# Reuse the fully specified polynomial guard mathematics, not unrelated searches.
guard_source=sections[12]
bern=re.search(r'^### 12\.4 (.*?)(?=^### 12\.5 )',guard_source,re.M|re.S).group(1).strip()
cross=re.search(r'^### 12\.5 (.*?)(?=^### 12\.6 )',guard_source,re.M|re.S).group(1).strip()
search+='\n\n### '+bern+'\n\n### '+cross+'\n'
search=search.replace('§7','the crossing subsection').replace('Catalogue B','Bank B')
search=search.replace('21×7×16×4 = 9,408',r'\(21\times7\times16\times4=9{,}408\)')

extra=r'''### Upper-sector additional guards

The selected 85-degree search adds 281 spin-14 guards at \(\mu_k=(1201+10k)/200\), \(k=0,\ldots,280\), with margin \(10^{-4}\). It does not use all spins of the larger cached bank as optimization guards. The stored JSON contains the actual binary64 grid values.

The 83.1759644845-degree and 86.4927421125-degree searches, and the direct upper-95-percent search, add the full 4,496-point bank: all even spins 10--40 at the same 281 intended energies, evaluated with source/angular orders 96, with requested finite margin .002. They also add spins 2 and 4 at \(\mu=4+10^{-8},4+10^{-6},4+10^{-4}\), with margin .001. Their tail margin is .01 and source margin .0001. The upper-83 and upper-86 refinements use Hessian eigenvalue floor \(10^{-9}\); the selected upper-95 search uses \(10^{-7}\).

The selected upper-\(g_2\) vector is an intermediate captured coefficient snapshot of a robust search, not that search's later final vector. Its parent search has threshold margin .001, tail margin .02, finite-bank margin .02, and source margin .0001, using the dense 192/96 cache. The capture record identifies the exact checkpoint and hash. The full 4,496-point bank and the six threshold guards are imposed by this robust driver.

Each row in the saved optimization metadata includes the actual run arguments, known warm start, cache and guard paths, recorded status, and hash. Empty fields mean unrecorded information, not an inferred default. The ordinary left/minimum-\(g_0\) search uses tail margin .002; later 259.294118509649-, 261.3814-, 265-, 268-, 263.2747-, and 269.5-degree selected runs use their recorded .005 tail margin. The tabulated machine-readable run arguments control if a later continuation differs.

### Solver choices and limitations of historical reconstruction

The Bank A constrained Newton method uses a diagonally equilibrated Hessian, eigenvalue floor \(\max(10^{-8}\lambda_{\max},10^{-12})\), trust radius 10, and CLARABEL quadratic subproblems. The recorded stopping budget is 60 accepted iterations; its status is an iteration limit. The retained original tolerances are function tolerance \(3\times10^{-11}\), gradient tolerance \(10^{-6}\), QP gap tolerances \(10^{-10}\), QP feasibility tolerance \(10^{-11}\), and QP iteration limit 250. The line search uses a .995 feasible step fraction, Armijo factor \(10^{-4}\), up to 70 halvings, and residual acceptance \(2\times10^{-8}\). These are algorithm settings and do not enter the proof as amplitude assumptions. The initial warm-start filename for this historical Bank A search was not retained in its compact run record; it is explicitly unavailable. The final exact witness is fully retained.

The Bank B searches use diagonal equilibration and whitening, the run-specific eigenvalue floor, a trust-region Newton quadratic model solved by CLARABEL, and a scaled \(L^1\) feasibility projection solved by HiGHS. Constraint generation starts from nearby rows, adds violations found against the complete bank, and checks the full bank before accepting the step. The inward QP cushion is \(10^{-9}\) in the normalized constraints. A persistent active set changes computational efficiency, not the final constraint set. The saved scripts and argument records specify iteration limits and warm starts. Multiple refinements and the initial cache construction must not be conflated with one timing number.

In the supplied Bank B driver the default trust radius is 10 in whitened step coordinates, not a bound on the multipliers. It increases by a factor 1.4 up to 100 after a step fraction exceeding .9, and decreases by a factor .5 down to .05 when the fraction is below .1. The QP has iteration cap 250, gap tolerances \(10^{-10}\), and feasibility tolerance \(10^{-11}\); the supplied implementation requests two solver threads. The .995 feasible fraction, Armijo factor \(10^{-4}\), and at most 70 halvings are retained. For \(\varepsilon=10^{-8}\), the full-grid line-search residual tolerance is \(\max(.01\varepsilon,10^{-13})=10^{-10}\). Default function and descent tolerances are \(10^{-11}\) and \(10^{-7}\), respectively. The equilibrated Hessian eigenvalues are floored at the greater of the run-specific relative floor times the largest eigenvalue and \(10^{-12}\). The source and kernel caches, complete guard matrices, and original argument blocks specify the finite problem, while coefficient hashes specify the independently checkable result.

Reoptimizing may produce a different coefficient vector because the problem is ill-conditioned and the driver evolved. The archive supports reproduction of the mathematical finite search and exact verification of the supplied fixed witness; it does not assert a bit-identical solver trajectory or optimal multipliers. The proof does not depend on recreating that trajectory.

### Exact derived witnesses and mixtures

Several inequalities are constructed from searched witnesses without a new optimization. For nonnegative exact \(t_i\), define
\[
n=\sum_i t_i n_i,\qquad \lambda=\sum_i t_i\lambda_i,
\quad U=\sum_i t_iU_i,\quad V=\sum_i t_iV_i.
\]
Convex positive homogeneity of disk support implies
\[
H\!\left(\sum_i t_iU_i,\sum_i t_iV_i\right)
\le\sum_i t_iH(U_i,V_i).
\]
The exact parent weights need not sum to one. A sum of parent certified costs is a valid initial bound. Reintegrating selected spins of the combined functional can reduce this cost; only those parent components are replaced. This nonlinear support improvement uses the same physical identities. The paired parent records and exact rational weights are retained, so rounding an exported mixture back to binary64 is not permitted.

Every selected Bank A inequality keeps the original multiplier vector and changes the target to \((1,b)\), with \(b\le0\). Thus \(V_b=V_0+bp\) and \(U\) is unchanged. In the saved indexing, cut C06 at slope \(-1/4\) uses one half of the certified slope-\(-1/2\) inequality plus one half of the slope-zero inequality; the effective multipliers are again the same original vector. All other selected Bank A records use their direct objective-shift witness. No selected Bank A inequality changes the \(X_4\) coefficient. The shifted 269.04041748833-degree Bank B record likewise keeps the 268-degree coefficients and the same \(g_0\) component, replacing only the \(g_2\) component by a more negative exact dyadic value. Monotonicity may transfer valid pieces of an existing certificate, but the complete record must still account for the new objective and every support contribution.

The upper-83 witness uses binary64 coefficients obtained from .999 times its refined vector plus .001 times its earlier vector; the frozen repaired vector, rather than a decimal recomputation of that mixture, defines the certificate. Its exact crossing sign was checked after repair. The direct upper-95 witness is byte-identical to its selected search candidate; trial repair files are not used. The upper-55 and upper-85.8424675 records and the lower-15, lower-45, and lower-60 records are exact rational combinations of their indicated parents, with selected spin contributions reintegrated. The registry distinguishes these constructions from independent searches.

### Exhaustiveness and inputs excluded from the optimization

All higher EFT coefficients remain free. The selected witness banks contain no fractional-power rows, odd-spin partial waves, odd-\(J\) Roy tests, additional generic nonforward harmonic tests, crossing rows beyond degree eight, independent sign inequalities on amplitude derivatives, prescribed threshold scattering lengths, elastic equality, hard spectral cutoffs, externally imposed \(g_2\) caps, or truncated primal amplitude ansatz. The observed primal and dual reference curves are used only to compare or place plotted points; they are not optimization constraints or proof inputs. None of the approximate graphical continuations is part of a certificate.
'''

text='\n\n'.join([intro.strip(),formulas,search,extra])+'\n'
text=text.replace('The user PDF is','The source manuscript is')

# Preserve exact run records while keeping the appendix readable.
run_files=[
'output/angles_resummed/angle_000p000.json',
'output/full_direction_20261007/left/angle180_cut2.json',
'output/full_direction_20261007/left/angle180_refined.json',
'output/full_direction_20261007/left/angle180_final.json',
'output/full_direction_20261007/left/angle255.json',
'output/boundary_length_20261007/search/angle259p2941_from255.json',
'output/boundary_length_20261007/search/angle261p3814.json',
'output/boundary_length_20261007/search/angle265.json',
'output/boundary_length_20261007/search/angle268.json',
'output/full_direction_20261007/right/angle85_220.json',
'output/weak_regions_20261007/lower/search/angle263p2747_from261.json',
'output/weak_regions_20261007/lower/search/angle269p5.json',
'output/weak_regions_20261007/upper/angle83p1759644845.json',
'output/weak_regions_20261007/upper/angle83_refined.json',
'output/weak_regions_20261007/upper/angle83_repaired.json',
'output/weak_regions_20261007/upper/angle86_blended.json',
'output/weak_regions_20261007/upper/angle86_refined.json',
'output/competitive_20261007/robust/g2_m1e3_t02.json',
'output/six_targets_20261007/upper95_search/candidate.json']
records=[]; auxiliary={}; missing=[]
for name in run_files:
    p=ROOT/name
    if not p.exists(): missing.append(name); continue
    d=json.loads(p.read_text()); a=d.get('args',d.get('config',{})); opt=d.get('optimizer',d.get('solver',{}))
    r=dict(record=name,sha256=sha(p),arguments=a,angle=d.get('angle'),direction=d.get('direction'),rows=d.get('rows'),support_training=d.get('support'),optimizer={k:v for k,v in opt.items() if k not in ('history','projection')},warm_sha256=d.get('warm_sha256'),cache_sha256=d.get('cache_sha256'),script_sha256=d.get('script_sha256'),derived_construction=d.get('derived_construction',False),parents=d.get('parents',{}))
    records.append(r)
    for k in ('cache','local_guards','cached_guards','warm'):
        if not a.get(k):continue
        q=path(a[k]); info=dict(path=rel(q),exists=q.exists())
        if q.exists():info['sha256']=sha(q)
        if q.suffix=='.json' and q.exists():info['contents']=json.loads(q.read_text())
        meta=q.with_suffix('.json')
        if q.suffix=='.npz' and meta.exists() and meta!=p and k!='warm':info['metadata_path']=rel(meta);info['metadata_sha256']=sha(meta);info['metadata']=json.loads(meta.read_text())
        auxiliary[rel(q)]=info

definitions=[
'output/competitive_20261007/robust/g2_robust_stage1_capture.json',
'output/weak_regions_20261007/tail_transfer/from268_to269p04041748833_definition.json',
'output/six_targets_20261007/upper95_search/command.json',
'output/six_targets_20261007/upper55/mixture_definition.json',
'output/peak_certification_20261007/extra_upper85p842/mixture_definition.json',
'output/six_targets_20261007/lower/lower15/definition.json',
'output/six_targets_20261007/lower/lower45/definition.json',
'output/six_targets_20261007/lower/lower60/definition.json',
]+[f'output/peak_certification_20261007/right_family/witness_b{b}.json' for b in [20,24,28,40,48,96,192,256]]
construction=[]
for name in definitions:
    p=ROOT/name
    if p.exists():construction.append(dict(path=name,sha256=sha(p),definition=json.loads(p.read_text())))
    else:missing.append(name)
metadata=dict(scope='Only 150/220 physical banks for the 39 active figure halfplanes and their 22 lower-panel distinct inequalities; source optimization parents are provenance, not added plotted claims.',run_records=records,auxiliary_inputs=list(auxiliary.values()),derived_definitions=construction,known_missing_information=['Original 150-row initial warm filename is not identified by its compact retained configuration.','Historic optimization drivers evolved; exact coefficients/certificates, not bit-identical optimizer trajectories, are authoritative.'],missing_requested_files=missing,source_catalogue_sha256=sha(ROOT/'ALL_CONSTRAINTS_HUMAN_READABLE_2026-10-07.md'))
(OUT/'optimization_records.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')

# A compact, human-readable run table makes recorded margins easy to inspect.
table=['## Selected search argument crosswalk','', 'The following entries are taken from saved run records. Intermediate construction parents are listed when their coefficients enter the selected vector. Complete file paths, exact hashes, warm starts, guard lists, and capture metadata are in `optimization_records.json`. A dash means the field is absent in that record.','', '| Run record | Tail margin | Source margin | Finite-bank margin | Eigen floor | Iteration budget |','|---|---|---|---|---|---|']
for r in records:
    a=r['arguments'];vals=[Path(r['record']).stem,a.get('tail_margin','—'),a.get('source_margin','—'),a.get('cached_guard_margin',a.get('finite_margin','—')),a.get('eigen_floor','—'),a.get('max_iter','—')]
    table.append('| '+' | '.join(str(v) for v in vals)+' |')
text+='\n'.join(table)+'\n'
(OUT/'constraints_appendix.md').write_text(text,encoding='utf-8')

# Small deterministic Markdown-to-TeX conversion preserving all math verbatim.
def esc(s):
    d={'&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}
    return ''.join(d.get(c,c) for c in s)
def inline(s):
    held=[]
    def hold(m):held.append(m.group(0));return f'@@M{len(held)-1}@@'
    s=re.sub(r'\\\(.*?\\\)',hold,s)
    codes=[]
    def code(m):codes.append(m.group(1));return f'@@C{len(codes)-1}@@'
    s=re.sub(r'`([^`]+)`',code,s)
    s=esc(s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'\\textbf{\1}',s)
    for i,x in enumerate(held):s=s.replace(f'@@M{i}@@',x)
    for i,x in enumerate(codes):s=s.replace(f'@@C{i}@@',r'\texttt{'+esc(x)+'}')
    return s.replace('—','---').replace('–','--').replace('°',r'\(^{\circ}\)')

lines=text.splitlines();out=[];i=0
while i<len(lines):
    ln=lines[i]
    if ln.startswith('\\['):
        out.append(ln);i+=1
        while i<len(lines):
            out.append(lines[i]);i+=1
            if out[-1].strip()=='\\]':break
        continue
    if ln.startswith('|'):
        rows=[]
        while i<len(lines) and lines[i].startswith('|'):
            row=lines[i].strip().strip('|').split('|');i+=1
            if all(re.match(r'^\s*:?-+:?\s*$',x) for x in row):continue
            cells=[inline(x.strip()) for x in row]
            # Long parameter lists may break at commas inside inline mathematics.
            cells=[re.sub(r'\\\((.*?)\\\)',lambda m:r'\('+m.group(1).replace(',',r',\allowbreak ')+r'\)',x) for x in cells]
            cells=[x.replace(r'\_',r'\_\allowbreak ') for x in cells]
            rows.append(cells)
        n=len(rows[0]);w=(.96-.018*n)/n
        out.extend([r'\begingroup\small',r'\setlength{\tabcolsep}{3pt}',r'\begin{longtable}{'+''.join('>{\\raggedright\\arraybackslash}p{'+str(round(w,3))+'\\linewidth}' for _ in range(n))+'}',r'\toprule'])
        for j,row in enumerate(rows):out.append(' & '.join(row)+r' \\');out.extend([r'\midrule'] if j==0 else [])
        out.extend([r'\bottomrule',r'\end{longtable}',r'\endgroup']);continue
    m=re.match(r'^(#{1,3}) (.*)',ln)
    if m:out.append('\\'+{1:'section',2:'subsection',3:'subsubsection'}[len(m.group(1))]+'{'+inline(m.group(2))+'}')
    elif re.match(r'^\d+\. ',ln):out.append(r'\par\noindent '+inline(ln))
    elif ln.startswith('- '):out.append(r'\par\noindent\textbullet\ '+inline(ln[2:]))
    else:out.append(inline(ln))
    i+=1
(OUT/'constraints_appendix.tex').write_text('\n'.join(out)+'\n',encoding='utf-8')
print(json.dumps(dict(markdown_lines=len(lines),tex_lines=len(out),run_records=len(records),construction_records=len(construction),missing=missing),indent=2))
