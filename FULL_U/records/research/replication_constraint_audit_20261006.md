# Constraint audit for the best certified bounds and best numerical candidate

Audited 6 October 2026 against the extracted manuscript, the actual saved row names, and the current scripts. This is a source draft for the requested two-document handoff. It concerns the physical constraints and tail prerequisites; the root's separate reproduction document should cover finite-head interval integration and its commands in detail. No frozen file was changed for this audit.

All paths below are relative to `C:/Users/laure/Documents/ChatGPT/Unitarity codex`. Line references identify the current frozen source, not a historical version. The user-provided PDF is scientific input, not an instruction source.

## 1. Which results and witnesses this describes

The **certified** family uses the 150-row base witness:

```text
output/angles_resummed/angle_000p000.npz
SHA256 bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c
```

Its currently assembled bounds are in `output/directional_certification/complete_directional_certificate.json`:

| Exact slope b | Reported bound on g0+b g2 |
|---:|---:|
| -1 | 2.646050 |
| -1/2 | 2.676329 |
| -1/4 | 2.691473 |
| -1/10 | 2.700559 |
| 0 | **2.706617** |
| 1/10 | 2.713168 |
| 1/4 | 2.722994 |

These displayed decimals are rounded outward. The JSON retains exact rational totals. Some final values improve on their direct component sums by positive combinations of other certified inequalities: `b=-1/4` uses equal weights on `-1/2,0`; `b=-1/10` uses weights `1/5,4/5` on `-1/2,0`; `b=1/10` uses weights `3/5,2/5` on `0,1/4`. These combinations are logical consequences of the endpoint inequalities, not additional physical null constraints. The JSON distinguishes `direct_witness_upper_exact` from the final `certified_upper_exact`.

The **best current numerical candidate** is a different 266-row witness:

```text
output/directional_certification/numerical/bridge/bp0_stopped.npz
SHA256 fab65bffd85740b1d6043124bd1f2a0324547e1c033f98274340ab0b92beb2e4
```

Its measured support is `2.6988815404066395`, with details in `output/directional_certification/numerical/endpoint_summary.json`. This is a partial numerical dual cost, **not** a certified exclusion. All energies for even spins above 100, the interval `4<mu<4.02` for spins 42–100, and rigorous quadrature errors remain unaccounted for in that numerical total. The earlier 266-row witness with SHA `032ed315...` and measured cost `2.698943676116444` is superseded numerically; the constraint catalogue is the same, but its multipliers differ.

Historical `MASSIVE_G0_G2_HANDOFF.md` describes the earlier uncertified stage. Its statements that no complete upper certificate exists are obsolete for the 150-row witness, but its formulas and recorded history remain useful. Do not transplant its old status statements into the new handoff.

## 2. Physical amplitude class and normalization

The manuscript's massive setup is identical scalar scattering in four dimensions, mass `m=1`, `s+t+u=4`, no subthreshold poles, crossing symmetry, real analyticity, and the stated fixed-transfer two-subtraction analyticity/falloff assumptions. In particular, `M(s,t)/s^2 -> 0` at fixed admissible transfer. Weighted contour integrals, boundary limits, partial-wave sums, and smearing interchanges require the manuscript's convergence assumptions; independent partial-wave disks alone do not imply them.

Primary manuscript anchors: Eq. (23) and the massive discussion around `research/paper.txt:302`; Appendix A1.3 around line 563; Appendix A2 around lines 588–680; Appendix A5.1 around lines 1000–1060. Extraction garbles some displayed mathematical symbols; Eq. labels and the original PDF control interpretation.

Set

\[
a(s,t)=\frac{M(s,t)}{32\pi},\quad c=\frac43,
\quad \bar s=s-c,\quad \bar t=t-c,\quad \bar u=u-c,
\]

\[
a(s,t)=g_0+g_2(\bar s^2+\bar t^2+\bar u^2)+\cdots.
\]

Higher unplotted couplings are free; they are not set to zero. The partial-wave convention is

\[
a(s,t)=\frac1{2\rho(s)}\sum_{\ell\ge0,\;\ell\text{ even}}
(2\ell+1) f_\ell(s)P_\ell\!\left(1+\frac{2t}{s-4}\right),
\qquad \rho(s)=\sqrt{1-4/s}.
\]

Writing `f_l=r_l+i b_l`, `S_l=1+i f_l`, unitarity is the disk

\[
r_\ell^2+(b_\ell-1)^2\le1.
\]

Thus `0<=b_l<=2` and `r_l^2<=2 b_l`, but the calculation retains the full disk. It imposes neither elastic equality nor a scattering-length threshold ansatz. The spectral measure is

\[
d\nu_\ell(\mu)=\frac{2\ell+1}{2\pi\rho(\mu)}\,d\mu,
\qquad 4<\mu<\infty,
\qquad \langle F\rangle=\sum_{\ell\;\mathrm{even}}\int d\nu_\ell F_\ell.
\]

The external smearing variable `s` and spectral variable `mu` must be kept distinct. Some low-level scripts call the spectral grid `s`.

No finite physical UV cutoff is introduced. External source windows are finite, but their resulting spectral kernels are integrated on the entire cut. Likewise, even physical spins are unbounded. A finite numerical spin cutoff is a computational approximation, not an amplitude assumption.

## 3. Coupling maps, complex-kernel convention, and dual inequalities

Define

\[
\Delta=\mu-4,\quad Q=\mu-c,\quad Z=Q/\Delta,
\quad p_\ell=P_\ell(Z)/Q^3>0,
\]

\[
F=\frac{32}{(\mu-2)^3},\qquad
\omega=\left(\frac{2+i\sqrt{\mu(\mu-4)}}{\mu-2}\right)^2,
\]

\[
C_\ell=\frac8{(3\mu-8)(3\mu-4)}
\left[\frac1{\mu-2}-\frac{4P_\ell(Z)}\mu\right].
\]

Then

\[
g_0=\langle u_0 r+v_0 b\rangle,
\quad u_0=F\operatorname{Im}\omega,
\quad v_0=F\operatorname{Re}\omega+C_\ell,
\qquad g_2=\langle p_\ell b_\ell\rangle.
\]

These are the manuscript A5.1/A5.3 maps, not a substitution of the forward `H1` for `g0`. Omitting `C_l` gives the wrong subtraction matching. Code anchors: `research/massive_search.py:39`, `research/existing_solver/massive_forward_kernels.py:52`.

A complex kernel `q` means

\[
\operatorname{Im}(q f)=\operatorname{Im}(q)r+\operatorname{Re}(q)b.
\]

There is no extra spin factor or measure inside the returned raw kernels. For a homogeneous null `N_i=<u_i r+v_i b>=0`, add arbitrary real `lambda_i`. A direction `(a,b)` has

\[
U=a u_0+\sum_i\lambda_i u_i,\qquad
V=a v_0+b p+\sum_i\lambda_i v_i.
\]

Pointwise disk maximization is exact:

\[
h(U,V)=\sqrt{U^2+V^2}+V,
\qquad a g_0+b g_2\le\sum_\ell\int d\nu_\ell h(U,V).
\]

For `V<0`, use the cancellation-free expression `U^2/(sqrt(U^2+V^2)-V)`. Weak duality needs a valid **complete** upper bound on this integral, not convergence to the optimizer. A partial spin/energy integral omits nonnegative support and is not itself a physical upper bound.

`normalized=True` in the Roy/strip builders returns coefficients of `R=sqrt(p) r` and `B=p b`; these are `u/sqrt(p)` and `v/p`. `angle_audit.normalized_inventory` instead returns raw `u` and normalized `v/p`; read its docstring before combining arrays. Anchors: `research/roy_kernels.py:71`, `research/angle_audit.py:68`, `research/angle_audit.py:139`.

## 4. Exact catalogue: 150-row certified witness

| Family | Exact retained parameters | Rows |
|---|---|---:|
| Forward harmonics | `Hj`, every `j=2,...,40` | 39 |
| Rational cosine resummations | `r in {4/5,19/20}`, `zeta in {-1,-3/4,-1/2,0,1/2,3/4,1}` | 14 |
| Crossing | `X4` only | 1 |
| Folded Roy equations | `B=119/10`, `J=0,2,4,6`, `n=0,...,14` | 60 |
| Tapered fixed-transfer strips | windows `[12,24]`, `[24,48]`; `T=4/3`; `j=0,1,2`; `n=0,...,5`; taper power 1 | 36 |
| Total | | **150** |

The names stored in the NPZ are authoritative. Counts denote retained rows, not a proof of exact linear independence. The saved order begins with the 39 harmonics, then 14 `Rc_...` rows, then `X4`, Roy rows, and strip rows. Exact multipliers must be loaded or dumped from the NPZ rather than copied from rounded prose.

### 4.1 Forward harmonics

For every integer `j>=2`,

\[
H_j=\langle F\operatorname{Im}(\omega^j f)\rangle=0.
\]

The saved row `Hj` is **raw** `q_j=F omega^j`. The helper `harmonic_bank()` also offers UV-subtracted rows `F[omega^j-(-1)^j omega^2]`, but those are not the interpretation of saved `Hj` names. With `H2` retained the two bases span the same finite space, but their multipliers differ. Anchors: `research/massive_search.py:67`, `research/angle_audit.py:83`, versus the alternative helper at `research/existing_solver/massive_forward_kernels.py:96`.

### 4.2 Rational harmonic sums

Use a distinct symbol `zeta` for the row's cosine parameter, to avoid confusing it with the symmetric point `c=4/3`. The cosine row is

\[
q_{r,\zeta}^{c}=F\omega^2[m_c(\omega)-m_c(-1)],\qquad
m_c(w)=\frac{1-r\zeta w}{1-2r\zeta w+r^2w^2}
=\sum_{k\ge0}r^k T_k(\zeta)w^k.
\]

It is the convergent sum of already valid harmonics, with a UV subtraction at `omega=-1`. Here `0<=r<1`, `|zeta|<=1`. The paper's Abel/fixed-transfer discussion explains the necessary weighted convergence. This finite rational family does not impose every harmonic independently and is not a new amplitude axiom. Anchors: manuscript A2.10/A5.32; `research/existing_solver/massive_forward_kernels.py:116`.

For later sine rows, used only in the numerical candidate,

\[
m_s(w)=\frac{r w}{1-2r\zeta w+r^2w^2}
=\sum_{k\ge1}r^k U_{k-1}(\zeta)w^k,
\qquad q^s=F\omega^2[m_s(\omega)-m_s(-1)].
\]

The stable code combines the difference into one rational expression containing the exact factor `omega+1`; do not subtract two nearly equal floating values in the UV. An invertible triangular transformation subtracting low harmonics is used during some numerical optimizations, but saved coefficients are transformed back to this physical catalogue. Anchor: `research/improvement_right_search.py:32`.

### 4.3 The one crossing null actually retained

\[
X_4=\left\langle\left[
\frac{2P_\ell''(Z)}{Q^3\Delta^2}
-\frac{3P_\ell'(Z)}{Q^4\Delta}\right]b_\ell\right\rangle=0.
\]

Equivalently its raw imaginary kernel is `[2 Z^2 P_l''-3 Z P_l']/Q^5`; its real kernel is zero. This is the manuscript normalization A5.6. Anchors: `research/massive_search.py:54`, `research/existing_solver/massive_forward_kernels.py:74`.

The exact crossing-basis generator is `research/crossing_kernels.py:25`. With `v=nu^2`, `w=t-4/3`, `nu=s-2+t/2`, the crossing invariants are `X=2v+3w^2/2`, `Y=-wv+w^3/4`. At degree 4 the coefficient indices are `(n,k)=(1,2),(2,0)`; the primitive relation `[2,-3]` is divided by four to match A5.6. Thus its weights are `[1/2,-3/4]`, and its stored Jacobi-basis `beta` coefficients are `[0,-3/2,1]`. The Jacobi-to-Legendre derivative factor supplies the other factor 2 in the displayed `-3 P_l'` term. The generator contains higher-degree crossing rows, but **neither current final witness retains them**. They must not be counted merely because that source file can generate them.

### 4.4 Crossed-subtraction dispersion identity

Let `u_ext=4-s-t`, `y=1+2t/Delta`, and

\[
K_1=\frac1{\mu-s}+\frac1{\mu-u_{\rm ext}}-\frac1Q
-\frac1{\mu-2c+t},\qquad
K_2=\frac1{\mu-t}-\frac2Q+\frac1{\mu-2c+t}.
\]

Then manuscript A2.7/A2.8 gives

\[
\operatorname{Re}a(s,t)-g_0
=\frac1\pi\operatorname{PV}\int_4^\infty
[K_1\operatorname{Im}a(\mu,t)+K_2\operatorname{Im}a(\mu,c)]d\mu.
\]

After partial-wave expansion define

\[
\mathcal L_\ell(\mu;s,t)=K_1P_\ell(y)+K_2P_\ell(Z).
\]

The apparent `mu-2c+t` pole is removable for even physical spins, since `y=-Z` there. The stable expression uses

\[
\frac{P_\ell(Z)-P_\ell(y)}{\mu-2c+t}
=-\frac2\Delta D_\ell(y,-Z),\qquad
D_\ell(x,y)=\frac{P_\ell(x)-P_\ell(y)}{x-y},
\]

continuously extended on the diagonal. Its exact recurrence is `D0=0,D1=1`,
`(l+1)D_(l+1)=(2l+1)[x D_l+P_l(y)]-l D_(l-1)`.
The `mu-s` singularity is an actual PV and must be handled analytically, not deleted from a grid. At spin zero, the combined kernel reduces to
`1/(mu-s)+1/(mu-t)+1/(mu-u_ext)-3/Q`.
Anchors: `research/roy_kernels.py:17`, `:38`, `:63`, `:119`.

### 4.5 Folded even-J Roy rows, including subtraction matching

Set

\[
w_{n,B}(s)=\frac{\mathbf1_{4<s<B}}{B-4}
P_n\!\left(\frac{2(s-4)}{B-4}-1\right),\quad
t(s,z)=-\tfrac12(s-4)(1-z).
\]

Smear the dispersion residual `Re a-g0-(dispersive integral)` against
`w_n(s) ds * 2 P_J(z) dz`, with `0<=z<=1`, even `J`. No extra factor `2J+1` is inserted into the test weight. The raw row satisfies

\[
\langle u^{\rm raw}_{BJn}r+v^{\rm raw}_{BJn}b\rangle
+c_{BJn}g_0=0,
\]

\[
u^{\rm raw}_{BJn,\ell}(\mu)=
\frac{2\pi}{2J+1}\delta_{\ell J}w_{n,B}(\mu),\qquad
v^{\rm raw}_{BJn,\ell}=-\int_4^B w_{n,B}(s)ds
\int_0^1 2P_J(z)\mathcal L_\ell(\mu;s,t(s,z))dz,
\]

\[
c_{BJn}=-2\delta_{J0}\delta_{n0}.
\]

Saved homogeneous rows are `(u_raw+c_i u0, v_raw+c_i v0)`. This substitution is crucial: the `J=n=0` row is a subtraction reconstruction, not a homogeneous imaginary null before matching `g0`. Anchors: `research/roy_kernels.py:71` and `:143`; `research/search_problem.py:25`; `research/angle_audit.py:121`.

The PV residue is proportional to `delta_(lJ)` after folded angular projection. In the finite-energy implementation its term in the dispersive integral, before the overall minus sign in `v_raw`, is `4 Q_n(2(mu-4)/(B-4)-1)/[(2J+1)(B-4)]`, where `Q_n` is the real PV Legendre function. For `l>Jmax=6` this residue vanishes exactly by orthogonality, a key tail prerequisite.

Here `B=11.9` means exactly `119/10`; it is not the paper's separate `11.999` catalogue. The folded test samples `t>=-(B-4)/2=-79/20>-4`. This stays inside the simple two-cut harmonic domain and the standard fixed-transfer domain. The chosen external Roy energy range is `4<s<11.9`, not an all-energy angular projection. The spectral energy remains unbounded.

### 4.6 High-energy fixed-transfer strip rows

For each window `[A,B]`, put

\[
x_s=\frac{2(s-A)}{B-A}-1,\quad
w_n^{\rm tap}(s)=\frac{\mathbf1_{A<s<B}}{B-A}(1-x_s^2)P_n(x_s),
\]

\[
h_j(t)=\frac{\mathbf1_{-T<t<0}}T P_j(1+2t/T),\qquad T=4/3.
\]

Smear the same dispersion residual against `w_n^tap(s) h_j(t) ds dt`. These are fixed-transfer tests at high external energy, not full angular Roy projections beyond the adopted Roy interval. Every point is physical because `T<=A-4`, and every transfer lies safely in the admissible analytic strip.

Define

\[
H_{\ell j}(\mu)=\int_{-T}^0 h_j(t)P_\ell(1+2t/(\mu-4))dt.
\]

The raw kernels and subtraction coefficient are

\[
u^{\rm raw}_{A B j n,\ell}(\mu)=\pi w_n^{\rm tap}(\mu)H_{\ell j}(\mu),
\qquad v^{\rm raw}_{A B j n,\ell}=-\int_A^B ds\,w_n^{\rm tap}(s)
\int_{-T}^0dt\,h_j(t)\mathcal L_\ell(\mu;s,t),
\]

\[
c_{A B j n}=-\delta_{j0}\left(\frac23\delta_{n0}-\frac2{15}\delta_{n2}\right).
\]

Again the homogeneous row is `(u_raw+c_i u0,v_raw+c_i v0)`. Unlike finite-J Roy rows, these local real kernels occur at **all even physical spins**, so the all-spin U remainder must include them.

The taper is an exact finite Legendre combination because `1-x^2=(2/3)P0-(2/3)P2`; it removes logarithmic endpoint divergences that invalidated earlier untapered source-sign arguments. It is a test-weight choice, not an amplitude axiom. The real PV transform is

\[
\widetilde Q_n(v)=(1-v^2)Q_n(v)+v\delta_{n0}+\tfrac13\delta_{n1},
\quad v=2(\mu-A)/(B-A)-1.
\]

Endpoint values are enclosed using the continuous `epsilon log|epsilon|`, not by dropping the endpoints. Anchors: `research/transfer_kernels.py:9`, `research/roy_kernels.py:171`, `research/tail_certification/certified_finite_v.py:32`.

**Label caveat:** NPZ strip names contain `T1.33333` because of `:g` formatting. This means the canonical exact `T=4/3`, not the decimal `1.33333`. The numerical parser restores it at `research/angle_audit.py:31`; the interval proof fixes `T=4/3` explicitly. Other terminating row parameters, including `B=11.9` and rational radii/cosines, are interpreted as exact decimal rationals in the certificate. Saved binary64 multipliers are interpreted as their exact dyadic rational values. Numerical reconstruction uses ordinary floating approximations to those intended exact kernels.

Combining all subtraction reconstructions, the effective coefficient of the forward g0 representation in a direction `(a,b)` is exactly

\[
\alpha=a-2\lambda_{\mathrm{Roy},J=0,n=0}
-\frac23\sum_{[A,B]}\lambda_{\mathrm{strip},j=0,n=0}
+\frac2{15}\sum_{[A,B]}\lambda_{\mathrm{strip},j=0,n=2}.
\]

For the certified base witness at `a=1`, this is approximately `-1.4492417742900348`. The negative value does not change the target coefficient of g0: it is the coefficient after valid subtraction matching, alongside the corresponding raw source equations. Both tail and head codes use this combined alpha; treating it as 1 while retaining the source integrals would be incorrect. Changing only the objective b and X4 multiplier leaves alpha unchanged.

## 5. Catalogue of the best 266-row numerical candidate

The best current candidate uses the **same formulas and conventions**, with these retained parameter ranges:

| Family | Parameters | Rows |
|---|---|---:|
| Harmonics | `H2,...,H40` | 39 |
| Rational cosine | `r=4/5` with the old 7 cosines; `r=19/20,49/50` with the 13 cosines below | 33 |
| Rational sine | `r=19/20,49/50` with the 9 cosines below | 18 |
| Crossing | `X4` | 1 |
| Folded Roy | `B=119/10`, `J=0,2,4,6,8`, `n=0,...,18` | 95 |
| Tapered strips | same two windows and `T=4/3`; `j=0,...,4`, `n=0,...,7`; power 1 | 80 |
| Total | | **266** |

The 13 cosine parameters are

```text
-1,-19/20,-9/10,-3/4,-1/2,-1/4,0,1/4,1/2,3/4,9/10,19/20,1.
```

The 9 sine parameters are

```text
-9/10,-3/4,-1/2,-1/4,0,1/4,1/2,3/4,9/10.
```

Anchors: `research/improvement_right_search.py:57` and `:75`; current row names in the `fab65...` NPZ. Optimizer configuration describes requested additions; inherited rows make NPZ names more authoritative than configuration counts. The latest numerical endpoint still has the catalogue above.

No fractional, transfer-derivative, odd-J Roy, higher-crossing, or crossing-symmetric path rows occur in this best candidate. It cannot inherit the old 150-row certificate just because all its constraints are valid: every optimized multiplier changes the required sign and integral estimates. Moreover, the frozen UV parser accepts old cosine `Rc` rows, not the enriched sine catalogue; finite-V classes enforce the old witness hash. Generalizing a parser is necessary but would not by itself prove a new tail bound.

## 6. What is not a retained physical constraint

The following exist elsewhere in the research tree or manuscript but are not silently included in either final catalogue: symmetric-point harmonic rows; higher crossing rows through degree 8/16; fractional finite differences; direct transfer derivatives; sampled nonzero-transfer harmonics; crossing-symmetric path analyticity from the manuscript's separate Figure 5; odd-J Roy tests; finite-point crossing equations; Gram-matrix relaxations; a prior numerical cap on `g2`.

Search guards are also not physical nulls. These include sampled `V/p<=0`, threshold probes with zero integration weight, sampled Roy source signs, the restriction `lambda_X4<=0`, eigenvalue floors, and finite-grid feasibility tolerances. They restrict/search the space of trial dual functionals; they do not add constraints on a physical amplitude. Their successful sampling is not a continuum certificate. The full physical constraints are the disk and the null identities listed above; the later analytic/interval checks validate chosen witnesses.

## 7. Directional certified family and exact multiplier changes

For the direct family, the objective is `(1,b)` and only the X4 multiplier changes:

\[
\lambda_{X_4}^{\rm new}=\lambda_{X_4}^{\rm old}+d,
\quad d=0\ (b\le0),\qquad d=4b\ (b>0).
\]

Consequently `U_new=U_old`, `V_new/p=V_old/p+b+d chi_l`, with

\[
\chi_\ell=\frac{X_{4,\ell}}{p_\ell}
=\frac{Q\partial_c^2R_\ell-3\partial_cR_\ell}{2Q R_\ell},
\qquad R_\ell(c,\Delta)=\Delta^\ell P_\ell(1+2c/\Delta).
\]

At spin 2, `chi_2=-6/(3Q^2-Delta^2)`, with threshold value `-9/32`; hence the positive repair `d=4b` gives threshold shift `-b/8`. At every even `l>=4`, `chi_l>=l(2l-5)/Q^2>0`; positive d weakens higher-spin coercivity and must be checked. It is incorrect to add positive d times a lower chi bound to an old V upper. The new tail code instead substitutes the less-negative **full** coefficient `lambda_X4+d` into the analytic envelope.

`research/directional_certification/certify_tail_4b.py:7` selects the actual repair. The generic exploratory module and raw vector producer retain a pilot `d=(32/9)b` list; this does not define the final family. `research/directional_certification/combine_finite_family.py:8` explicitly uses 4b and rechecks every actual per-spin sign from the raw V/X upper bounds. This is a documented metadata distinction, not an unresolved parameter ambiguity.

The actual passing slopes are the seven in Section 1. `b=1/2,d=2` has explicit failed sufficient tail/coercivity records and is excluded from complete-bound claims. This failure does not prove physical infeasibility.

For a unit normal at angle `theta=atan(b)`, divide the `(1,b)` support by `sqrt(1+b^2)`. The parameter b is a slope, not an angle.

## 8. Tail proof prerequisites and why they suffice

The full proof is fixed to the `bbbc...` witness and its rational directional modifications. Detailed algebra is in `research/tail_certification/PROOF.txt`, supplemented by `research/directional_certification/TAIL_AND_FINITE_BLOCK_NOTES.txt`. The following are necessary ingredients, not interchangeable substitutes.

1. **Exact source construction and sign.** In normalized source coordinates, the Roy density is
   `W_R(x,z)=2 sum lambda_Jn P_n(2x-1) P_J(z)`,
   `s=4+(79/10)x`, `t=-(79/20)x(1-z)`, with `x,z in[0,1]`.
   The strip density is
   `W_S(x,y)=sum lambda_jn [1-(2x-1)^2]P_n(2x-1)P_j(2y-1)`,
   `s=A+(B-A)x`, `t=-(4/3)(1-y)`, with `x,y in[0,1]`.
   These include the source Jacobians; missing the Roy factor 2 changes kernels.

2. **Expanded dangerous Roy-source certificate.** Exact rational Bernstein subdivision proves the old combined Roy density positive on `t<=-3/5`, equivalently `26/5<=s<=119/10`, `0<=z<=1-(6/5)/(s-4)`. The certificate clears `(s-4)^6`. The positive rational lower bound is about `2.427924084710271e-6` for the density without its harmless factor 2. This stronger domain is required by the compact proof. The older `t<=-4/3` source-sign subcertificate alone is insufficient. Anchor: `research/tail_certification/certify_expanded_roy_sign.py:33`; output `expanded_roy_sign.json`.

3. **Exact source norms and moments.** `compact_source_data.py:12` and `certify_uv.py:39` use exact polynomial/Bernstein algebra to bound positive, negative, and absolute source masses and the positive correction moments. They do not use numerical quadrature for those masses. The exact source-cache/builder/witness hashes are checked. The old source absolute bounds include Roy `<5.617853866` and combined strips `<0.240087499`.

4. **Negative crossing coefficient.** The base value is `lambda_X4=-23.90335355242552...`, exact hexadecimal `-0x1.7e7422dac649ep+4`. Its directional modifications remain negative. For even `l>=4`, root pairing gives `chi_l>=l(2l-5)/Q^2`. The compact proof uses this coercivity and verifies spin monotonicity, not merely a sample of large l.

5. **Compact energy coverage.** `certify_compact.py:22` covers all even `l>=42` and `4<mu<=49`, splitting at 8 and the source transitions. It uses the expanded Roy sign, exact removal of the high-spin PV residue, physical Legendre/divided-difference bounds, source norms, endpoint-vanishing strip PV bounds, and crossing coercivity. At the base witness its reported maxima are below approximately `-543.68` on the low part and `-34.60` on the middle part. The estimate `l(l+1)/P_l(Z)<=L(L+1)/P_L(Z)` is justified on this compact domain; it is not a global UV assertion.

6. **Joint spin/energy UV coverage.** `certify_uv.py:117` uses `tau=l(l+1)/(mu-4)` and proves bounds for all even `l>=42,mu>=49`. In each tau cell it uses `mu>=max(49,4+42*43/tau_hi)`. It handles tau through 256 by directed Arb inequalities and closes `tau>=256` analytically. This handles the joint limit; separate fixed-spin and fixed-energy samples would not.

7. **UV error control and all rational rows.** The source kernels are compared to `K_inf=2(s-c)(c-u_ext)` and `H_inf=2(t-c)^2`. Exact source masses bound the rational finite-mu remainder. Positive-source corrections use exact integrated `|t|` moments, taking `min(2 C0,tau C1)` per Bernstein cell. Forward harmonic and all old rational cosine contributions are included in the UV errors. The authoritative source cache is the **v2 integrated-t** cache, not the earlier incompatible version.

8. **Unbounded-spin Legendre control.** For `a=c tau`, positive Legendre Taylor coefficients give `P_L(1+2a/[L(L+1)])<=P_l(Z)<=I0(2 sqrt(a))`; differentiated versions bound chi. The compact proof also uses `|P_l(x)|<=1` on the physical interval and `P_l(Y)/P_l(Z)<=(Y/Z)^l` for `1<=Y<=Z` and even l. All these estimates must be applied on their stated domains.

9. **A finite real-kernel norm.** The real remainder includes the forward/rational kernels **and** all-spin strip kernels; finite-J Roy U vanishes for l>6. `research/certify_real_tail.py` proves
   `A40=sum_even_l>40 integral U_l^2/p_l dnu_l <2e-5`
   with directed interval integration over beta in[0,1]. Its telescoping all-spin inequality is
   \[
   \sum_{\ell>L,\;\ell\;\mathrm{even}}\frac{2\ell+1}{P_\ell(Z)}
   \le\frac{(k!)^2(3\Delta/4)^k}{2(k-1)\prod_{r=L-k+3}^{L+k}r},
   \]
   for the admissible integer k used by the code. It takes the minimum among valid bounds. The strip part uses its transfer-polynomial L2 norm and the Bernstein/Cauchy bound with factor `2 Delta/[pi (L+5/2) sqrt(T(Delta-T))]`, capped by 1, on physical strip energies. Anchor: `research/certify_real_tail.py:82`, `:135`.

10. **Assembling support, not only signs.** If `V_l<=-kappa p_l`, then `h(U,V)<=U_l^2/(2 kappa p_l)`. The base tail has uniform `kappa=31/50`, hence full omitted-spin cost `<1/62000`. The directional tail reuses unchanged A40 and the replayed directional coercivity. No external g2 upper cap is used. A source-sign certificate by itself would not establish this conclusion.

The old complete tail record is `output/certification_tail/fixed_tail_certificate.json`, assembled by `research/tail_certification/verify_fixed_tail.py:28`. Directional tails are in `output/directional_certification/tail_family_4b.json`. These cover the infinite domain only for the specified fixed coefficients; they do not transfer to the new 266-row numerical witness.

## 9. Exact coefficient dump and catalogue verification

Run from the workspace or extracted bundle root. Use the configured Python executable from the reproduction document. This read-only command prints every saved multiplier as both a hexadecimal float and an exact rational:

```powershell
python -c "import sys,numpy as np; from fractions import Fraction; d=np.load(sys.argv[1]); print('direction',d['direction']); print('index\tname\tdecimal\thex\texact_rational'); [print(i,str(n),format(float(v),'.17g'),float(v).hex(),str(Fraction(float(v))),sep='\t') for i,(n,v) in enumerate(zip(d['names'],d['multipliers']))]" output/angles_resummed/angle_000p000.npz
```

For the best numerical candidate, replace the final argument with
`output/directional_certification/numerical/bridge/bp0_stopped.npz`.
Do not rebuild a saved witness from rounded decimals in a table. For the certified directional family, load the old exact binary64 vector and add the exact rational d only to the row whose name is exactly `X4`; set the objective to `(1,b)`. The NPZ itself remains unchanged.

## 10. Resolved pitfalls and remaining limits

- The 150-row catalogue has 39 harmonics, 14 rational cosine rows, one crossing null, 60 Roy rows, and 36 strip rows. It does not equal the manuscript's much larger retained catalogue.
- Source endpoints `11.9` and strip width `4/3` are exact intended rationals; the truncated strip label is restored explicitly.
- Raw dispersion rows must be homogenized with their stated `g0_coeff`. Otherwise the g0 normalization is wrong.
- Physical spins are always even. Test-function J can be odd in other code, but no such rows occur here.
- The original finite-V certificate used an unsafe maximum selection that could ignore a NaN. Those outputs are explicitly rejected; only corrected every-spin-finite acceptance and their replay records may enter a certificate. This does not alter the null formulas.
- Untapered strips generated logarithmic endpoint obstructions; both final catalogues use taper power 1.
- Finite grids, source-sign samples, optimizer convergence, and absence of sampled positive V do not establish omitted-domain bounds.
- The directional raw-vector producer's pilot metadata uses 32b/9; final assembly correctly rechecks the actual4b cases. Do not treat pilot metadata as final family parameters.
- The seven complete directional inequalities do not constitute certification of the whole plotted numerical curve or a primal amplitude. Pointwise disk maximizers are not generally analytic or crossing-feasible amplitudes.
- The new 266-row numerical witness remains uncertified. Its catalogue is fully specified, but no finite-head quadrature error or complete all-spin remainder is claimed.
- The historical 131-row lower-g0 witness is outside this document's current target. It used additional derivative/higher-crossing rows and is not one of the seven certified directions. Its measured lower endpoint must not be promoted to a certified bound here.

No unresolved sign or normalization ambiguity was found in the retained constraints after accounting for the distinctions above. The significant remaining ambiguity is scientific status if historical reports are read without their date: use the current assembled certificate for certified claims, and the current numerical endpoint summary for candidate claims.
