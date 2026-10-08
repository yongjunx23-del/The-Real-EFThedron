# Constraints and construction of the dual witnesses

This appendix defines the constraints used by the certified boundary in the final two-panel figure. The 39 active half-planes of the upper panel include the 22 distinct inequalities represented in the lower panel (the minimum-\(g_0\) inequality is displayed at both ends of the perimeter). Each half-plane is an exclusion inequality, not a constructed scattering amplitude. The accompanying registry ties every inequality and plotted position to its exact normal, multiplier vector, and selected proof packet. The smooth connecting segments in the lower panel are display guides and introduce no new certified inequalities.

There are two physical row banks, of sizes 150 and 220. The analytic inequality \(g_2\ge0\) is also included. The counts specify retained identities, not linear independence or a complete characterization of physical amplitudes. No independent Gram or SOS coupling constraint is used. The rational source-polynomial inequalities below constrain a dual witness during its search; they do not impose new assumptions on amplitudes.

The certificates use the frozen coefficients exactly. A binary64 coefficient denotes its exact dyadic rational, while a rational combination of witnesses retains its exact rational weights. Nominal angles identify records only: the rational normal supplied with the inequality is authoritative. Multipliers in the data file are expressed in the physical rows defined below, after undoing every numerical conditioning transformation.

## Physical assumptions and conventions

The amplitude is identical-scalar scattering in four dimensions with mass squared one, \(s+t+u=4\), crossing symmetry, real analyticity, no subthreshold poles, and the manuscript's fixed-transfer analyticity and twice-subtracted dispersion assumptions. In particular, the assumed fixed-transfer falloff is \(M(s,t)/s^2\to0\). Contour deformation, boundary limits, partial-wave expansion, differentiation, and smearing require the convergence conditions of the underlying identities. These are assumptions of the amplitude class; they do not follow merely from partial-wave unitarity.

Write

\[
\mathcal A(s,t)=\frac{M(s,t)}{32\pi},\qquad c=\frac43,
\qquad \bar s=s-c,\quad\bar t=t-c,\quad\bar u=u-c,
\]

\[
\mathcal A=g_0+g_2(\bar s^2+\bar t^2+\bar u^2)+\cdots.
\]

Higher EFT coefficients remain free. They are not set to zero. The partial-wave normalization is

\[
\mathcal A(s,t)=\frac1{2\rho(s)}
\sum_{\ell=0,2,4,\ldots}(2\ell+1) f_\ell(s)
P_\ell\!\left(1+\frac{2t}{s-4}\right),\qquad
\rho(s)=\sqrt{1-\frac4s}.
\]

Let \(f_\ell=r_\ell+i b_\ell\) and \(S_\ell=1+i f_\ell\). At every \(\mu>4\) and every nonnegative even spin the unitarity constraint is the full disk

\[
\boxed{r_\ell(\mu)^2+(b_\ell(\mu)-1)^2\le1.}
\]

This includes \(0\le b_\ell\le2\) and \(r_\ell^2\le2b_\ell\), but the full disk is used in the dual support. Elastic equality is not imposed. Neither a finite physical spin cutoff nor a finite physical energy cutoff is assumed.

Define the common spectral measure and bracket by

\[
d\nu_\ell(\mu)=\frac{2\ell+1}{2\pi\rho(\mu)}d\mu,
\qquad
\langle F\rangle=\sum_{\ell=0,2,\ldots}\int_4^\infty d\nu_\ell(\mu)F_\ell(\mu).
\]

The spectral variable \(\mu\) is distinct from the external source variable \(s\) used in Roy and strip smearing. Returned kernel pairs contain no extra measure or \(2\ell+1\) factor.

## Coupling identities and the dual inequality

Put

\[
\Delta=\mu-4,\qquad Q=\mu-c,\qquad Z=Q/\Delta,
\qquad p_\ell(\mu)=\frac{P_\ell(Z)}{Q^3}>0,
\]

\[
F(\mu)=\frac{32}{(\mu-2)^3},\qquad
\omega(\mu)=\left(\frac{2+i\sqrt{\mu(\mu-4)}}{\mu-2}\right)^2,
\]

\[
C_\ell(\mu)=\frac8{(3\mu-8)(3\mu-4)}
\left[\frac1{\mu-2}-\frac{4P_\ell(Z)}\mu\right].
\]

The actual coupling maps used are

\[
\boxed{g_0=\langle u_0r+v_0b\rangle,
\quad u_0=F\operatorname{Im}\omega,
\quad v_0=F\operatorname{Re}\omega+C_\ell,}
\]

\[
\boxed{g_2=\langle p_\ell b_\ell\rangle.}
\]

Thus \(g_2\ge0\) follows analytically and supplies the exact bottom cut. The \(C_\ell\) term in \(g_0\) is essential. Replacing the complete \(g_0\) map by the forward first harmonic changes the problem.

Our complex-kernel convention is

\[
\operatorname{Im}(q f)=\operatorname{Im}(q)r+\operatorname{Re}(q)b.
\]

A homogeneous identity \(N_i=\langle u_i r+v_i b\rangle=0\) is multiplied by a free real \(\lambda_i\). Define

\[
U=a u_0+\sum_i\lambda_i u_i,\qquad
V=a v_0+b p+\sum_i\lambda_i v_i.
\]

The disk implies the pointwise, exact support inequality

\[
Ur+Vb\le H(U,V):=\sqrt{U^2+V^2}+V,
\]

and therefore

\[
\boxed{a g_0+b g_2\le\langle H(U,V)\rangle.}
\]

When \(V<0\), evaluating \(H=U^2/(\sqrt{U^2+V^2}-V)\) avoids cancellation. For a proof region where \(V\le-\kappa p\), \(\kappa>0\),

\[
H(U,V)\le \frac{U^2}{2\kappa p}.
\]

The numerical cutoff, quadrature grid and source-sign screens do not change these identities. A full upper bound on the nonnegative integral is necessary for a certificate.

## The two retained physical row banks

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


## Forward harmonic constraints

For each retained \(j\ge2\),

\[
\boxed{H_j=\langle F\operatorname{Im}(\omega^j f_\ell)\rangle=0.}
\]

Thus \(u_j=F\operatorname{Im}\omega^j\), \(v_j=F\operatorname{Re}\omega^j\). Saved names `H2` through `H40` mean these raw harmonics. The optional conditioning basis \(F[\omega^j-(-1)^j\omega^2]\) spans the same space when \(H_2\) is retained, but it is not the basis of the exported saved multipliers.

## Resummed constraints

Use \(r\) for the radius and \(\zeta\) for the cosine parameter, distinct from the symmetric point \(c\). For \(0\le r<1\), define

\[
m_c(w)=\frac{1-r\zeta w}{1-2r\zeta w+r^2w^2}
=\sum_{k\ge0}r^kT_k(\zeta)w^k,
\]

\[
m_s(w)=\frac{r w}{1-2r\zeta w+r^2w^2}
=\sum_{k\ge1}r^kU_{k-1}(\zeta)w^k,
\]

where \(T_k,U_k\) are Chebyshev polynomials. Each retained cosine or sine row is

\[
q_{r,\zeta}^{c/s}=F\omega^2[m_{c/s}(\omega)-m_{c/s}(-1)],
\qquad \boxed{\langle\operatorname{Im}(q_{r,\zeta}^{c/s} f_\ell)\rangle=0.}
\]

It is a convergent linear combination of valid harmonics, including its exact subtraction at the UV endpoint \(\omega=-1\). It does not impose each term in the infinite sum as an independent row. The stable exact differences used numerically are

\[
m_c(\omega)-m_c(-1)=
\frac{r(\omega+1)[\zeta+r-r(1+r\zeta)\omega]}
{(1-2r\zeta\omega+r^2\omega^2)(1+2r\zeta+r^2)},
\]

\[
m_s(\omega)-m_s(-1)=
\frac{r(\omega+1)(1+r^2\omega)}
{(1-2r\zeta\omega+r^2\omega^2)(1+2r\zeta+r^2)}.
\]

`Rc_r_zeta` and `rational_cos_r..._c...` are naming variants for the same cosine definition. Rational-row conditioning subtracts retained low harmonics during optimization; exported coefficients are transformed back to the definitions above.

## All seven crossing constraints, including their exact normalization

Every crossing row has \(u_i=0\). Let

\[
\mathcal X_{D,q;\ell}(\mu)=\frac1{Q^{D+1}}
\sum_{j=1}^{D-2}a_{D,q,j}Z^j P_\ell^{(j)}(Z).
\]

The corresponding identity is

\[
\boxed{X_{D,q}=\langle\mathcal X_{D,q;\ell}b_\ell\rangle=0.}
\]

These are **ordinary derivatives**, not derivatives divided by factorial. The full coefficient table is:

| Row | \(D\) | \(a_1\) | \(a_2\) | \(a_3\) | \(a_4\) | \(a_5\) | \(a_6\) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `X4` | 4 | \(-3\) | \(2\) | 0 | 0 | 0 | 0 |
| `X5` | 5 | \(8\) | \(-12\) | \(16/3\) | 0 | 0 | 0 |
| `X6` | 6 | \(-160\) | \(128\) | \(-64\) | \(64/3\) | 0 | 0 |
| `X7_1` | 7 | \(30\) | \(-70\) | \(72\) | \(-40\) | \(32/3\) | 0 |
| `X7_2` | 7 | \(48\) | \(-80\) | \(64\) | \(-32\) | \(128/15\) | 0 |
| `X8_1` | 8 | \(112\) | \(-144\) | \(160/3\) | \(32/3\) | \(-64/5\) | \(128/45\) |
| `X8_2` | 8 | \(-784\) | \(720\) | \(-1120/3\) | \(352/3\) | \(-128/5\) | \(256/45\) |

In particular,

\[
X_4=\left\langle\left[\frac{2P_\ell''(Z)}{Q^3\Delta^2}
-\frac{3P_\ell'(Z)}{Q^4\Delta}\right]b_\ell\right\rangle=0.
\]

For an independent exact construction, set \(\nu=s-2+t/2\), \(v=\nu^2\), \(w=t-4/3\). Crossing-invariant polynomials are generated by

\[
X=2v+\tfrac32w^2,\qquad Y=-wv+\tfrac14w^3.
\]

At degree \(D\), use coefficient slots \((n,k)=(n,D-2n)\), \(n=1,\ldots,\lfloor D/2\rfloor\). The exact null vectors in that order are:

| Row | Coefficient-slot weights |
|---|---|
| `X4` | \((1/2,-3/4)\) |
| `X5` | \((2,-1)\) |
| `X6` | \((16,8,-45)\) |
| `X7_1` | \((20,-3,0)\) |
| `X7_2` | \((16,0,-3)\) |
| `X8_1` | \((16,-16,9,0)\) |
| `X8_2` | \((32,40,0,-189)\) |

They annihilate every \(X^aY^b\) with \(2a+3b=D\). Their dispersive coefficient formula is

\[
\beta_j=2\sum_{(n,k):k\ge j}
q_{n,k}(-1)^{k-j}2^{-(k-j)}\binom{D-j}{k-j},
\qquad a_j=\frac{2^j}{j!}\beta_j.
\]

Equivalently the implementation uses

\[
\mathcal X_{D,q;\ell}=
\sum_j\beta_j\frac{(\ell+1)_j}{j!}
\frac{P_{\ell-j}^{(j,j)}(Z)}{\Delta^jQ^{D+1-j}}.
\]

Terms with \(j>\ell\) vanish. Here the last \(P\) is a Jacobi polynomial. The scalar \(j=0\) coefficient cancels exactly. The first row's primitive integer vector has been divided by four to match the manuscript's \(X_4\) normalization; omitting that division changes all quoted \(X_4\) multipliers.

## Dispersion identity underlying all Roy and strip rows

For external \((s,t)\), let \(u_e=4-s-t\) and \(y=1+2t/\Delta\). Define

\[
K_1(\mu;s,t)=\frac1{\mu-s}+\frac1{\mu-u_e}
-\frac1Q-\frac1{\mu-2c+t},
\]

\[
K_2(\mu;s,t)=\frac1{\mu-t}-\frac2Q+\frac1{\mu-2c+t},
\]

\[
\mathcal L_\ell(\mu;s,t)=K_1P_\ell(y)+K_2P_\ell(Z).
\]

The crossed-subtraction dispersion relation is

\[
\operatorname{Re}\mathcal A(s,t)-g_0
=\frac1\pi\operatorname{PV}\int_4^\infty
\big[K_1\operatorname{Im}\mathcal A(\mu,t)
+K_2\operatorname{Im}\mathcal A(\mu,c)\big]d\mu.
\]

The apparent pole at \(\mu-2c+t=0\) is removable for even spins. With

\[
D_\ell(x,y)=\frac{P_\ell(x)-P_\ell(y)}{x-y},
\]

continuously extended to \(x=y\),

\[
\frac{P_\ell(Z)-P_\ell(y)}{\mu-2c+t}=-\frac2\Delta D_\ell(y,-Z).
\]

The recurrence \(D_0=0,D_1=1\),
\((\ell+1)D_{\ell+1}=(2\ell+1)[xD_\ell+P_\ell(y)]-\ell D_{\ell-1}\), avoids subtractive division. The physical \(\mu=s\) pole is an actual principal value and must be retained analytically. Deleting its grid point is not a valid PV prescription.

## Folded Roy constraints

Set

\[
w_{Bn}(s)=\frac{\mathbf1_{4<s<B}}{B-4}
P_n\!\left(\frac{2(s-4)}{B-4}-1\right),
\qquad t(s,z)=-\frac{s-4}{2}(1-z).
\]

Smear the preceding dispersion residual against
\(w_{Bn}(s)ds\,2P_J(z)dz\), \(z\in[0,1]\), \(J\) even. Its raw identity is

\[
\langle u^{\rm raw}_{BJn}r+v^{\rm raw}_{BJn}b\rangle
+c_{BJn}g_0=0,
\]

\[
u^{\rm raw}_{BJn,\ell}(\mu)=
\frac{2\pi}{2J+1}\delta_{\ell J}w_{Bn}(\mu),
\]

\[
v^{\rm raw}_{BJn,\ell}(\mu)=
-\int_4^B ds\,w_{Bn}(s)\int_0^1dz\,2P_J(z)
\mathcal L_\ell(\mu;s,t(s,z)),
\]

\[
c_{BJn}=-2\delta_{J0}\delta_{n0}.
\]

The exported homogeneous row is

\[
\boxed{(u_{BJn},v_{BJn})=
(u^{\rm raw}_{BJn}+c_{BJn}u_0,
 v^{\rm raw}_{BJn}+c_{BJn}v_0).}
\]

Thus the \(J=n=0\) row explicitly reconstructs the subtraction constant. It is not a pure imaginary null before this matching. The factor two from folded angular integration and the factor \(1/(2J+1)\) in the local real kernel are both necessary. No extra \(2J+1\) multiplies the test function.

All retained endpoints are exactly \(B=119/10\), not \(11.999\). The external domain \(4<s<11.9\), \(0<z<1\) has \(t\ge-79/20>-4\). It is only the external test domain; spectral \(\mu\) still extends to infinity.

Using the real PV Legendre function

\[
Q_n^{\rm PV}(x)=\frac12\operatorname{PV}\int_{-1}^{1}\frac{P_n(y)}{x-y}dy,
\]

the projected PV term in the dispersive integral, before the minus sign in \(v^{\rm raw}\), is

\[
\delta_{\ell J}\frac{4Q_n^{\rm PV}(2(\mu-4)/(B-4)-1)}{(2J+1)(B-4)}.
\]

It vanishes exactly at spins above the largest retained \(J\). This fact is used in the tail estimates.

## Tapered fixed-transfer strip constraints

For each retained \([A,B]\), define

\[
x_s=\frac{2(s-A)}{B-A}-1,\qquad
w_{ABn}^{\rm tap}(s)=\frac{\mathbf1_{A<s<B}}{B-A}(1-x_s^2)P_n(x_s),
\]

\[
h_j(t)=\frac{\mathbf1_{-T<t<0}}{T}P_j(1+2t/T),\qquad T=4/3.
\]

These tests smear the same dispersion residual over external \(s\) and \(t\); they do not extrapolate a full-angle Roy projection to arbitrary external energy. All tested points are physical because \(T\le A-4\).

Let

\[
H_{\ell j}(\mu)=\int_{-T}^0 h_j(t)P_\ell(1+2t/(\mu-4))dt.
\]

Then

\[
u^{\rm raw}_{ABjn,\ell}(\mu)=\pi w_{ABn}^{\rm tap}(\mu)H_{\ell j}(\mu),
\]

\[
v^{\rm raw}_{ABjn,\ell}(\mu)=
-\int_A^B ds\,w_{ABn}^{\rm tap}(s)
\int_{-T}^0dt\,h_j(t)\mathcal L_\ell(\mu;s,t),
\]

\[
c_{ABjn}=-\delta_{j0}
\left(\frac23\delta_{n0}-\frac2{15}\delta_{n2}\right).
\]

The retained homogeneous row is again

\[
\boxed{(u_{ABjn},v_{ABjn})=(u^{\rm raw}_{ABjn}+c_{ABjn}u_0,
 v^{\rm raw}_{ABjn}+c_{ABjn}v_0).}
\]

Unlike the finite-\(J\) Roy local kernel, the strip local real kernel occurs at arbitrarily high even spin. The real-part infinite-spin remainder must include it.

The taper is exactly \(1-x^2=(2/3)P_0(x)-(2/3)P_2(x)\). Its PV transform is

\[
\widetilde Q_n(v)=(1-v^2)Q_n^{\rm PV}(v)+v\delta_{n0}+\tfrac13\delta_{n1}.
\]

Endpoint limits use the continuous limit of \(\epsilon\log|\epsilon|\); endpoints are not silently discarded. The legacy row label `T1.33333` denotes **exactly \(4/3\)**. Interpreting it as the decimal \(1.33333\) would reconstruct a different kernel.

After combining all subtraction reconstructions, the coefficient of the bare forward \(g_0\) kernel is

\[
\alpha=a-2\lambda_{\mathrm{Roy},0,0}
-\frac23\sum_{[A,B]}\lambda_{\mathrm{strip},0,0}
+\frac2{15}\sum_{[A,B]}\lambda_{\mathrm{strip},0,2}.
\]

The target is still \(a g_0+b g_2\). The effective \(\alpha\) is not a replacement target normal.

## Regular forward-transfer derivative constraints

These 36 rows occur only in Bank B; the derivative-order symbol \(p\) in this subsection is distinct from the positive kernel \(p_\ell\). At fixed spectral \(\mu>4\), let

\[
A_t=2+t/2,\qquad Q_t=\mu-2+t/2,\qquad
B_t=\sqrt{(\mu-4)(\mu+t)},
\]

\[
\omega_t=\left(\frac{A_t+iB_t}{Q_t}\right)^2,
\qquad F_t=\frac{8A_t^2}{Q_t^3},
\qquad P_t=P_\ell(1+2t/\Delta).
\]

For \(p=1,2,3\), \(j=2,\ldots,12\), the complex kernel is

\[
q_{p,j}(\mu,\ell)=
\left.\frac{d^p}{dt^p}
\left[F_tP_t\omega_t^j(1-\omega_t)^{2p-2}\right]\right|_{t=0},
\]

and the row is

\[
\boxed{\mathrm{Dt}p\mathrm{H}_j:
\quad\langle\operatorname{Im}(q_{p,j}f_\ell)\rangle=0.}
\]

The factor \((1-\omega_t)^{2p-2}\) is inserted **before differentiation**. It is a finite linear combination of higher harmonics that regularizes the threshold behavior; dropping it would define a different and potentially divergent dual functional.

For the remaining row at each \(p\), define the real subtraction kernel

\[
S_\ell(\mu,t)=P_t
\left[\frac2{Q_t}-\frac1Q-\frac1{\mu-8/3+t}\right]
+P_\ell(Z)
\left[\frac1{\mu-t}-\frac2Q+\frac1{\mu-8/3+t}\right].
\]

It equals the dispersive kernel evaluated at the crossing-symmetric external energy \(s=2-t/2\). Define

\[
q_{p,G}=\left.\frac{d^p}{dt^p}
\left[F_tP_t\omega_t(1-\omega_t)^{2p-2}
-S_\ell(\mu,t)\right]\right|_{t=0}.
\]

Then

\[
\boxed{\mathrm{Dt}p\mathrm{G}_{\mathrm{regular}}:
\quad\langle\operatorname{Im}(q_{p,G}f_\ell)\rangle=0.}
\]

The subtraction \(S_\ell\) is real, so it changes only the \(b_\ell\) coefficient. It is **not multiplied** by \((1-\omega_t)^{2p-2}\). The \(j=1\) term has coefficient one in the regularized harmonic combination; the other terms are higher homogeneous harmonics. The derivative removes the remaining constant \(g_0\) subtraction. The code evaluates ordinary derivatives through exact truncated Taylor algebra and multiplies the Taylor coefficient by \(p!\); these are not finite-difference constraints.

## Numerical optimization and the complete search restrictions

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


### What the Bernstein source inequalities mean

Let \(y=s-4\), \(Y=B-4=79/10\), \(y_0=8/3\), and write

\[
y=y_0+(Y-y_0)\xi,\qquad
z=q\left(1-\frac{y_0}{y}\right),\qquad 0\le\xi,q\le1.
\]

For the 220-row Roy degrees, \(J_{\max}=6\), \(n_{\max}=14\). Multiply the source by \((y/Y)^6\):

\[
\mathcal P(\xi,q)=\left(\frac yY\right)^6
W_R(4+y,q(1-y_0/y)).
\]

This is a polynomial with maximum degrees \((20,6)\). On each of 16 equal \(\xi\)-intervals and four equal \(q\)-intervals expand both \(\mathcal P\) and \(\mathcal R=(y/Y)^6\) in the tensor Bernstein basis of degree \((20,6)\). Every coefficient of \(\mathcal R\) is positive. The \(21\times7\times16\times4=9{,}408\) inequalities are

\[
\frac{\mathcal P_{I,J}}{\mathcal R_{I,J}}\ge m_R
\]

with the requested source margin and the additional optimizer cushion. They are linear in the multipliers. If checked exactly after freezing the coefficients, they prove \(W_R\ge m_R\) throughout \(t\le-4/3\) in the stated external Roy domain. Floating feasibility alone is not that certificate. Particular compact proofs may require and separately prove a larger source-sign domain; the fixed proof records specify it.

### What the crossing-series inequalities mean

With the exact crossing coefficients \(\beta_{i,j}\) from the crossing subsection define

\[
\mathcal F(k,z)=\sum_i\lambda_i\sum_j
\beta_{i,j}\frac{k(k-1)\cdots(k-j+1)}{j!c^j}
z^{D_i-2-j},\qquad c=4/3.
\]

For \(z=1/Q\in[0,3/8]\), the combined normalized crossing kernel is a positive weighted average of \(\mathcal F(k,z)\). Indeed, with \(x=c/\Delta\),

\[
P_\ell(1+2x)=\sum_{k=0}^{\ell}
\frac{(\ell+k)!}{(\ell-k)!(k!)^2}x^k,
\]

\[
\frac{\sum_i\lambda_i\mathcal X_{i,\ell}}p
=\frac{\sum_{k=0}^{\ell}\frac{(\ell+k)!}{(\ell-k)!(k!)^2}x^k\mathcal F(k,1/Q)}
{\sum_{k=0}^{\ell}\frac{(\ell+k)!}{(\ell-k)!(k!)^2}x^k}.
\]

Substitute \(k=10+y\), expand in \(y\), and put each coefficient polynomial in the degree-five Bernstein basis on \(z\in[0,3/8]\). Require all Bernstein coefficients of the \(y^0\) polynomial to be at most \(-1\) and those of every higher power to be nonpositive. The 42 nonzero constraints prove the sufficient search property \(\mathcal F(k,z)\le-1\) for every real \(k\ge10\). The optimizer's special base offset cancels the common \(10^{-8}\) cushion here. Low coefficients \(k<10\) still need to be included in the actual finite-spin and tail proofs. This component alone is not the full \(V/p\) sign certificate.


### Upper-sector additional guards

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

## Selected search argument crosswalk

The following entries are taken from saved run records. Intermediate construction parents are listed when their coefficients enter the selected vector. Complete file paths, exact hashes, warm starts, guard lists, and capture metadata are in `optimization_records.json`. A dash means the field is absent in that record.

| Run record | Tail margin | Source margin | Finite-bank margin | Eigen floor | Iteration budget |
|---|---|---|---|---|---|
| angle_000p000 | — | — | — | — | — |
| angle180_cut2 | 0.002 | 0.0001 | — | 1e-07 | 8 |
| angle180_refined | 0.002 | 0.0001 | — | 1e-07 | 80 |
| angle180_final | 0.002 | 0.0001 | — | 1e-07 | 40 |
| angle255 | 0.002 | 0.0001 | — | 1e-07 | 28 |
| angle259p2941_from255 | 0.005 | 0.0001 | — | 1e-07 | 30 |
| angle261p3814 | 0.005 | 0.0001 | — | 1e-07 | 40 |
| angle265 | 0.005 | 0.0001 | — | 1e-07 | 40 |
| angle268 | 0.005 | 0.0001 | — | 1e-07 | 40 |
| angle85_220 | 0.002 | 0.0001 | — | 1e-07 | 20 |
| angle263p2747_from261 | 0.005 | 0.0001 | — | 1e-07 | 40 |
| angle269p5 | 0.005 | 0.0001 | — | 1e-07 | 40 |
| angle83p1759644845 | 0.01 | 0.0001 | 0.002 | 1e-07 | 20 |
| angle83_refined | 0.01 | 0.0001 | 0.002 | 1e-09 | 20 |
| angle83_repaired | — | — | — | — | — |
| angle86_blended | 0.01 | 0.0001 | 0.002 | 1e-07 | 20 |
| angle86_refined | 0.01 | 0.0001 | 0.002 | 1e-09 | 40 |
| g2_m1e3_t02 | 0.02 | 0.0001 | 0.02 | — | 40 |
| candidate | 0.01 | 0.0001 | 0.002 | 1e-07 | 20 |
