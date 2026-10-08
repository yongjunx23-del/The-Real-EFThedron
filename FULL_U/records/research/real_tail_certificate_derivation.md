# Directed certificate for the real-kernel tail

This is one component of a dual certificate. It does not certify the sign of the absorptive kernel \(V\), or the finite-spin support integral.

The executable is `research/certify_real_tail.py`. The fixed witness is `output/angles_resummed/angle_000p000.npz`, SHA-256 `bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c`. Coefficients are the exact binary64 numbers stored in that file; decimal row parameters are interpreted as exact rationals, and the canonical transfer width is exactly \(T=4/3\). This convention agrees with the absorptive-tail certificate. Ordinary floating search/replay approximates these intended parameters.

Write \(d=\mu-4\), \(Q=\mu-4/3\), \(Z=Q/d\), \(p_l=P_l(Z)/Q^3\), and \(d\nu_l=(2l+1)d\mu/(2\pi\sqrt{1-4/\mu})\). We certify an upper bound on

\[
A_L=\sum_{l=L+2,L+4,\ldots}\int_4^\infty U_l^2/p_l\,d\nu_l.
\]

For even \(L\) above all local Roy spins, the real kernel contains only a spin-independent forward part and physical-angle strip projections. Crossing rows have no real kernel. The subtraction-elimination coefficient is

\[
\alpha=a-2\lambda_{\mathrm{Roy},J=0,n=0}
-\frac23\sum_{\mathrm{windows}}\lambda_{j=0,n=0}
+\frac2{15}\sum_{\mathrm{windows}}\lambda_{j=0,n=2}.
\]

The final term is essential for power-one tapered strips.

## Exact infinite-spin summation inequality

For every integer \(2\le k\le L+2\), positivity of the Taylor coefficients of \(P_l(1+x)\) gives

\[
P_l(Z)\ge\frac{(l+k)!}{(l-k)!(k!)^2}\left(\frac4{3d}\right)^k.
\]

The remaining even-spin sum telescopes exactly:

\[
\sum_{l=L+2,L+4,\ldots}(2l+1)\frac{(l-k)!}{(l+k)!}
=\frac1{2(k-1)\prod_{r=L-k+3}^{L+k}r}.
\]

Therefore

\[
\sum_{l>L,\mathrm{even}}\frac{2l+1}{P_l(Z)}\le B_k(\mu)=C_kd^k,
\quad C_k=\frac{(k!)^2(3/4)^k}{2(k-1)\prod_{r=L-k+3}^{L+k}r}.
\]

All \(C_k\) are exact rational numbers. Any one \(k\) is valid; choosing a different \(k\) in each integration cell preserves the upper bound. No physical spin cutoff is assumed.

## Forward kernel with its UV zero factored out

Set \(\beta=\sqrt{1-4/\mu}\in[0,1]\),

\[
x=\Re\omega=\frac{1-6\beta^2+\beta^4}{(1+\beta^2)^2},\qquad
\Im\omega=\frac{4\beta(1-\beta^2)}{(1+\beta^2)^2}.
\]

If the complex forward combination is \(Fq(\omega)\), \(F=32/(\mu-2)^3\), define the real rational function \(K(x)=\Im q(\omega)/\Im\omega\), extended continuously at the endpoints. Then

\[
U_{\mathrm f}(\beta)=\frac{16\beta(1-\beta^2)^4}{(1+\beta^2)^5}K(x).
\]

For harmonics, \(K=\alpha+\sum_j\lambda_j U_{j-1}(x)\), with the rational generators added explicitly. For a rational generator \(m(\omega)=n(\omega)/d(\omega)\), write \(n=n_r+i y n_i\), \(d=d_r+i y d_i\), \(y=\Im\omega\). Its contribution after multiplication by \(\omega^2\) and subtraction of \(m(-1)\) is

\[
(2x^2-1)\frac{n_i d_r-n_r d_i}{|d|^2}
+2x\left(\frac{n_r d_r+(1-x^2)n_i d_i}{|d|^2}-m(-1)\right).
\]

For \(d=1-2rc\omega+r^2\omega^2\), the positive denominator is evaluated as

\[
|d|^2=[2rx-c(1+r^2)]^2+(1-c^2)(1-r^2)^2.
\]

After including \(Q^3\), the spectral measure and \(B_k\), the forward density is

\[
D_k(\beta)=\frac{65536}{27\pi}C_k4^k
\frac{\beta^{2k+2}(2+\beta^2)^3(1-\beta^2)^{3-k}}
{(1+\beta^2)^{10}}K(x)^2.
\]

The last cell touching \(\beta=1\) uses only \(k=2,3\), so the infinite-energy endpoint is enclosed without a divergent intermediate expression. There is no discarded UV interval.

For stable interval evaluation of the harmonic polynomial, a degree-five Taylor expansion about the interval midpoint is used. Its sixth-order remainder follows from

\[
\max_{|x|\le1}\frac{|U_n^{(p)}(x)|}{p!}
=\frac{2^p(n+p+1)!}{(2p+1)!(n-p)!},\quad n\ge p,
\]

summed with absolute multiplier weights. The identity follows from the derivative relation for Gegenbauer polynomials; their maximum magnitude on \([-1,1]\) is attained at 1 for these positive parameters. Coefficients below the remainder order are evaluated by differentiated recurrence in Arb. The exact monotonicity of \(x(\beta)=\cos(4\arctan\beta)\) supplies endpoint enclosures, including the extrema \(x=1,-1\).

## Certified physical strip envelope

On one window \([A,B]\), define \(\xi=2(\mu-A)/(B-A)-1\) and

\[
W_\mu(z)=\sum_j A_j(\xi)P_j(z),\quad
A_j(\xi)=\sum_n\lambda_{jn}P_n(\xi),\qquad
\|W_\mu\|_2^2=\sum_j\frac{A_j(\xi)^2}{2j+1}.
\]

Here the norm uses the normalized measure \(dz/2\). The taper is exactly \(1-\xi^2\). The elementary bound \(|P_l|\le1\) gives an envelope with factor 1. A sharper one follows from [Wang, Lemma 2.1, Eq. (2.2)](https://arxiv.org/pdf/1803.00336):

\[
|P_l(\cos\theta)|\le\sqrt{\frac{2}{\pi(l+1/2)\sin\theta}}.
\]

Write \(t=-Tv\), \(0\le v\le1\), so \(\sin^2\theta=4Tv/d\,(1-Tv/d)\). Cauchyâ€“Schwarz and \(\int_0^1v^{-1/2}dv=2\) yield

\[
|U_{\mathrm{strip},l}(\mu)|\le\frac\pi{B-A}\(1-\xi^2\)\|W_\mu\|_2
\sqrt{\min\left(1,\frac{2d}{\pi(L+5/2)\sqrt{T(d-T)}}\right)}.
\]

This holds for every \(l\ge L+2\). The source windows are physical: \(d\ge A-4>T\). The factor is evaluated with directed interval bounds; the spectral integral contains the square of \(|U_{\mathrm f}|+\sum|U_{\mathrm{strip}}|\). Cells crossing a window endpoint include that strip conservatively. The windows in the supplied witness are disjoint, \([12,24]\) and \([24,48]\).

The optional older route uses a constant strip bound plus Minkowski; it is much weaker. The `--combined-strips` route above is used for the final component records.

## Directed integration and saved evidence

The integration is an upper Riemann sum, not a Gauss rule with an estimated error. Every dyadic cell in \([0,1]\) is evaluated as an Arb ball at 160-bit precision. Its upper density is multiplied by the exact cell width. Cells with the largest contributions are bisected until a requested upper target or cell-count limit is reached. A failed accuracy target still leaves a valid, looser upper sum, provided every cell is finite.

The output records all cell endpoints, exact dyadic upper masses, chosen \(k\), witness hash, precision, and the total dyadic upper endpoint. The partition and its mass sum were independently checked using exact rational arithmetic. The component proof was independently reviewed against the displayed inequalities. This is executable ball-arithmetic certification, not a proof-assistant formalization.

For the original upper witness, the retained directed bounds are:

- \(A_{40}<2\times10^{-5}\), from `output/certification_20261006/right_A40_arb.json`.
- \(A_{100}<3\times10^{-8}\), from `output/certification_20261006/right_A100_arb.json`.

These alone do not imply a bound on \(g_0\). If a separate proof gives \(V\le-\kappa p\) on the same omitted domain, its full pointwise disk support is at most \(A_L/(2\kappa)\). If only \(V/p\le\eta\) and an independently valid \(g_2\le G\) are available, the physical tail contribution is at most \(\sqrt{2GA_L}+G\max(\eta,0)\). Finite-head integration must also be controlled.

## Reproduction

Install `python-flint==0.9.0` into an isolated environment or `tmp/certification_deps`; the research environment already contains NumPy. Run from the workspace root:

```text
python research/certify_real_tail.py output/angles_resummed/angle_000p000.npz --L 40 --combined-strips --target 2e-5 --max-cells 8000 --out output/certification_20261006/right_A40_arb.json
python research/certify_real_tail.py output/angles_resummed/angle_000p000.npz --L 100 --combined-strips --target 3e-8 --max-cells 8000 --out output/certification_20261006/right_A100_arb.json
```

Do not apply this certificate to a different coefficient file, to transfer-derivative rows, or to overlapping source windows without extending the proof and implementation. A fresh candidate requires its own hash-linked run.
