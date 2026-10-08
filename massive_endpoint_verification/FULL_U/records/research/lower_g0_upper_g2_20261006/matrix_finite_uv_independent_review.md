# Review of the finite-spin ultraviolet sign certificate

Reviewed `tail_finite_uv_region.py` and the stated domain of `tail_blend90_finite_uv_L10000_0_256_256.json` on 6 October 2026. This is an independent derivation/code review and partition check, not a second full interval-kernel replay.

No substantive mathematical error was found for the recorded witness and domain: every even `ell>=10000`, with `0<=tau=ell(ell+1)/(mu-4)<=256`.

The implementation uses `eta=1/[ell(ell+1)]`, `w=tau*eta=1/(mu-4)`. Its third-order complex jets implement the finite-energy forward and derivative kernels exactly. The transfer Legendre coefficients `1,tau,tau^2(1-2eta)/4,tau^3(1-2eta)(1-6eta)/36` are the exact first four Taylor coefficients. I separately checked the two rational, energy-factored subtraction kernels `Ka` and `Kb` with symbolic algebra; both agree identically with `Q^3` times their original definitions.

The UV-subtracted rational harmonic expressions also factor correctly. For the cosine generator the numerator of `m(omega)-m(-1)` is

`r(omega+1)[c+r-r(1+rc)omega]`;

for the sine generator it is `r(omega+1)(1+r^2 omega)`. The code retains their finite-energy contributions using interval complex arithmetic, rather than merely setting them to their zero limiting values.

For the source kernel, put `a=s-c`, `b=c-u`, with `c=4/3`. Its exact finite-energy multiplier relative to `K1_infinity=2ab` is

\[
R_1=\frac{Q^2[Q+(t-c)/2]}
{(\mu-s)(\mu-u)(\mu-4+t+c)}.
\]

On the actual source domains, `s<=48`, `-4<=t<=0`, `u<=0`. The numerator-minus-denominator expansion in `1/mu` has first coefficient `2-3t/2<=8`; the second coefficient has absolute value below `2153`; the cubic coefficient has absolute value at most `2112*(20/3)+64/9<14088`. The denominator is bounded below by `(1-48/mu)(1-20/(3mu))`. This validates the implemented `eps1` bound.

All physical-source Legendre arguments lie in `[-1,1]` on this UV domain. Consequently

`|R1 P_l-J0| <= |R1-1| + |P_l-J0|`.

The first term uses `|P_l|<=1`; there is no omitted product of the two errors. The Mehler–Dirichlet comparison correctly bounds the second term by

`2 sqrt(b)[sqrt(1+1/(4T))/sqrt(1-beta)-1] + 4 beta/(6-10 beta)`,

where `b=(-t)tau`, `beta=b/T`, and `T=ell(ell+1)`. The domain has `beta<3/5` by a wide margin, making every denominator positive. The bound follows from the unit-mass Mehler measure, the Lipschitz bound for cosine, monotonicity of sinc, `asin(x)/x<=1/sqrt(1-x^2)`, and `sin(theta)>=theta-theta^3/6`.

For the second source kernel,

`H_finite/H_infinity = 1/[1-((t-c)/Q)^2]`.

The implemented additive error `norm2*r^2/(1-r^2)` is therefore valid. The source norms are taken from the exact Bernstein integral bounds for the same witness hash.

The denominator `P_l(1+2c tau/T)` has nonnegative Taylor coefficients. Its first 64 terms evaluated with the smallest allowed `T` and the left cell endpoint provide a lower bound; `0F1(;1;c tau_hi)` provides an upper bound. All omitted crossing numerator terms have index at least 64 and are nonpositive by the exact coefficient certificate, whose cutoff is 10. No finite crossing term is omitted with an uncontrolled sign. Since `ell>=10000`, all retained 64 terms belong to the actual finite Legendre polynomial.

The saved intervals form a gapless partition of `[0,256]` and the record reports no failed cells. Thus this is a finite-spin, finite-energy sign component, stronger in scope than a joint-limit-only calculation.

It still covers only the stated portion of the spin-energy plane. In particular, it does not cover the complementary region `tau>256`, intermediate omitted spins, or the real-kernel support cost, and does not certify the final numerical coupling bound.
