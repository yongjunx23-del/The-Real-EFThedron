For the massive symmetric-point weight, write Q=s-4/3, delta=s-4, Z=Q/delta=1+8/(3 delta), and p_l=P_l(Z)/Q^3. Let L be an even head cutoff. Positivity of every Taylor coefficient of P_l(1+y) gives, for any integer 2<=k<=L+2,

P_l(Z) >= [(l+k)!/((l-k)!(k!)^2)] * [4/(3 delta)]^k, for every l>=L+2.

The remaining rational spin sum telescopes in steps of two:

sum_{l=L+2,L+4,...} (2l+1) (l-k)!/(l+k)!
 = 1/[2(k-1) product_{r=L-k+3}^{L+k} r].

Consequently the omitted inverse-Legendre sum is at most

B_k(s) = (k!)^2 [3 delta/4]^k / [2(k-1) product_{r=L-k+3}^{L+k} r].

The minimum of these bounds at each energy remains an upper bound. Higher k is especially valuable on the finite strip windows; k=2 or 3 controls the far ultraviolet. This does not impose a spin cutoff on the physical amplitude.

For spins exceeding all retained angular Roy test spins, U_l is the spin-independent forward real kernel plus the physical-angle projection of the fixed-transfer strip weights. The elementary |P_l(x)|<=1 for -1<=x<=1 gives a common envelope Ubar(s). Then

A_tail = sum_{l>L} integral U_l^2/p_l dnu_l
 <= integral_4^infinity [Ubar(s)^2 Q^3/(2 pi rho)] min_k B_k(s) ds.

A sharper strip envelope uses the classical Bernstein bound

|P_l(cos(theta))| <= sqrt(2/[pi(l+1/2)sin(theta)]).

The exact l+1/2 form is stated as Lemma 2.1, Eq. (2.2), in Haiyong Wang, [A new and sharper bound for Legendre expansion of differentiable functions](https://arxiv.org/pdf/1803.00336), with the endpoint statement understood by continuity of the weighted inequality.

On a strip t=-T v, 0<v<1, sin(theta)^2=4 T v/(s-4) * (1-T v/(s-4)). The resulting v^(-1/4) endpoint weight is integrable. Replacing l by its smallest tail value L+2 gives a common envelope for every omitted spin. The script evaluates this angular integral by Gauss-Jacobi quadrature and the scalar energy integral adaptively; these evaluations are not directed interval enclosures.

If an independently valid global g2 bound is G and the omitted absorptive kernel obeys V_l/p_l<=eta, Cauchy-Schwarz and r_l^2<=2 b_l imply

tail objective <= sqrt(2 G A_tail) + G max(eta,0).

Alternatively, a proved V_l<=-kappa p_l would give the tighter disk support tail <=A_tail/(2 kappa). Neither V inequality is established by the A-tail calculation. Sampled negative V values are diagnostics, not a proof of the infinite-spin inequality.
