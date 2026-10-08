# A route from the joint limit to finite-spin control

This note gives analytic inequalities useful for a future uniform error
bound. It does not contain a completed finite-spin certificate.

Let T=l(l+1), tau=T/(mu-4). For a physical source transfer t<=0, put
b=(-t)tau and y=1-2b/T. Suppose0<=b/T<1, so y is physical. Write
theta=2 asin sqrt(b/T). The [Mehler–Dirichlet integral in DLMF14.12.1](https://dlmf.nist.gov/14.12.E1), with order0, gives

    P_l(cos theta) = (2/pi) integral_0^1
      cos((l+1/2)theta u) / sqrt(1-u^2)
      /sqrt(sinc(theta(1+u)/2) sinc(theta(1-u)/2)) du.

The corresponding integral with unit sinc factor and phase2sqrt(b) is
J_0(2sqrt(b)). For0<=theta<pi, sinc is positive and decreasing, hence

    |P_l(1-2b/T)-J_0(2sqrt(b))|
       <= |(l+1/2)theta-2sqrt(b)| + 1/sinc(theta)-1.

This avoids the exponentially large absolute-series error for oscillatory
Legendre factors. If beta=b/T and beta<3/5, elementary sine and arcsine
inequalities give the explicit convenient upper bound

    2sqrt(b) [sqrt(1+1/(4T))/sqrt(1-beta)-1]
       + 4beta/(6-10beta).

It is O((b^(3/2)+b)/T) on any fixed bounded b interval. More precisely,
the first term uses asin(x)/x <=1/sqrt(1-x^2), and the second follows
from theta^2<=4beta/(1-beta) and sin(theta)>=theta-theta^3/6.

The denominator at Z=1+2a/T, a=(4/3)tau, has the positive expansion

    P_l(Z)=sum_k a^k/(k!)^2 product_{j=0}^{k-1}(1-j(j+1)/T),

where coefficients vanish beyond k=l. The product inequality1-prod(1-x_j)
<=sum x_j gives

    0 <= I_0(2sqrt(a))-P_l(Z)
       <= [(D^3-D) I_0(2sqrt(a))]/(3T), D=a d/da.

For k>l the same bound follows because the sum of j(j+1)/T through
j=k-1 is already at least1. Since D^2 I_0=a I_0,

    (D^3-D) I_0 = a I_0 + (a-1)sqrt(a) I_1,

and0<=I_1/I_0<=1, the relative denominator error is at most

    [a+max(a-1,0)sqrt(a)]/(3T).

For example at tau<=256 and l>=1000 the latter is below about0.00221.
Around the shallow minimum tau~2.3 it is much smaller. These estimates
suggest that a uniform finite-spin proof with a cutoff in the thousands
may be practical, but that cutoff is not yet established.

To finish, one still needs:

- exact integral bounds for source|WK| with the t-dependent Mehler error;
- finite-mu corrections to Q^3 K1 and Q^3 K2, keeping source signs;
- corresponding finite-spin corrections to crossing and differentiated
  harmonics, where cancellation near the shallow limiting maximum matters;
- separate large-tau/compact-energy estimates and finite-spin integration.

The exact crossing-series certificate can help directly: its coefficient
polynomial is nonpositive beyond k=10 for the old witness, so an arbitrarily
long negative series tail can be dropped with the correct inequality
direction. These estimates must be regenerated for changed coefficients.
