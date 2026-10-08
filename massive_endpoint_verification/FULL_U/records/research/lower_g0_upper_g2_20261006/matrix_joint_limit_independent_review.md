# Independent review of the joint-limit certificates

Read-only mathematical/code review of `tail_certify_joint_limit.py` and `tail_candidate_joint_certificate.py`, 6 October 2026. This did not independently rerun every interval cell or reimplement the spectral kernels.

No substantive error was found in the reviewed formulas or interval coverage.

For the source function `S(tau)=integral WK J_0(2 sqrt((-t) tau))`, the exact moments are `M_k=integral WK (-t)^k`. Expanding `J_0` gives `S=sum (-1)^k M_k tau^k/(k!)^2`. The code's centered Taylor coefficient of order `q` is therefore correctly

`sum_k (-1)^(k+q) M_(k+q) mid^k/[k! (k+q)! q!]`.

The absolute moment majorant follows from termwise absolute values of the source polynomial on the positive unit square. The omitted series begins at `k=N`; the ratio of its positive majorizing terms is at most `a/[(N+1)(N+q+1)]`, as implemented.

The sixth-order remainder over a `tau` cell is bounded by `absMom[6]*h^6/(6!)^2`. One factorial comes from the derivative identity for `0F1`, and one from Taylor's theorem. The estimate `|0F1(;q+1;-a)|<=1` for real `a>=0` follows from the standard Poisson integral for the corresponding Bessel function; all cells have nonnegative `tau`.

The derivative-null expressions have the correct product-rule factors from differentiating `(32+16t+2t^2) 0F1(;1;t tau)`. The `G_regular` subtraction additionally contains the correct polynomial derivative terms and the constants `-16/3`, `4`, and `0` for derivative orders one, two, and three.

The new rational generators have `r=0.8<1` and the UV subtraction `m(omega)-m(-1)`. They are regular at `omega=-1`, so this difference vanishes as `omega+1`. Their normalized contribution to the joint limiting absorptive kernel is exactly zero. This conclusion relies on the implemented subtraction and on poles staying away from `omega=-1`; it is not true for arbitrary rational generators without that subtraction.

For large `tau`, the exact crossing coefficient certificate makes all Legendre-series terms above its cutoff nonpositive. Dropping those terms, dropping negative coefficients from the remaining numerator polynomial, and using `|S|<=C` yields a positive polynomial majorant. One positive Taylor term of `I_0` provides the denominator lower bound. The chosen denominator degree exceeds every numerator degree, making each resulting power decrease for `tau>=256`. The exact rational endpoint evaluation therefore gives a valid uniform large-`tau` bound.

The compact certificate sorts its exact rational partition and checks adjacency from `0` through `256`. The large-`tau` argument covers `[256,infinity)`. The replay checks the witness hash and rebuilds the exact moment data. Together these establish the claimed sign for the **joint limiting function**.

They do not establish a uniform error between the finite-spin kernel and that limit, nor evaluate the physical real-kernel remainder or the full dual objective. The scripts' explicit scope statement correctly retains these missing obligations.
