# Lower-g0 joint UV limit: new certificate component

This is a certificate for a limiting kernel, **not yet a certified lower bound
on g0**. It concerns only the fixed old witness
`output/angle_constraints/ming0_scaled_threshold_repaired.npz`, SHA256
`c7efae0020c9883d6949978ff924d11c3b828f470f10e24ff813f993e47cc06a`.
It does not automatically apply to any newly optimized witness.

## Why lower g0 is the easier target

The old numerically audited lower-g0 candidate has measured head cost
7.781765419411301 through spin40 and an additional restricted-band tail
0.0010173390549025035 through spin100. Their sum7.782782758466204 is below
the published graphical dual support8.014557983981451 by0.231775225515247.
However, omitted energies/spins and numerical integration errors remain
unbounded, so the sum is a candidate value, not a certified inequality.
Existing upper-g2 searches are about1.1323 versus published graphical1.08168
and have additional large-spin feasibility problems.

## Joint limit

Set c=4/3, T=l(l+1), tau=T/(mu-4), and hold tau fixed while l and mu grow.
Then P_l(1+2c/(mu-4)) tends to

    P(tau) = I_0(2 sqrt(c tau)) = 0F1(;1;c tau).

For physical source t<=0, P_l(1+2t/(mu-4)) tends to
J_0(2 sqrt(-t tau)). For the finite source density W used by the witness,
set K=2(s-c)(c-u), H=2(t-c)^2, u=4-s-t. With alpha equal to the total
forward-g0 coefficient, its combined source/forward contribution is

    -32 alpha/9 - integral W H
    + [ -280 alpha/9 + 32 sum_j lambda_Hj (-1)^j
        - integral W K J_0(2 sqrt(-t tau)) ] / P(tau).

All source Jacobians are included in W. For Roy, x,z in[0,1],
W=2 sum lambda_Jn P_n(2x-1)P_J(z), s=4+(119/10-4)x,
t=-(119/10-4)x(1-z)/2. For a tapered strip,
W=sum lambda_jn[1-(2x-1)^2]P_n(2x-1)P_j(2z-1),
s=A+(B-A)x, t=-(4/3)(1-z).

Each degree-D crossing row contributes only its j=D-2 leading derivative:

    lambda_D beta_j tau^j 0F1(;j+1;c tau) / [(j!)^2 P(tau)].

For derivative harmonics of order p, define

    C_p(tau)=32 tau^p/p! +16p tau^(p-1)/(p-1)!
             +2p(p-1) tau^(p-2)/(p-2)! ,

with absent negative powers omitted. A Dt_p H_j row contributes
(-1)^j 4^(p-1) C_p/P. A Dt_p G_regular row has j=1 and subtracts
the p-th derivative at t=0 of

    [(-8/9+4t/3-t^2/2) I_0(2 sqrt(t tau))
       + P(tau)(32/9-16t/3+2t^2)] / P(tau).

The implementation retains all these terms; it does not discard the
derivative harmonics when taking the joint limit.

## Rigorous source integration and compact tau arithmetic

`tail_certify_joint_limit.py` uses binary64 multipliers as exact binary
rationals and exact rational source endpoints. Expand W K in monomials
x^a z^b. The moments M_k=integral W K (-t)^k are exact because

    integral_0^1 z^b(1-z)^k dz=b!k!/(b+k+1)!,

and the x integral is1/(a+k+1) for Roy or1/(a+1) for strips.
The source function is sum_k (-tau)^k M_k/(k!)^2. At each cell midpoint,
its derivatives through order5 are evaluated with Arb and200 series terms,
including an explicit geometric upper bound for the omitted series tail.
The cell enclosure uses a fifth-order centered Taylor polynomial plus
the sixth-order remainder

    integral |WK| |t|^6 * h^6 / (6!)^2.

This follows from the Poisson-integral bound
|J_q(x)| <= (|x|/2)^q/q! for real x, integer q>=0. Absolute source norms
are enclosed by exact integrals of absolute monomial coefficients.
This norm estimate can be loose but has the correct inequality direction.

The first256-bit run certified every point in[0,256] using53,206 gapless
accepted dyadic cells (106,411 evaluations). Numerical diagnostics give
a maximum about-0.010809035 near tau2.29594. That diagnostic maximum
is not itself a certified optimum.

A separate320-bit replay rebuilt all207 exact moments, matched them
exactly, and freshly reevaluated all53,206 cells. A subsequent320-bit
refinement certified the explicit uniform margin

    limiting V/p <= -1/200

on[0,256], with56,321 accepted cells and59,436 evaluations. Together
with the large-tau proof below this gives the same uniform margin for
**all tau>=0**, within the joint limiting function.

## Exact closure for tau>=256

The independent crossing-series certificate writes crossing/P as

    sum_k (c tau)^k F(k,0)/(k!)^2 / P(tau),

and proves F(k,0)<=0 for k>=10 using exact rational polynomial signs.
Drop that nonpositive tail. Bound the remaining source numerator using
|J_0|<=1 and retain only positive coefficients of the degree<=9 numerator.
Since

    P(tau) >= (c tau)^18/(18!)^2,

every retained ratio is decreasing for tau>=256. Exact rational evaluation
at256 gives the uniform upper bound

    limiting V/p <= -10.246672513611358...   (tau>=256).

The rational value, rather than the printed decimal, is stored in
`tail_joint_large_tau_certificate.json`.

## Scope and missing work

These statements certify the analytic joint limiting expression and its
floating-point implementation with exact moment arithmetic. The320-bit
replay freshly reconstructs the exact source moments and reevaluates every
compact cell, using the same kernel implementation; it is not an
independently written verification program.

They do **not** establish a uniform finite-l error to this limit. A limiting
sign alone does not rule out islands whose position changes with l. The
old witness already has integrable positive-V islands around spins60–70,
which are visible in its numerical tail audit. To certify an actual g0
inequality one must still control finite spins over all energies, all
remaining spins over compact energies, and a uniform joint-limit error,
then bound every support integral and quadrature error. These certificate
components cannot justify relabeling the candidate7.78278 as rigorous.
