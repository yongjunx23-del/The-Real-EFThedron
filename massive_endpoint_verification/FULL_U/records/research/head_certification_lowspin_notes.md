# Low-spin component certificate

This note concerns only the fixed old 150-row right-point witness
`output/angles_resummed/angle_000p000.npz`, SHA256
`bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c`.
The implementation is `research/head_certification_lowspin.py`. It treats
spins 2, 4, and 6 over the entire spectral interval, including both endpoints.

Write c=4/3, d=mu-4, Q=mu-c, beta=sqrt(1-4/mu), e=1-beta^2,
and R_l(d,t)=sum_j binom(l,j)binom(l+j,j)d^(l-j)t^j. Thus
Rstar=R_l(d,c)>0, Pstar=Rstar/d^l, and p=Pstar/Q^3. The support
integral in beta is `4(2l+1)/pi * integral H(U/e^2,V/e^2) d beta`,
where H(U,V)=sqrt(U^2+V^2)+V. The stable negative-V evaluation and
convex interpolation implementation are shared with the spin-zero certificate.

The complete combined source kernel has the exact polynomial identity

    d^l L = R_l(d,t)/(mu-s) + R_l(d,t)/(mu-u) - R_l(d,t)/Q
            + B_l(d,t) + Rstar/(mu-t) - 2 Rstar/Q,
    B_l(d,t)=(Rstar-R_l(d,t))/(d+c+t).

For even l, B_l is a polynomial because R_l(d,-d-c)=Rstar. SymPy exact
rational polynomial division asserts that its remainder is zero. Its source
integral, and that of R_l, use the certified T_j polynomial moments.

For a Roy source, put r=s-4 and b=(1-z)/2, so t=-rb. The physical PV
denominator is removed by

    R_l(d,-rb)/(d-r) = d^l P_l(z)/(d-r)
      + sum_{j=1}^l (-1)^(j+1) C_lj b^j
          sum_{a=0}^{j-1} d^(l-1-a) r^a.

The second term's source moments are integrated as exact rationals. For the
first term, half-range orthogonality of even Legendre polynomials isolates
only source angular order J=l. If W_l(y)=sum_n lambda_{l,n}P_n(y), its
contribution to V/e^2 is -4 Q_W(y)/[(2l+1)(B-4)e^2], and its real local
contribution is 2pi W_l(y)/[(2l+1)(B-4)e^2]. For a tapered strip define
H_j(d)=integral_0^1 P_j(2z-1) R_l(d,c(z-1)) dz and
W_j(y)=(1-y^2)sum_n lambda_{jn}P_n(y). Its corresponding factors are
-2 H_j(d) Q_W(y)/[(B-A)d^l e^2] and
pi H_j(d) W_j(y)/[(B-A)d^l e^2].

All Q_W functions use the exact identity

    Q_W(y) = W(y)/2 log|(y+1)/(y-1)| - R_W(y),
    R_W(y)=sum_m a_m sum_{h=0}^{m-1} y^(m-1-h) mean(z^h).

Here mean(z^h)=0 for odd h, and 1/(h+1) for even h. When |y|>1.5,
the code instead uses N=64 terms of its inverse-y moment series and adds
the rigorous remainder norm*|y|^(-N)/(|y|-1). The norm is the sum of
absolute Legendre coefficients, valid also with the taper factor <=1.

The crossed denominators use certified polynomial moments of (u+h_u)^k
and (t+h_t)^k. Their centers and radii enclose the entire source windows.
For each denominator the omitted geometric series is bounded by
norm*q^N/(denominator-radius), times a positive-coefficient bound on
|R_l(d,t)| where appropriate. These are added as uniform errors, not
differentiated as constants. Root independently recomputed all moment tables
at Arb340 and checked fresh intervals were contained in every saved input
interval (`output/certification_head/all_moment_tables_verified.json`).

For mu>=128 the code uses factored UV arrays A_j(x), B(x), x=1/mu,
provided by `certified_uv_moments.py`. The first two denominator expansion
coefficients cancel exactly using s+t+u=4. The scaled source contribution is

    -x/16 [sum_j C_lj (x/(1-4x))^j A_j(x) + Pstar B(x)].

The retained 32 terms carry the explicit all-source geometric remainder
4 norm B_source^2 (B_source*x)^32/(1-B_source*x), multiplied by the
sum of absolute Legendre factors. This expression is finite at beta=1.

The threshold cap uses beta<=2^-14 for l=2 and beta<=2^-9 for l=4,6.
Its interval calculation proves V/p<0 throughout the cap. The only singular
PV term is bounded using d^l|log d|<=dmax^l|log dmax|, valid on these caps.
Then H<=U^2/(-2V), with an interval upper bound for 1/p, gives the cap mass.
The cap includes beta=0 without numerically evaluating d^-l. Ordinary bulk
cells below beta=1/2 can use the same normalized V/p calculation when its
upper endpoint is negative; this is much more stable than raw V intervals.

Small gaps around source endpoints mu=11.9,12,24,48 use exact rational
beta endpoints. The code separates polynomial log coefficients and bounds
integrated absolute logs analytically. Elsewhere the true local-source branch
in the second-derivative series is determined from exact rational beta
endpoints. Inflated interval images of mu are never used to choose that
analytic branch. Plain interval fields retain the union with zero at a
possibly intersected source endpoint.

Bulk support integration takes the minimum of two independently valid upper
bounds: an interval Riemann bound and convex interpolation using interval
second derivatives of U,V. The uniform geometric truncation error contributes
4*h*error to the unweighted interpolation integral. Exact rational endpoints,
directed Arb arithmetic, and saved dyadic cell upper endpoints make the
final sum checkable. Midpoint estimates only prioritize subdivisions; they
are not used as bounds or acceptance criteria.

Point validation against independently assembled 96x96 numerical kernels was
run for all three spins at mu=4.0001,4.02,4.5,5,7,10,11.8,11.99,12.1,18,
23.9,24.1,36,47.9,48.1,80,127,128,160,1000,1e6. Away from the threshold
large powers, scaled field errors were normally below 1e-12; threshold
relative errors were about 1e-9 or smaller, compatible with the floating
reference's conditioning. This check is validation only; the certificate
uses exact identities and explicit remainder bounds.

Commands (PowerShell, workspace root; the bundled Python path is local):

    $py='C:\Users\laure\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
    & $py research/head_certification_lowspin.py --spin 2 --validate
    & $py research/head_certification_lowspin.py --spin 4 --validate
    & $py research/head_certification_lowspin.py --spin 6 --validate
    & $py research/head_certification_lowspin.py --spin 2 --target .0116 --out output/certification_head/spin2_tight_component.json
    & $py research/head_certification_lowspin.py --spin 4 --target .0021 --out output/certification_head/spin4_tight_component.json
    & $py research/head_certification_lowspin.py --spin 6 --target .0009 --out output/certification_head/spin6_tight_component.json

The earlier files `spin2_component.json`, `spin4_component.json`, and
`spin6_component.json` predate the independently reviewed analytic-branch
correction and are superseded by the `*_tight_component.json` files. They
must not be used in the final certificate.

Final corrected results (Arb224; full spectral interval):

| Spin | Certified upper | Cells |
|---|---:|---:|
| 2 | 0.011599994410977523 | 4654 |
| 4 | 0.002099895888628153 | 4361 |
| 6 | 0.000899754702700977 | 3258 |

The decimal table is rounded upward for display. Authoritative upper endpoints
are the exact dyadic values in the JSON files. The sum is less than 0.0146.
The corrected evaluator SHA256 is
`bef7c1b73ad31e21cc56f0400cbf75aefc20807efd479b4a1e46a1ecce3443b5`.

The independent reviewer checked the polynomial division, PV factors,
crossing numerator, normalized-support inequality, threshold logarithm bound,
and UV remainder. A full fixed-partition replay at Arb288 also passed for all
three spins, including exact partition coverage and all dependency/table hashes.
The replay output files append `_replay.json` to the certificate stem. Run:

    & $py research/replay_lowspin_certificate.py output/certification_head/spin2_tight_component.json
    & $py research/replay_lowspin_certificate.py output/certification_head/spin4_tight_component.json
    & $py research/replay_lowspin_certificate.py output/certification_head/spin6_tight_component.json

The replayed uppers were respectively 0.01159999441097389514,
0.00209989588870513151, and 0.00089975470278891457. Differences in the last
few decimal places arise from interval radius rounding; all are independent
valid upper bounds below the chosen rational targets .0116, .0021, .0009.
