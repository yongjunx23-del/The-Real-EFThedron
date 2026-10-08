"""Exact all-even-spin Gram majorants for real kernels with derivatives.

For U_l=sum_{q=0}^d a_q [l(l+1)]^q, returns a matrix G such that
sum_even_l>L (2l+1) U_l^2/P_l(1+8/(3 delta)) <= |a|^T G |a|.
This is a component lemma, not a bound on a coupling or an energy integral.
"""
from __future__ import annotations
import argparse, json, math
from fractions import Fraction as F
from pathlib import Path


def reference_beta_integral(m, order, y):
    """Exact integral_0^y h^(order-m-2) (1-h)^m dh."""
    return sum((F((-1)**j*math.comb(m,j), order-m-1+j)*y**(order-m-1+j)
                for j in range(m+1)), F(0))


def legendre_reference_parameter(L, order):
    if L < order or L % 2:
        raise ValueError('Even L >= polynomial denominator order required')
    Tmin = (L+2)*(L+3)
    coefficients = [F(1)]
    for k in range(1, order+1):
        coefficients.append(coefficients[-1]*F(Tmin-(k-1)*k,Tmin*k*k)*F(4,3))
    # Exact rational choice; no floating roots or rounded comparisons.
    def passes(a):
        return all(F(math.comb(order,k))*a**k <= coefficients[k]
                   for k in range(order+1))
    a = F(1,4096)
    while passes(a+F(1,4096)):
        a += F(1,4096)
    if not passes(a):
        raise RuntimeError('Reference positivity failed')
    return a, coefficients


def gram_matrix(delta, L=40, degree=3, order=None):
    """Return exact rational G and its lower-denominator certificate."""
    delta = F(delta)
    if delta <= 0:
        raise ValueError('mu-4 must be positive')
    order = 2*degree+2 if order is None else order
    if order <= 2*degree+1:
        raise ValueError('Insufficient denominator degree for the infinite sum')
    a, coefficients = legendre_reference_parameter(L,order)
    TL = L*(L+1)
    ratio = F((L+2)*(L+3),TL)
    factor = F(2*L+5, 2*(2*L+1))
    y = 1/(1+a*TL/delta)
    moments = [factor*ratio**m*(delta/a)**(m+1)*reference_beta_integral(m,order,y)
               for m in range(2*degree+1)]
    return [[moments[i+j] for j in range(degree+1)] for i in range(degree+1)], {
        'L': L, 'degree': degree, 'denominator_order': order,
        'a': str(a), 'legendre_scaled_coefficients': list(map(str,coefficients)),
        'lower_bound': 'P_l(1+8/(3 delta)) >= (1+a*l*(l+1)/delta)^order, all even l>L',
        'quadratic_form_uses_absolute_monomial_coefficients': True,
    }


def discrete_beta_difference(arg):
    """digamma(arg)+EulerGamma+(0 or 2 log2), as an exact rational."""
    arg=F(arg)
    if arg.denominator == 1:
        return sum((F(1,j) for j in range(1,int(arg))), F(0))
    if arg.denominator != 2 or arg <= 0:
        raise ValueError('Positive integer/half-integer required')
    n=(arg.numerator-1)//2
    return 2*sum((F(1,j) for j in range(1,2*n+1)),F(0))-sum((F(1,j) for j in range(1,n+1)),F(0))


def exact_taylor_gram(L=40, degree=3, order=8):
    """Exact sum G_qr=sum_even_l>L(2l+1) B_q(l)B_r(l)/B_order(l).

    Useful on compact energy intervals. Using this bound alone in the far UV
    can fail integrability; use gram_matrix there.
    """
    if L%2 or L+2 < order or order <= 2*degree+1:
        raise ValueError('Insufficient spin cutoff or denominator order')
    roots=[j*(j+1) for j in range(order)]
    n0=L//2+1
    result=[]
    for q in range(degree+1):
        row=[]
        for r in range(degree+1):
            residues=[]
            for i,root in enumerate(roots):
                numerator=math.factorial(order)**2*math.prod(root-roots[j] for j in range(q))*math.prod(root-roots[j] for j in range(r))
                denominator=(math.factorial(q)*math.factorial(r))**2*math.prod(root-other for j,other in enumerate(roots) if i!=j)
                residues.append(F(numerator,denominator))
            assert sum(residues,F(0)) == 0
            value=-sum((a*(discrete_beta_difference(F(n0)-F(i,2))+discrete_beta_difference(F(n0)+F(i+1,2)))/2
                        for i,a in enumerate(residues)),F(0))
            assert value>0
            row.append(value)
        result.append(row)
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    import numpy as np
    from scipy.special import eval_legendre
    L=40;degree=3
    compact=exact_taylor_gram(L,degree)
    examples=[]
    for delta in [F(1),F(10),F(100),F(1000),F(10000)]:
        G,certificate=gram_matrix(delta,L,degree)
        # Example with exactly controlled polynomial real kernel.
        coeff=[F(1),-1/delta,F(1,10)/delta**2,-F(1,1000)/delta**3]
        upper=sum((abs(coeff[i]*coeff[j])*G[i][j] for i in range(4) for j in range(4)),F(0))
        spins=np.arange(L+2,2002,2);T=spins*(spins+1)
        U=sum(float(a)*T.astype(float)**i for i,a in enumerate(coeff))
        partial=float(np.sum((2*spins+1)*U**2/eval_legendre(spins,1+8/(3*float(delta)))))
        assert partial <= float(upper)
        examples.append({'delta':str(delta),'finite_spin_partial_sum':partial,'exact_all_spin_upper':str(upper),
                         'upper_approx':float(upper),'reference_certificate':certificate})
    output={'scope':'Exact all-spin real-kernel Gram component lemma; no energy integral or coupling bound yet',
            'compact_reference_matrix':[[str(x) for x in row] for row in compact],
            'examples':examples}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2),encoding='utf8')
    print(json.dumps({'scope':output['scope'],'examples':[{'delta':e['delta'],'partial':e['finite_spin_partial_sum'],'upper':e['upper_approx']} for e in examples]},indent=2))


if __name__=='__main__':
    main()
