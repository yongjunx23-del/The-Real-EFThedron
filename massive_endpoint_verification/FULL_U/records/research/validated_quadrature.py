"""Arb-enclosed Gauss-Legendre rules, with independently bracketed roots.

SciPy supplies guesses only. A rule is accepted after all n disjoint intervals
strictly inside (-1,1) bracket a root of P_n. Degree n then proves completeness
and uniqueness. Weights are evaluated from the exact derivative formula.
The rule alone is NOT a quadrature-error certificate for nonpolynomials.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tmp/certification_deps'))
from flint import arb, ctx
from scipy.special import roots_legendre

_cache = {}

def legendre_pair(n, x):
    """P_n(x), P_(n-1)(x), in ball arithmetic."""
    # Arb's special-function implementation internally guards precision.
    # A naive ball recurrence unnecessarily amplifies rounding near +/-1.
    return x.legendre_p(n), x.legendre_p(n-1) if n else arb(0)

def legendre_derivative(n, x):
    p, q = legendre_pair(n, x)
    return n*(x*p-q)/(x*x-1)

def gauss_legendre(n, precision=None):
    """Return immutable tuples (node_balls, weight_balls), certified on [-1,1]."""
    if precision is not None and precision != ctx.prec:
        raise ValueError('Set flint.ctx.prec before constructing/using the rule')
    if n < 1:
        raise ValueError('n must be positive')
    key = (n, ctx.prec)
    if key in _cache:
        return _cache[key]
    roots = []
    radius = arb(2)**(-ctx.prec+30)
    for guess in roots_legendre(n)[0]:
        x = arb(float(guess))
        for _ in range(8):
            p, q = legendre_pair(n, x)
            dp = n*(x*p-q)/(x*x-1)
            x = (x-p/dp).mid()
        rad = radius
        for _ in range(32):
            a, b = (x-rad).lower(), (x+rad).upper()
            pa, pb = legendre_pair(n,a)[0], legendre_pair(n,b)[0]
            if pa*pb < 0:
                root = arb((a+b)/2, (b-a)/2)
                assert root > -1 and root < 1
                if roots:
                    assert roots[-1].upper() < root.lower()
                roots.append(root)
                break
            rad *= 2
        else:
            raise ArithmeticError('Failed to sign-bracket a Gauss root')
    weights = tuple(2/((1-x*x)*legendre_derivative(n,x)**2) for x in roots)
    assert all(w > 0 for w in weights)
    assert sum(weights, arb(0)).contains(2)
    result = (tuple(roots), weights)
    _cache[key] = result
    return result

def analytic_error_bound(n, rho, magnitude, interval_length=2):
    """GL error from a Bernstein-ellipse bound |f|<=magnitude, rho>1.

    Chebyshev coefficients satisfy |a_k|<=2 M rho^-k; exactness below
    degree 2n and positivity of weights give |I(T_k)-Q(T_k)|<=4.
    Map an interval of the supplied length to [-1,1].
    Caller must independently certify the ellipse analyticity and magnitude.
    """
    rho, magnitude = arb(rho), arb(magnitude)
    assert rho > 1 and magnitude >= 0
    return 4*interval_length*magnitude/rho**(2*n)/(1-1/rho)

def geometric_remainder_bound(nterms, q, numerator_bound, center_denominator,
                              domain_mass=1):
    """Bound integral MINUS positive quadrature via two remainder bounds.

    For numerator/(d-e), certify |e/d|<=q<1 and |numerator|<=N on
    the real source domain, d>=center_denominator>0. Truncate the geometric
    series before q**nterms. The caller MUST check that the truncated source
    polynomial degree is within exactness of every tensor Gauss rule.
    domain_mass is the integral of the positive reference source measure.
    """
    q, den = arb(q), arb(center_denominator)
    assert q >= 0 and q < 1 and den > 0 and nterms >= 0
    return 2*domain_mass*numerator_bound*q**nterms/(den*(1-q))

def self_test():
    ctx.prec = 160
    records = []
    for n in (1, 2, 8, 24, 48):
        x, w = gauss_legendre(n)
        powers = [arb(1) for _ in x]
        for k in range(2*n):
            value = sum((wi*pi for wi,pi in zip(w,powers)),arb(0))
            exact = arb(2)/(k+1) if k%2 == 0 else arb(0)
            assert value.overlaps(exact), (n,k,value,exact)
            powers = [pi*xi for pi,xi in zip(powers,x)]
        records.append({'order':n,'root_brackets':len(x),'moments_verified':2*n})
    # Nonpolynomial example with an exact geometric remainder bound.
    n = 8
    x,w = gauss_legendre(n)
    approx = sum((wi/(3-xi) for xi,wi in zip(x,w)),arb(0))
    error = geometric_remainder_bound(2*n, arb(1)/3, arb(1), arb(3), 2)
    exact = arb(2).log()
    assert arb(approx, error.upper()).overlaps(exact)
    return {'precision_bits':ctx.prec,'rules':records,'geometric_example_passed':True}

if __name__ == '__main__':
    import json
    print(json.dumps(self_test(),indent=2))
