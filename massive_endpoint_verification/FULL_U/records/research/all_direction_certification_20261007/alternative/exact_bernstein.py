"""Exact FLINT polynomial arithmetic for rational Bernstein matrix certificates.

This supplies a third extraction implementation: affine polynomial composition
followed by one exact factorial convolution. No floating operation accepts PSD.
"""
from pathlib import Path
from fractions import Fraction
from math import factorial, comb
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tmp/certification_deps'))
from flint import fmpq, fmpq_poly

def q(value):
    if isinstance(value, fmpq):return value
    value=Fraction(value)
    return fmpq(value.numerator,value.denominator)

def fractions(poly,n=None):
    if n is None:n=len(poly)
    return [Fraction(str(poly[i])) for i in range(n)]

def bernstein(p,lo,hi):
    """Exact degree-len(p)-1 Bernstein coefficients on [lo,hi]."""
    n=len(p)-1
    lo,hi=q(lo),q(hi)
    affine=fmpq_poly(list(map(q,p)))(fmpq_poly([lo,hi-lo]))
    scaled=fmpq_poly([affine[k]*factorial(n-k) for k in range(n+1)])
    inverse=fmpq_poly([fmpq(1,factorial(k)) for k in range(n+1)])
    convolution=scaled*inverse
    return [convolution[j]*fmpq(factorial(j),factorial(n)) for j in range(n+1)]

def verify_matrix(arr,lo,hi,H):
    A,B,C,D=arr;H=q(H)
    aa=bernstein([q(a)*H+q(b) for a,b in zip(A,B)],lo,hi)
    cc=bernstein(C,lo,hi);dd=bernstein([q(d)*H for d in D],lo,hi)
    return all(a>=0 and d>=0 and a*d>=c*c for a,c,d in zip(aa,cc,dd))

def self_test():
    F=Fraction
    for n in range(15):
        p=[F((-1)**k*(k+2),k+1) for k in range(n+1)]
        for lo,hi in [(F(0),F(1)),(F(1,7),F(3,5)),(F(31,32),F(1))]:
            affine=[sum(p[j]*comb(j,k)*lo**(j-k)*(hi-lo)**k for j in range(k,n+1)) for k in range(n+1)]
            expected=[sum(affine[k]*F(comb(j,k),comb(n,k)) for k in range(j+1)) for j in range(n+1)]
            assert list(map(str,bernstein(p,lo,hi)))==list(map(str,expected))
    return True

if __name__=='__main__':
    print('PASS' if self_test() else 'FAIL')
