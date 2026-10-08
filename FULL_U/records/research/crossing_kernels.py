"""Massive scalar crossing nulls from exact coefficient elimination.

For each total degree D use v=nu**2, w=t-4/3, and the invariant generators
X=2*v+3*w**2/2 and Y=-w*v+w**3/4. Coefficient rows with v**n, n>=1,
are dispersive. Rational vectors orthogonal to all invariant polynomials
give crossing nulls. The returned v kernels multiply Im f_l in the measure
dmu_l=(2*l+1)*dmu/(2*pi*rho). There is no real-part or subtraction term.

No quadrature, finite spin sum, or numerical search here proves a bound.
"""
from __future__ import annotations
from functools import lru_cache
from math import factorial, comb, gcd
import json
from pathlib import Path
import numpy as np
from scipy.special import eval_jacobi, poch, eval_legendre, lqmn
import sympy as sy

C=4./3.


@lru_cache(None)
def crossing_catalogue(max_degree=8):
    """Return immutable canonical rational null metadata, through max_degree."""
    if max_degree<4: return ()
    v,w=sy.symbols('v w')
    X=2*v+sy.Rational(3,2)*w**2
    Y=-w*v+sy.Rational(1,4)*w**3
    rows=[]
    for D in range(4,max_degree+1):
        indices=tuple((n,D-2*n) for n in range(1,D//2+1))
        invariants=tuple((a,b) for b in range(D//3+1) for a in range(D//2+1) if 2*a+3*b==D)
        basis=[]
        for a,b in invariants:
            polynomial=sy.Poly(sy.expand(X**a*Y**b),v,w)
            basis.append([polynomial.coeff_monomial(v**n*w**k) for n,k in indices])
        matrix=sy.Matrix(basis)
        nulls=matrix.nullspace()
        for number,vector in enumerate(nulls):
            denominator=sy.ilcm(*[q.q for q in vector])
            integers=[int(q*denominator) for q in vector]
            common=0
            for q in integers: common=gcd(common,abs(q))
            integers=[q//common for q in integers]
            if next(q for q in integers if q)!=abs(next(q for q in integers if q)):
                integers=[-q for q in integers]
            weights=tuple(sy.Rational(q) for q in integers)
            # Primitive first relation [2,-3] yields four times Eq.(A5.6).
            if D==4: weights=tuple(q/4 for q in weights)
            assert matrix*sy.Matrix(weights)==sy.zeros(len(invariants),1)
            # Exact combination BEFORE evaluation cancels the scalar P_l term.
            # beta_j multiplies (l+1)_j/j! * P_(l-j)^(j,j)(Z)
            #                       / [(mu-4)^j Q^(D+1-j)].
            beta=[]
            for j in range(D-1):
                coefficient=sy.Rational(0)
                for weight,(n,k) in zip(weights,indices):
                    if k>=j:
                        coefficient += 2*weight*(-1)**(k-j)*sy.binomial(D-j,k-j)/2**(k-j)
                beta.append(coefficient)
            assert beta[0]==0, (D,weights,beta)
            name=f'X{D}' if len(nulls)==1 else f'X{D}_{number+1}'
            rows.append({'name':name,'degree':D,'coefficient_indices':indices,
                         'invariant_powers':invariants,'weights':weights,'beta':tuple(beta)})
    return tuple(rows)


def build_crossing_bank(mu,spins,max_degree=8,normalized=False):
    """Return {labels,u,v,g0_coeff,metadata}, arrays (rows,spins,mu).

    normalized=False: raw coefficients of r,b. normalized=True: v/p, the
    coefficients of B=p*b; u remains zero. Use no additional spin measure.
    The first row is exactly the existing X4 normalization.
    """
    mu=np.asarray(mu,dtype=float);spins=np.asarray(spins,dtype=int)
    if mu.ndim!=1 or np.any(mu<=4) or np.any(~np.isfinite(mu)):
        raise ValueError('mu must be finite, one-dimensional and strictly >4')
    if spins.ndim!=1 or np.any(spins<0) or np.any(spins%2):
        raise ValueError('Physical spins must be nonnegative even integers')
    catalogue=crossing_catalogue(max_degree)
    d=mu-4.;Q=mu-C;Z=Q/d
    values=np.zeros((len(catalogue),len(spins),len(mu)))
    legendre=np.array([eval_legendre(int(l),Z) for l in spins]) if normalized else None
    for row,entry in enumerate(catalogue):
        D=entry['degree']
        for j,beta in enumerate(entry['beta']):
            if beta==0: continue
            for il,l in enumerate(spins):
                if l<j: continue
                # This is the Jacobi form of d^j P_l/dZ^j, with the
                # 2^-j canceled by differentiating Z+2*w/d.
                jacobi=eval_jacobi(int(l)-j,j,j,Z)
                if normalized:
                    # Form the ratio before multiplying threshold powers;
                    # a raw high-spin kernel may overflow even if v/p is small.
                    term=(float(beta)*poch(int(l)+1,j)/factorial(j)
                          *(jacobi/legendre[il])/d**j/Q**(D-2-j))
                else:
                    term=(float(beta)*poch(int(l)+1,j)/factorial(j)
                          *jacobi/d**j/Q**(D+1-j))
                values[row,il]+=term
    metadata=[{key:([str(x) for x in val] if key in ('weights','beta') else val)
               for key,val in entry.items()} for entry in catalogue]
    return {'labels':np.array([entry['name'] for entry in catalogue]),
            'u':np.zeros_like(values),'v':values,'g0_coeff':np.zeros(len(catalogue)),
            'metadata':metadata}


def crossing_rows(mu,spins,max_degree=8):
    """Convenience list of (name,u,v) ready for massive_search rows."""
    bank=build_crossing_bank(mu,spins,max_degree,False)
    return list(zip(bank['labels'],bank['u'],bank['v']))


def _self_test():
    from massive_search import first_crossing
    assert len(crossing_catalogue(8))==7
    assert len(crossing_catalogue(16))==35
    mu=np.array([4.001,4.3,5.,8.,20.,1000.])
    spins=np.arange(0,14,2)
    bank=build_crossing_bank(mu,spins)
    baseline=first_crossing(mu,spins)
    relative=np.max(np.abs(bank['v'][0]-baseline)/(1+abs(baseline)))
    np.testing.assert_allclose(bank['v'][0],baseline,rtol=3e-13,atol=2e-13)
    assert np.max(abs(bank['v'][:,0]))==0
    p=np.array([eval_legendre(int(l),(mu-C)/(mu-4)) for l in spins])/(mu-C)**3
    np.testing.assert_allclose(build_crossing_bank(mu,spins,normalized=True)['v'],
                               bank['v']/p[None,:,:],rtol=2e-12,atol=1e-12)
    # Independent direct polynomial-in-w coefficient convolution, without the
    # symbolic precombination used by build_crossing_bank.
    point=8.;degree=8;l=10;d=point-4.;Q=point-C;Z=Q/d
    direct=[]
    for entry in crossing_catalogue(degree):
        value=0.
        for weight,(n,k) in zip(entry['weights'],entry['coefficient_indices']):
            one=0.
            for j in range(k+1):
                if j>l: continue
                pcoef=poch(l+1,j)*eval_jacobi(l-j,j,j,Z)/factorial(j)/d**j
                dcoef=(-1)**(k-j)*comb(2*n+k-j,k-j)/2**(k-j)/Q**(2*n+1+k-j)
                one += 2*pcoef*dcoef
            value += float(weight)*one
        direct.append(value)
    evaluated=build_crossing_bank([point],[l])['v'][:,0,0]
    np.testing.assert_allclose(direct,evaluated,rtol=1e-12,atol=1e-12)
    # Nontrivial crossing test: A=1/[(M-s)(M-t)(M-u)]. This is symmetric
    # and decreases as 1/s^2 at fixed t. Its s-pole residue has infinitely
    # many even partial waves. Their exact coefficients use Legendre Q_l.
    # An isolated spin>=2 exchange would fail the two-subtraction falloff
    # and therefore is deliberately NOT used as a zero-null test.
    pole_errors=[]
    for mass_squared in [6.,8.,20.]:
        d=mass_squared-4.;AA=mass_squared+d/2;BB=d/2
        pole_spins=np.arange(0,42,2)
        q=lqmn(0,int(pole_spins[-1]),AA/BB)[0][0,pole_spins]
        spectral_coeff=(2*pole_spins+1)*q/(AA*BB)
        kernels=build_crossing_bank([mass_squared],pole_spins)['v'][:,:,0]
        nulls=kernels@spectral_coeff
        relative_pole=np.max(abs(nulls)/(np.abs(kernels)@spectral_coeff+1e-100))
        assert relative_pole<2e-10,(mass_squared,relative_pole,nulls)
        pole_errors.append({'mass_squared':mass_squared,'relative_null_residual':float(relative_pole),
                            'maximum_absolute_null_residual':float(max(abs(nulls)))})
    result={'through_degree8_count':7,'through_degree16_count':35,
            'first_X4_relative_difference':float(relative),'spin_zero_exactly_zero':True,
            'pole_exchange_checks':pole_errors,'rows':bank['metadata']}
    Path(__file__).with_name('crossing_kernel_audit.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':
    _self_test()
