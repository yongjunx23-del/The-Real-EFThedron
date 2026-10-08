"""All-energy fixed-transfer strip rows extending finite-energy Roy rows.

Public raw-default wrapper around the cancellation-safe implementation.
No finite-spin construction by itself is an outer exclusion bound.
"""
from roy_kernels import build_fixed_transfer_bank


def build_tapered_transfer_bank(mu,spins,specs,power=1,normalized=False):
    """Use endpoint-vanishing (1-x^2)^power P_n(x) energy tests.

    Exact finite Legendre combinations; removes logarithmic endpoint poles for
    power>=1. This is a restriction of valid tests, not a new physical axiom.
    """
    import numpy as np
    from numpy.polynomial.legendre import legmul
    full=[]; combinations=[]
    for A,B,T,j,n in specs:
        coeff=np.eye(1,n+1,n)[0]
        for _ in range(power): coeff=legmul(coeff,[2/3,0.,-2/3])
        terms=[]
        for m,value in enumerate(coeff):
            if value:
                spec=(A,B,T,j,m)
                if spec not in full: full.append(spec)
                terms.append((full.index(spec),value))
        combinations.append(terms)
    bank=build_transfer_bank(mu,spins,full,normalized=normalized)
    transform=np.zeros((len(specs),len(full)))
    for i,terms in enumerate(combinations):
        for k,v in terms: transform[i,k]=v
    return dict(u=np.einsum('ab,bij->aij',transform,bank['u']),
                v=np.einsum('ab,bij->aij',transform,bank['v']),
                g0_coeff=transform@bank['g0_coeff'],p=bank['p'],specs=specs,
                labels=[f"Taper{power}_Strip_A{A:g}_B{B:g}_T{T:g}_j{j}_n{n}"
                        for A,B,T,j,n in specs],normalized=normalized)


def build_transfer_bank(mu, spins, specs=None, energy_order=40,
                        transfer_order=32, normalized=False):
    """Rows obey <u*r+v*b>+g0_coeff*g0=0 unless normalized=True.

    Default bank: windows [12,24], [24,48], t in [-4/3,0],
    Legendre energy degrees n=0..5 and transfer degrees j=0..2.
    Explicit specs consist of (A,B,T,j,n) tuples.
    Output u,v shape (nrows,nspins,nmu), no measure/quadrature weights included.
    Set normalized=True for coefficients of R=sqrt(p)r and B=p*b.
    """
    if specs is None:
        specs=[(A,B,4/3,j,n) for A,B in [(12,24),(24,48)]
               for j in range(3) for n in range(6)]
    return build_fixed_transfer_bank(mu,spins,specs,energy_order,
                                    transfer_order,normalized)


def _self_test():
    import numpy as np
    from scipy.integrate import quad
    from scipy.special import eval_legendre
    from numpy.polynomial.legendre import leggauss
    from roy_kernels import real_legendre_q, C
    A,B,T=12.,24.,4/3
    mus=np.array([4.01,8.,11.7,12.1,17.13,23.8,24.5,26.,100.,1e5])
    specs=[(A,B,T,j,n) for j in range(3) for n in range(4)]
    bank=build_transfer_bank(mus,[0],specs,70,64)
    te,tw=leggauss(90)
    t=-T+T*(te+1)/2
    for k,(_,_,_,j,n) in enumerate(specs):
        hw=tw/2*eval_legendre(j,te)
        for im,m in enumerate(mus[:-1]):
            def regular(s):
                u=4-s-t
                return np.sum(hw*(1/(m-t)+1/(m-u)-3/(m-C)))
            w=lambda s:eval_legendre(n,2*(s-A)/(B-A)-1)/(B-A)
            val=quad(lambda s:w(s)*regular(s),A,B,epsabs=2e-11)[0]
            if j==0:
                val+=2*float(real_legendre_q(n,2*(m-A)/(B-A)-1))/(B-A)
            np.testing.assert_allclose(bank['v'][k,0,im],-val,atol=3e-11,rtol=2e-8)
            local=np.pi*w(m) if A<m<B and j==0 else 0
            np.testing.assert_allclose(bank['u'][k,0,im],local,atol=2e-12,rtol=1e-10)
    # Genuine higher spins and threshold/UV span: independent quadrature orders.
    spins=range(0,22,2)
    low=build_transfer_bank(mus,spins,specs,40,32,True)
    high=build_transfer_bank(mus,spins,specs,80,70,True)
    err=np.max(np.abs(low['v']-high['v'])/(1+np.abs(high['v'])))
    assert err<2e-9,err
    print(f"transfer_kernels self-test passed: independent scalar PV; l<=20 refinement {err:.3g}")


if __name__=='__main__':
    _self_test()
