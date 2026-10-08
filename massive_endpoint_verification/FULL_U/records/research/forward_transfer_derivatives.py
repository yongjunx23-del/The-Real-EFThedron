"""Threshold-regular forward transfer derivatives of massive harmonic nulls.

Imposing raw nonforward harmonics in disk support can diverge at high spin.
These derivative combinations insert (1-omega)^(2p-2) BEFORE differentiation,
as in manuscript A5.7. Derivatives are evaluated by exact truncated Taylor
algebra, not finite differences. No numeric cutoff is introduced.
"""
from math import factorial,comb
import numpy as np
from scipy.special import eval_legendre


def mul(a,b):
    n=len(a)-1
    return np.array([sum((a[k]*b[j-k] for k in range(j+1)),np.zeros_like(a[0]*b[0]))
                     for j in range(n+1)])


def power(a,p):
    out=np.zeros_like(a);out[0]=1
    while p:
        if p%2:out=mul(out,a)
        p//=2
        if p:a=mul(a,a)
    return out


def build_derivative_bank(s,ell,max_order=3,hmax=20):
    s=np.asarray(s,dtype=float)[None,:];ell=np.asarray(ell,dtype=int)[:,None]
    if max_order not in (1,2,3):raise ValueError('Implemented orders1..3.')
    shape=np.broadcast_shapes(s.shape,ell.shape);N=max_order
    a=np.zeros((N+1,)+shape,dtype=complex);a[0]=2.;a[1]=.5
    q=s-2;delta=s-4;Q=s-4/3
    b=np.zeros_like(a);b0=np.sqrt(s*delta)
    bc=1.
    for k in range(N+1):
        if k:bc*=((.5-(k-1))/k)
        b[k]=b0*bc/s**k
    qp={}
    for pp in [1,2,3]:
        qp[pp]=np.array([np.broadcast_to((-1)**k*comb(pp+k-1,k)*(.5/q)**k/q**pp,shape)
                         for k in range(N+1)],dtype=complex)
    numerator=a+1j*b
    omega=mul(power(numerator,2),qp[2])
    # Stable exact 1-omega=-2i sqrt((s-4)(s+t)) (A+i sqrt(...))/Q_t^2.
    difference=-2j*mul(mul(b,numerator),qp[2])
    F=8*mul(power(a,2),qp[3])
    Pt=np.zeros_like(a)
    for il,l in enumerate(ell[:,0]):
        for k in range(min(int(l),N)+1):
            Pt[k,il]=comb(int(l),k)*comb(int(l)+k,k)/delta[0]**k
    Pstar=eval_legendre(ell,Q/delta)
    Ka=np.array([2*qp[1][k]-np.broadcast_to((-1)**k/(s-8/3)**(k+1),shape)
                 for k in range(N+1)],dtype=complex)
    Kb=np.array([np.broadcast_to(1/s**(k+1)+(-1)**k/(s-8/3)**(k+1),shape)
                 for k in range(N+1)],dtype=complex)
    Ka[0]-=1/Q;Kb[0]-=2/Q
    sub=mul(Pt,Ka)+Pstar[None,:,:]*Kb
    base=mul(F,Pt)
    rows=[]
    for p in range(1,N+1):
        diffpow=power(difference,2*p-2)
        for j in range(2,hmax+1):
            qjet=mul(mul(base,power(omega,j)),diffpow)
            val=factorial(p)*qjet[p]
            rows.append((f'Dt{p}H_j{j}',val.imag,val.real))
        qjet=mul(mul(base,omega),diffpow)
        val=factorial(p)*qjet[p]
        vv=val.real-factorial(p)*sub[p].real
        rows.append((f'Dt{p}G_regular',val.imag,vv))
    return rows


def _self_test():
    import mpmath as mp
    mp.mp.dps=65
    ss=np.array([4.01,5.,12.,100.]);ells=np.array([0,2,8])
    rows=build_derivative_bank(ss,ells,3,4)
    lookup={name:(u,v) for name,u,v in rows}
    worst=0.
    for il,l in enumerate(ells):
        for im,s0 in enumerate(ss):
            sm=mp.mpf(str(s0));c=mp.mpf(4)/3;Q=sm-c
            ps=mp.legendre(int(l),Q/(sm-4))
            def val(t,j,p,g):
                A=2+t/2;Qt=sm-2+t/2;bb=mp.sqrt((sm-4)*(sm+t))
                om=((A+1j*bb)/Qt)**2;F=8*A*A/Qt**3
                pt=mp.legendre(int(l),1+2*t/(sm-4))
                z=F*pt*om**j*(1-om)**(2*p-2)
                if g:
                    ka=2/Qt-1/Q-1/(sm-mp.mpf(8)/3+t)
                    kb=1/(sm-t)+1/(sm-mp.mpf(8)/3+t)-2/Q
                    z-=pt*ka+ps*kb
                return z
            for p in range(1,4):
                for j in [2,4]:
                    ref=complex(mp.diff(lambda tt:val(tt,j,p,False),mp.mpf(0),p))
                    u,v=lookup[f'Dt{p}H_j{j}'];got=v[il,im]+1j*u[il,im]
                    worst=max(worst,abs(ref-got)/(1+abs(ref)))
                    np.testing.assert_allclose(got,ref,rtol=3e-11,atol=2e-9)
                ref=complex(mp.diff(lambda tt:val(tt,1,p,True),mp.mpf(0),p))
                u,v=lookup[f'Dt{p}G_regular'];got=v[il,im]+1j*u[il,im]
                worst=max(worst,abs(ref-got)/(1+abs(ref)))
                np.testing.assert_allclose(got,ref,rtol=3e-10,atol=2e-9)
    print('Forward transfer derivative validation max relative-scale error',worst)


if __name__=='__main__':_self_test()
