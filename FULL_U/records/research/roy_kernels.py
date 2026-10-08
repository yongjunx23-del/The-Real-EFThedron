"""Smeared massive Roy equations, with both removable poles cancelled.

Research numerical kernel builder, not an all-spin/interval certificate.
Normalization: dmu_l=(2l+1) ds/(2*pi*rho), f_l=r_l+i*b_l.
Rows obey <u*r+v*b> + g0_coeff*g0 = 0.
With normalized=True, rows instead multiply R=sqrt(p)*r and B=p*b.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_legendre, lqmn

C = 4.0 / 3.0


def legendre_divided_all(x, y, lmax):
    """Return P_l(x), P_l(y), and (P_l(x)-P_l(y))/(x-y), l=0..lmax.

    Recurrence is exact also at x=y; no cancellation/division by x-y occurs.
    This routine is intended for modest l and finite inputs.
    """
    x, y = np.broadcast_arrays(x, y)
    px = np.empty((lmax + 1,) + x.shape, dtype=float)
    py = np.empty_like(px)
    dd = np.empty_like(px)
    px[0], py[0], dd[0] = 1.0, 1.0, 0.0
    if lmax:
        px[1], py[1], dd[1] = x, y, 1.0
    for l in range(1, lmax):
        a, b = (2 * l + 1) / (l + 1), l / (l + 1)
        px[l + 1] = a * x * px[l] - b * px[l - 1]
        py[l + 1] = a * y * py[l] - b * py[l - 1]
        dd[l + 1] = a * (x * dd[l] + py[l]) - b * dd[l - 1]
    return px, py, dd


def real_legendre_q(n, y):
    """Real PV branch Q_n(y)=.5*PV int_-1^1 P_n(x)/(y-x) dx."""
    y = np.asarray(y, dtype=float)
    if np.any(np.abs(y) == 1):
        raise ValueError("Spectral node equals external smearing endpoint.")
    out = np.empty_like(y)
    inside = np.abs(y) < 1
    if np.any(inside):
        z = y[inside]
        q0 = 0.5 * (np.log1p(z) - np.log1p(-z))
        if n == 0:
            out[inside] = q0
        else:
            q1 = z * q0 - 1
            for k in range(1, n):
                q0, q1 = q1, ((2*k+1)*z*q1-k*q0)/(k+1)
            out[inside] = q1
    outside = ~inside
    if np.any(outside):
        z = y[outside]
        vals = lqmn(0, n, np.abs(z))[0][0, n]
        out[outside] = np.where(z < 0, (-1) ** (n + 1) * vals, vals)
    return out


def scalar_kernel_stable(mu, s, t):
    """l=0 combined crossed-subtraction kernel; avoids UV cancellation."""
    u = 4.0 - s - t
    q = mu - C
    return ((s-C)**2/(mu-s) + (t-C)**2/(mu-t)
            + (u-C)**2/(mu-u)) / q**2


def build_roy_bank(mu, spins, specs, energy_order=48, angular_order=40,
                   normalized=True):
    """Build folded even-J rows for (B,J,n) specifications.

    w_n(s)=P_n(2(s-4)/(B-4)-1)/(B-4), 4<s<B.
    For the currently justified numerical family require 4<B<12.
    All supplied spins and J must be nonnegative even integers.
    Keys u,v have shape (nrows,nspins,nmu). g0_coeff has shape (nrows,).
    p has shape (nspins,nmu). No quadrature weights are included.

    mu is the complete caller-specified spectral grid. This function does NOT
    discard/control additional spins or energies and thus alone proves no bound.
    """
    mu = np.asarray(mu, dtype=float)
    spins = np.asarray(spins, dtype=int)
    if mu.ndim != 1 or np.any(mu <= 4) or np.any(~np.isfinite(mu)):
        raise ValueError("mu must be finite, one-dimensional, and strictly >4.")
    if spins.ndim != 1 or np.any(spins < 0) or np.any(spins % 2):
        raise ValueError("Only nonnegative even spins are supported.")
    specs = [(float(B), int(J), int(n)) for B, J, n in specs]
    if any(not (4 < B < 12) or J < 0 or J % 2 or n < 0 for B,J,n in specs):
        raise ValueError("Require 4<B<12, even J>=0, and n>=0.")
    lmax = int(max(spins, default=0))
    delta = mu - 4.0
    q = mu - C
    Z = q / delta
    pz = np.array([eval_legendre(int(l), Z) for l in spins])
    p = pz / q[None, :]**3
    ur = np.zeros((len(specs), len(spins), len(mu)))
    vr = np.zeros_like(ur)
    gc = np.zeros(len(specs))
    ze, zw = leggauss(angular_order)
    z = (ze + 1.0) / 2.0
    # zw already includes the folded angular factor 2 dz.
    se, sw = leggauss(energy_order)
    for B in sorted(set(row[0] for row in specs)):
        indexes = [k for k,row in enumerate(specs) if row[0] == B]
        s = 4.0 + (B - 4.0) * (se + 1.0) / 2.0
        ss = s[:, None]
        zz = z[None, :]
        t = -(ss-4.0) * (1.0-zz) / 2.0
        u = 4.0-ss-t
        weights = {k: (sw/2.0 * eval_legendre(specs[k][2], se))[:, None]
                   * (zw * eval_legendre(specs[k][1], z))[None, :]
                   for k in indexes}
        for im, m in enumerate(mu):
            d, qm, zm = delta[im], q[im], Z[im]
            y = 1.0 + 2.0*t/d
            pz_here = pz[:, im, None, None]
            if m >= B + 1.0:
                # Direct nonsingular formula factored to exhibit its UV falloff.
                # Here y+Z is bounded away from zero. Evaluating only the
                # requested Legendre orders avoids a 0..lmax recurrence when
                # a few large-spin UV guards are requested. Any cancellation
                # in P_l(y)/P_l(Z)-1 affects only a vanishing UV correction.
                py=np.array([eval_legendre(int(l),y) for l in spins])
                dd_over_pz=(py/pz_here-1.)/(y+zm)[None,:,:]
                base = scalar_kernel_stable(m, ss, t)
                correction = (2.0/d)*(ss-C)*(C-u)*(2*m+t-4.0)/(
                    (m-ss)*qm*(m-u))
                integrand_over_pz = base[None, :, :] + correction[None,:,:]*dd_over_pz
            else:
                py_all, _, dd_minus = legendre_divided_all(y, -zm, lmax)
                py = py_all[spins]
                ddm = dd_minus[spins]
                # Cancel the external PV pole analytically before integration:
                # P_l(y)/(m-s)=P_l(z)/(m-s)+(1-z)/d D_l(y,z).
                _, _, dd_phys_all = legendre_divided_all(y, zz, lmax)
                dd_phys = dd_phys_all[spins]
                regular = ((1.0-zz)[None,:,:]/d*dd_phys
                           +py*(1.0/(m-u)-1.0/qm)[None,:,:]
                           +pz_here*(1.0/(m-t)-2.0/qm)[None,:,:]
                           -2.0/d*ddm)
                integrand_over_pz = regular/pz_here
            for k in indexes:
                _, J, n = specs[k]
                val = np.einsum('lsz,sz->l', integrand_over_pz, weights[k])
                if m < B + 1.0 and J in spins:
                    il = int(np.flatnonzero(spins == J)[0])
                    yy = 2.0*(m-4.0)/(B-4.0)-1.0
                    pv = 4.0*float(real_legendre_q(n, yy))/((2*J+1)*(B-4.0))
                    val[il] += pv/pz[il, im]
                # Equation contains -D_l b_l. Divide by p=Pz/Q^3.
                vr[k, :, im] = -q[im]**3 * val
        for k in indexes:
            _, J, n = specs[k]
            if J in spins:
                il = int(np.flatnonzero(spins == J)[0])
                active = (mu > 4) & (mu < B)
                w = eval_legendre(n, 2*(mu[active]-4)/(B-4)-1)/(B-4)
                ur[k, il, active] = 2*np.pi*w/(2*J+1)/np.sqrt(p[il,active])
            gc[k] = -2.0 if J == 0 and n == 0 else 0.0
    if not normalized:
        ur *= np.sqrt(p)[None, :, :]
        vr *= p[None, :, :]
    return dict(u=ur, v=vr, g0_coeff=gc, p=p, specs=specs,
                labels=[f"Roy_B{B:g}_J{J}_n{n}" for B,J,n in specs],
                normalized=normalized)


def build_fixed_transfer_bank(mu, spins, specs, energy_order=40,
                             transfer_order=32, normalized=True):
    """High-energy fixed-transfer strip equations, without full Roy projection.

    specs=(A,B,T,j,n), external s in [A,B], t in [-T,0], with
    w_n=P_n(mapped s)/(B-A), h_j=P_j(mapped t)/T.
    Require A>4, B>A and 0<T<=min(4/3,A-4). This keeps every external point
    physical and every fixed transfer safely in the conventional analytic strip.
    Real-part kernels involve ALL even spins, unlike the finite-J Roy family.
    Thus an all-spin remainder must account for these real kernels too.

    The same u,v,g0_coeff,p output convention as build_roy_bank is used.
    """
    mu = np.asarray(mu,dtype=float)
    spins = np.asarray(spins,dtype=int)
    if mu.ndim != 1 or np.any(mu<=4) or np.any(~np.isfinite(mu)):
        raise ValueError("Require finite one-dimensional mu>4.")
    if spins.ndim != 1 or np.any(spins<0) or np.any(spins%2):
        raise ValueError("Require nonnegative even spins.")
    specs=[(float(A),float(B),float(T),int(j),int(n)) for A,B,T,j,n in specs]
    if any(A<=4 or B<=A or T<=0 or T>min(C,A-4) or j<0 or n<0
           for A,B,T,j,n in specs):
        raise ValueError("Require A>4, B>A, 0<T<=min(4/3,A-4), j,n>=0.")
    lmax=int(max(spins,default=0))
    delta=mu-4
    q=mu-C
    Z=q/delta
    pz=np.array([eval_legendre(int(l),Z) for l in spins])
    p=pz/q[None,:]**3
    ur=np.zeros((len(specs),len(spins),len(mu)))
    vr=np.zeros_like(ur)
    gc=np.zeros(len(specs))
    se,sw=leggauss(energy_order)
    te,tw=leggauss(transfer_order)
    for k,(A,B,T,j,n) in enumerate(specs):
        width=B-A
        ss=(A+width*(se+1)/2)[:,None]
        tt=(-T+T*(te+1)/2)[None,:]
        uext=4-ss-tt
        hw=tw/2*eval_legendre(j,te)
        ew=sw/2*eval_legendre(n,se)
        all_weights=ew[:,None]*hw[None,:]
        # H_l(v)=integral h_j(t) P_l(1+2t/(v-4)) dt.
        yy=1+2*tt/delta[:,None]
        H=np.array([eval_legendre(int(l),yy)@hw for l in spins])
        active=(mu>A)&(mu<B)
        w=eval_legendre(n,2*(mu[active]-A)/width-1)/width
        ur[k,:,active]=(np.pi*H[:,active]*w[None,:]/np.sqrt(p[:,active])).T
        gc[k]=-1.0 if n==0 and j==0 else 0.0
        for im,m in enumerate(mu):
            dm,qm,zm=delta[im],q[im],Z[im]
            y=1+2*tt/dm
            pz_here=pz[:,im,None,None]
            if m<A-1 or m>B+1:
                # In the UV y+Z>0, so only selected Legendre orders are needed.
                # Below the window retain the cancellation-safe recurrence.
                if m>B+1:
                    py=np.array([eval_legendre(int(l),y) for l in spins])
                    dd_over_pz=(py/pz_here-1.)/(y+zm)[None,:,:]
                else:
                    py_all,_,dd_all=legendre_divided_all(y,-zm,lmax)
                    dd_over_pz=dd_all[spins]/pz_here
                base=scalar_kernel_stable(m,ss,tt)
                corr=(2/dm)*(ss-C)*(C-uext)*(2*m+tt-4)/(
                    (m-ss)*qm*(m-uext))
                integrand=base[None,:,:]+corr[None,:,:]*dd_over_pz
                val=np.einsum('lst,st->l',integrand,all_weights)
            else:
                py_all,_,dd_all=legendre_divided_all(y,-zm,lmax)
                py=py_all[spins]
                ddm=dd_all[spins]
                I1=2*float(real_legendre_q(n,2*(m-A)/width-1))/width
                yc=2*(4-m-tt-A)/width-1
                I2=-2*real_legendre_q(n,yc)/width
                regular=py*I2[None,:,:]
                if n==0:
                    regular=regular-py/qm+pz_here*(1/(m-tt)-2/qm)[None,:,:]-2/dm*ddm
                val=(H[:,im]*I1+np.einsum('lnt,t->l',regular,hw))/pz[:,im]
            vr[k,:,im]=-qm**3*val
    if not normalized:
        ur*=np.sqrt(p)[None,:,:]
        vr*=p[None,:,:]
    return dict(u=ur,v=vr,g0_coeff=gc,p=p,specs=specs,
                labels=[f"Strip_A{A:g}_B{B:g}_T{T:g}_j{j}_n{n}"
                        for A,B,T,j,n in specs], normalized=normalized)


def _self_test():
    """Numerical independent identities; not tail/continuum certification."""
    from scipy.integrate import quad
    rng = np.random.default_rng(291)
    for l in (0, 2, 4, 10, 20):
        x = rng.uniform(-2, 2, 20)
        y = rng.uniform(-2, 2, 20)
        px, py, dd = legendre_divided_all(x, y, l)
        ref = (eval_legendre(l,x)-eval_legendre(l,y))/(x-y)
        np.testing.assert_allclose(dd[l], ref, rtol=2e-10, atol=1e-8)
        _, _, diag = legendre_divided_all(x,x,l)
        if l:
            deriv = l*(x*eval_legendre(l,x)-eval_legendre(l-1,x))/(x*x-1)
            np.testing.assert_allclose(diag[l], deriv, rtol=1e-10, atol=1e-8)
    for n in range(7):
        for yy in (-2.,-.73,.12,.95,1.2,5.,1000.):
            if abs(yy)<1:
                # scipy cauchy weight integrates 1/(x-yy), hence minus sign.
                val = -quad(lambda x: eval_legendre(n,x),-1,1,
                            weight='cauchy',wvar=yy,epsabs=1e-11)[0]/2
            else:
                val = quad(lambda x: eval_legendre(n,x)/(yy-x),-1,1,
                           epsabs=1e-11)[0]/2
            np.testing.assert_allclose(real_legendre_q(n,yy),val,atol=2e-12,rtol=1e-9)
    # At l=0, compare a full smeared row with direct independent PV integration.
    B=11.999
    specs=[(B,0,n) for n in range(4)]
    mus=np.array([4.01,5.,7.13,11.8,12.5,15.,100.,1e6])
    bank=build_roy_bank(mus,[0],specs,64,64,normalized=False)
    zn, wn=leggauss(90)
    z=(zn+1)/2
    for n in range(4):
        for im,m in enumerate(mus[:-1]):
            def rest(s):
                t=-(s-4)*(1-z)/2
                uu=4-s-t
                return np.sum(wn*(1/(m-t)+1/(m-uu)-3/(m-C)))
            w=lambda s:eval_legendre(n,2*(s-4)/(B-4)-1)/(B-4)
            reg=quad(lambda s:w(s)*rest(s),4,B,epsabs=2e-11)[0]
            pv=4*float(real_legendre_q(n,2*(m-4)/(B-4)-1))/(B-4)
            np.testing.assert_allclose(bank['v'][n,0,im],-(reg+pv),atol=4e-11,rtol=2e-8)
    print("roy_kernels self-test passed: divided differences, Legendre-Q PV, scalar rows")


if __name__ == '__main__':
    _self_test()
