"""Massive d=4 forward kernels from The Real EFThedron, Appendix A5.

Conventions: m=1, f_l=r_l+i b_l, r_l**2+(b_l-1)**2<=1.
Every returned real kernel pair (u,v) means
    sum_{l even} (2l+1) integral_4^infty [u_l(s) r_l(s)+v_l(s)b_l(s)]
                                      ds/(2*pi*sqrt(1-4/s)).
Thus complex q means u=Im(q), v=Re(q), i.e. Im(q*f).
The spin/energy array layout is always (n_spin,n_energy), with a preceding
row index for banks. No quadrature weights or spin factors are included.

These are exact kernel identities under the paper's stated analytic assumptions.
They do not certify energy quadrature or an omitted-spin tail.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.special import eval_legendre


def _inputs(s, ell):
    s = np.atleast_1d(np.asarray(s, dtype=float))
    ell = np.atleast_1d(np.asarray(ell, dtype=int))
    if s.ndim != 1 or np.any(s <= 4) or not np.all(np.isfinite(s)):
        raise ValueError("s must be a finite one-dimensional array strictly above 4")
    if ell.ndim != 1 or np.any(ell < 0) or np.any(ell % 2):
        raise ValueError("ell must contain nonnegative even physical spins")
    return s, ell


def forward_geometry(s):
    """Return (F,omega,log_D), D=((s-2)/2)^2-1 >0.

    omega=(2+i*sqrt(s*(s-4)))**2/(s-2)**2 avoids inverse trig.
    log_D=log(s)+log(s-4)-log(4) remains accurate near threshold.
    """
    s = np.asarray(s, dtype=float)
    if np.any(s <= 4):
        raise ValueError("s must exceed 4")
    q = s-2
    omega = (2/q + 1j*np.sqrt(s*(s-4))/q)**2
    return 32/q**3, omega, np.log(s)+np.log(s-4)-math.log(4)


def _pair(q, ell):
    q = np.asarray(q)
    shape = (*q.shape[:-1], len(ell), q.shape[-1])
    return (np.broadcast_to(q.imag[..., None, :], shape).copy(),
            np.broadcast_to(q.real[..., None, :], shape).copy())


def coupling_kernels(s, ell):
    """Return (u,v) with row order (g0,g2), Eq. (A5.1),(A5.3).

    g0 = <F Im(omega*f) + C_l*b>,
    C_l = 8/[(3s-8)(3s-4)] * [1/(s-2) - 4 P_l(Z)/s],
    g2 = <P_l(Z)/Q^3 * b>, Q=s-4/3, Z=Q/(s-4).
    """
    s, ell = _inputs(s, ell)
    F, omega, _ = forward_geometry(s)
    Q = s-4/3
    P = eval_legendre(ell[:, None], (Q/(s-4))[None, :])
    C = 8/((3*s-8)*(3*s-4))*(1/(s-2)-4*P/s)
    u = np.zeros((2, len(ell), len(s)))
    v = np.empty_like(u)
    u[0] = F*omega.imag
    v[0] = F*omega.real+C
    v[1] = P/Q**3
    if not np.all(np.isfinite(v)):
        raise FloatingPointError("Legendre kernels overflow: reduce spin or rescale kernels")
    return u, v


def crossing_x4(s, ell):
    """First massive crossing null, Eq. (A5.6), in the common measure.

    X4=<[2 Z^2 P_l''(Z)-3 Z P_l'(Z)]/Q^5 * b>=0.
    Derivatives use the differentiated three-term recurrence, avoiding
    division by Z^2-1 in the ultraviolet.
    """
    s, ell = _inputs(s, ell)
    Q=s-4/3
    Z=Q/(s-4)
    maxell=int(np.max(ell)) if len(ell) else 0
    ps=[np.ones_like(s),Z.copy()]
    ds=[np.zeros_like(s),np.ones_like(s)]
    dds=[np.zeros_like(s),np.zeros_like(s)]
    for n in range(1,maxell):
        ps.append(((2*n+1)*Z*ps[-1]-n*ps[-2])/(n+1))
        ds.append(((2*n+1)*(ps[-2]+Z*ds[-1])-n*ds[-2])/(n+1))
        dds.append(((2*n+1)*(2*ds[-2]+Z*dds[-1])-n*dds[-2])/(n+1))
    v=np.stack([(2*Z**2*dds[l]-3*Z*ds[l])/Q**5 for l in ell])
    return np.zeros_like(v),v


def harmonic_bank(s, ell, maximum=64, include_h2=True):
    """Eq. (A5.32) polynomial bank, optionally retaining the H2 row.

    The improved UV rows F*(omega**j-(-1)**j*omega**2), j=3..maximum,
    are null because BOTH harmonics involved vanish. Retain H2 to span
    the same full bank H2,...,Hmaximum.
    """
    s, ell = _inputs(s, ell)
    F, omega, _ = forward_geometry(s)
    rows, names = [], []
    if include_h2:
        rows.append(F*omega**2)
        names.append("H2")
    for j in range(3, maximum+1):
        rows.append(F*(omega**j-(-1)**j*omega**2))
        names.append(f"H{j}-{'+' if j%2 else ''}{(-1)**j}H2")
    q = np.stack(rows) if rows else np.zeros((0, len(s)), complex)
    return (*_pair(q, ell), names)


def rational_generator(s, radius, cosine, kind="cos"):
    """UV-subtracted rational harmonic sum in Eq. (A5.32).

    kind='cos': m(z)=(1-r*c*z)/(1-2*r*c*z+r*r*z*z)
      = sum_{k>=0} r^k T_k(c) z^k.
    kind='sin': m(z)=r*z/(1-2*r*c*z+r*r*z*z)
      = sum_{k>=1} r^k U_{k-1}(c) z^k.
    The returned complex kernel is F*omega^2*(m(omega)-m(-1)).
    Any 0<=r<1 is admissible; near-unit radii need resolved quadrature.
    """
    r, c = float(radius), float(cosine)
    if not (0 <= r < 1 and -1 <= c <= 1):
        raise ValueError("require 0<=radius<1 and -1<=cosine<=1")
    F, omega, _ = forward_geometry(s)
    den = 1-2*r*c*omega+r*r*omega**2
    denm = 1+2*r*c+r*r
    # Evaluate the DIFFERENCE as one rational expression. Its omega+1
    # factor is exact, retaining ultraviolet suppression near omega=-1.
    if kind == "cos":
        diff = r*(omega+1)*(c+r-r*(1+r*c)*omega)/(den*denm)
    elif kind == "sin":
        diff = r*(omega+1)*(1+r*r*omega)/(den*denm)
    else:
        raise ValueError("kind must be 'cos' or 'sin'")
    return F*omega**2*diff


def rational_bank(s, ell, radii=(.5,.8,.95), cosines=(-1,-.75,-.5,-.25,0,.25,.5,.75,1)):
    """Return exact rational null rows; paper defaults give 48 rows."""
    s, ell = _inputs(s, ell)
    qs, names = [], []
    for r in radii:
        for c in cosines:
            for kind in (("cos", "sin") if abs(c) < 1 else ("cos",)):
                qs.append(rational_generator(s, r, c, kind))
                names.append(f"rational_{kind}_r{r:g}_c{c:g}")
    q = np.stack(qs) if qs else np.zeros((0,len(s)), complex)
    return (*_pair(q, ell), names)


PAPER_FRACTIONAL = [(1,-7/16,1/16),(1,-15/32,1/16),(1,-7/16,1/32)] + [
    (2,a,1/8) for a in (-3/8,-1/4,-1/8,0,1/8,1/4)]


def fractional_difference(s, n, alpha, step):
    """Forward Delta_h^(n+1) C_n(alpha), Eqs. (A1.7)-(A1.8).

    Evaluate the divided-difference combination analytically as
    z^(-2n-1)*exp(alpha*L)*expm1(h*L)^(n+1),
    z=(s-2)/2, L=log(z*z-1)-i*pi.
    This exactly sums the n+2 sampled fractional kernels and avoids
    cancellation from their alternating binomial coefficients.
    """
    if int(n) != n or n < 1 or step <= 0:
        raise ValueError("require integer n>=1 and positive step")
    if not (-.5 < alpha and alpha+(n+1)*step <= n-1+1e-14):
        raise ValueError("every sampled alpha must lie in (-1/2,n-1]")
    s = np.asarray(s, float)
    _, _, log_D = forward_geometry(s)
    L = log_D-1j*np.pi
    return np.exp(alpha*L-(2*n+1)*np.log((s-2)/2))*np.expm1(step*L)**(n+1)


def fractional_bank(s, ell, triples=PAPER_FRACTIONAL):
    s, ell = _inputs(s, ell)
    qs = [fractional_difference(s,n,a,h) for n,a,h in triples]
    names = [f"Delta{n+1}_C{n}_a{a:g}_h{h:g}" for n,a,h in triples]
    q = np.stack(qs) if qs else np.zeros((0,len(s)), complex)
    return (*_pair(q, ell), names)


def fractional_derivative(s, n, alpha, order=None):
    """Interior alpha derivatives of order >=n+1 annihilate C_n's EFT side.

    Strict upper admissibility alpha<n-1 avoids adding logarithms at the
    marginal UV endpoint. This is a new optional row, not in the 9-row bank.
    """
    order = n+1 if order is None else order
    if int(n) != n or n < 1 or int(order) != order or order < n+1:
        raise ValueError("require n>=1 and derivative order>=n+1")
    if not (-.5 < alpha < n-1):
        raise ValueError("derivative alpha must be strictly inside (-1/2,n-1)")
    s = np.asarray(s, float)
    _, _, log_D = forward_geometry(s)
    L = log_D-1j*np.pi
    return np.exp(alpha*L-(2*n+1)*np.log((s-2)/2))*L**order


def stable_disk(U,V):
    """Return support, extremal real part, absorption without V<0 cancellation."""
    U,V = np.broadcast_arrays(np.asarray(U,float), np.asarray(V,float))
    R = np.hypot(U,V)
    r = np.divide(U,R,out=np.zeros_like(R),where=R>0)
    negative = V < 0
    h = R+V
    np.divide(U*U,R-V,out=h,where=negative)
    b = np.divide(h,R,out=np.zeros_like(R),where=R>0)
    return h,r,b

