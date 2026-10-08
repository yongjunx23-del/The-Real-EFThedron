"""Shared, explicit row catalogue for directional massive searches."""
from __future__ import annotations
import numpy as np
import scipy.special as sp
from massive_search import spectral_grid, coupling_kernels, forward_rows, first_crossing
from roy_kernels import build_roy_bank
from transfer_kernels import build_transfer_bank, build_tapered_transfer_bank

def build_problem(n=24,lmax=20,hmax=20,roy_order=10,roy_spin=4,
                  transfer_order=5,transfer_taper=1,cross_degree=4,
                  radii=(),B=11.9,threshold_probe=True,fractional=False,
                  windows=((12.,24.),(24.,48.))):
    s,ell,W,p=spectral_grid(n,lmax,B)
    if threshold_probe:
        probe=np.array([4.+1e-10])
        pp=sp.eval_legendre(ell[:,None],(probe-4/3)/(probe-4))/(probe-4/3)**3
        s=np.r_[probe,s]; W=np.column_stack([np.zeros(len(ell)),W]);p=np.column_stack([pp,p])
    u0,v0,_=coupling_kernels(s,ell)
    rows=forward_rows(s,ell,hmax,radii)
    if cross_degree==4: rows.append(('X4',np.zeros_like(p),first_crossing(s,ell)))
    else:
        from crossing_kernels import crossing_rows
        rows.extend(crossing_rows(s,ell,cross_degree))
    if roy_order>=0:
        specs=[(B,J,k) for J in range(0,roy_spin+1,2) for k in range(roy_order+1)]
        bank=build_roy_bank(s,ell,specs,normalized=False)
        rows.extend((name,u+c*u0,v+c*v0) for name,u,v,c in zip(bank['labels'],bank['u'],bank['v'],bank['g0_coeff']))
    if transfer_order>=0:
        specs=[(A,Bw,4/3,j,k) for A,Bw in windows for j in range(3) for k in range(transfer_order+1)]
        bank=(build_transfer_bank(s,ell,specs,normalized=False) if transfer_taper==0 else
              build_tapered_transfer_bank(s,ell,specs,power=transfer_taper,normalized=False))
        rows.extend((name,u+c*u0,v+c*v0) for name,u,v,c in zip(bank['labels'],bank['u'],bank['v'],bank['g0_coeff']))
    if fractional:
        from existing_solver.massive_forward_kernels import fractional_bank
        uu,vv,names=fractional_bank(s,ell)
        rows.extend(zip(names,uu,vv))
    return dict(s=s,ell=ell,W=W,p=p,u0=u0,v0=v0,rows=rows)

