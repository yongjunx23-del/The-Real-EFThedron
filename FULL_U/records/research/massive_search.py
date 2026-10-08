"""Reconstruct massive 4D disk-moment searches from the supplied PDF.

Finite quadrature/spin SOCP results are exploration results, NOT continuum bounds.
The saved dual coefficients are candidates for independent continuum validation.
"""
from __future__ import annotations
import argparse
import json
import os
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
import scipy.special as sp
from scipy.optimize import minimize, linprog
import cvxpy as cp

ROOT = Path(__file__).resolve().parent.parent


def spectral_grid(n=24, lmax=12, B=11.9, smax=None):
    """Gauss quadrature in velocity removes the ds/rho threshold singularity."""
    stop=1. if smax is None else np.sqrt(1-4/smax)
    edges = [e for e in [0., .12, .35, np.sqrt(1-4/B), *[np.sqrt(1-4/k) for k in (12,24,48)], .95, .995] if e<stop]+[stop]
    edges = sorted(set(edges))
    x, w = sp.roots_legendre(n)
    beta = np.concatenate([(a+b)/2+(b-a)/2*x for a,b in zip(edges[:-1],edges[1:])])
    wb = np.concatenate([w*(b-a)/2 for a,b in zip(edges[:-1],edges[1:])])
    s = 4/(1-beta**2)
    ell = np.arange(0, lmax+1, 2)
    W = (2*ell[:,None]+1)*4/np.pi*wb/(1-beta**2)**2
    Q = s-4/3
    p = sp.eval_legendre(ell[:,None], Q/(s-4))/Q**3
    return s, ell, W, p


def coupling_kernels(s, ell):
    s = np.asarray(s)[None,:]
    ell = np.asarray(ell)[:,None]
    Q = s-4/3
    P = sp.eval_legendre(ell, Q/(s-4))
    F = 32/(s-2)**3
    omega = ((2+1j*np.sqrt(s*(s-4)))/(s-2))**2
    # Eq. A5.3: forward subtraction matching to the symmetric point.
    C = 8/((3*s-8)*(3*s-4))*(1/(s-2)-4*P/s)
    u0 = np.broadcast_to(F*omega.imag,P.shape).copy()
    v0 = F*omega.real+C
    p = P/Q**3
    return u0, v0, p


def first_crossing(s, ell):
    s=np.asarray(s)[None,:]
    ell=np.asarray(ell)[:,None]
    Q=s-4/3; d=s-4; Z=Q/d
    # Jacobi derivative formula: d^k P_l/dz^k=2^-k(l+1)_k P_(l-k)^(k,k).
    P1=np.zeros(np.broadcast_shapes(ell.shape,s.shape))
    P2=P1.copy()
    for i,l in enumerate(ell[:,0]):
        if l>=1: P1[i]=(l+1)/2*sp.eval_jacobi(int(l-1),1,1,Z[0])
        if l>=2: P2[i]=(l+1)*(l+2)/4*sp.eval_jacobi(int(l-2),2,2,Z[0])
    return 2*P2/(Q**3*d**2)-3*P1/(Q**4*d)


def forward_rows(s, ell, hmax=12, radii=()):
    s=np.asarray(s)[None,:]
    F=32/(s-2)**3
    omega=((2+1j*np.sqrt(s*(s-4)))/(s-2))**2
    rows=[]
    shape=(len(ell),s.size)
    for j in range(2,hmax+1):
        q=F*omega**j
        rows.append((f"H{j}",np.broadcast_to(q.imag,shape),np.broadcast_to(q.real,shape)))
    for r in radii:
        for c in [-1.,-.75,-.5,0.,.5,.75,1.]:
            den=1-2*r*c*omega+r*r*omega**2
            m=(1-r*c*omega)/den
            atminus=(1+r*c)/(1+2*r*c+r*r)
            q=F*omega**2*(m-atminus)
            rows.append((f"Rc_{r}_{c}",np.broadcast_to(q.imag,shape),np.broadcast_to(q.real,shape)))
    return rows


def solve_grid(s,ell,W,p,u0,v0,rows,direction=(1.,0.),verbose=False):
    """Cone scaling B=w*p*b, R=sqrt(w*p)*r retains exact grid disks."""
    shape=p.shape
    q=(W*p).ravel()
    R=cp.Variable(q.size)
    B=cp.Variable(q.size)
    def coefficients(u,v):
        return (np.sqrt(W/p)*u).ravel(), (v/p).ravel()
    a0,b0=coefficients(u0,v0)
    ar=[];br=[];scales=[];names=[]
    for name,u,v in rows:
        a,b=coefficients(u,v)
        scale=max(np.linalg.norm(a),np.linalg.norm(b),1e-100)
        ar.append(a/scale);br.append(b/scale);scales.append(scale);names.append(name)
    A=np.array(ar); C=np.array(br)
    # Equilibrate dependent/poorly-conditioned rows by an SVD of both components.
    if len(rows):
        ua,sv,_=np.linalg.svd(np.column_stack([A,C]),full_matrices=False)
        keep=sv>sv[0]*1e-10
        whitening=(ua[:,keep]/sv[keep]).T
        Awhite=whitening@A; Cwhite=whitening@C
    else:
        whitening=np.zeros((0,0)); Awhite=A;Cwhite=C;keep=[];sv=[]
    tref=np.maximum(np.minimum(np.sqrt(q),1.),1e-8)
    conic=cp.SOC(B+tref,cp.vstack([cp.multiply(np.sqrt(2*tref),R),cp.multiply(np.sqrt(2*tref/q),B),B-tref]),axis=0)
    constraints=[conic]
    if len(rows): constraints.append(Awhite@R+Cwhite@B==0)
    g0=a0@R+b0@B
    g2=cp.sum(B)
    prob=cp.Problem(cp.Maximize(direction[0]*g0+direction[1]*g2),constraints)
    start=time.time()
    prob.solve(solver="CLARABEL",max_iter=250,tol_gap_abs=1e-8,tol_gap_rel=1e-8,
               tol_feas=1e-9,verbose=verbose)
    if R.value is None:
        return {"status":prob.status,"support":None,"seconds":time.time()-start}, None
    if len(rows):
        lam=-(whitening.T@constraints[1].dual_value)/np.array(scales)
    else: lam=np.array([])
    U=direction[0]*u0+sum((la*r[1] for la,r in zip(lam,rows)),np.zeros(shape))
    V=direction[0]*v0+direction[1]*p+sum((la*r[2] for la,r in zip(lam,rows)),np.zeros(shape))
    hypot=np.hypot(U,V)
    h=np.where(V<0,U*(U/np.maximum(hypot-V,1e-300)),hypot+V)
    dual=float(np.sum(W*h))
    rv=R.value/np.sqrt(q);bv=B.value/q
    result={"status":prob.status,"support":float(prob.value),"g0":float(g0.value),
            "g2":float(g2.value),"dual_support_grid":dual,"grid_dual_gap":dual-float(prob.value),
            "max_disk_violation":float(np.max(rv**2+(bv-1)**2-1)),
            "max_scaled_null_residual":float(np.max(np.abs(Awhite@R.value+Cwhite@B.value))) if len(rows) else 0.,
            "row_count":len(rows),"row_rank":int(sum(keep)),"points":q.size,
            "seconds":time.time()-start,"claim":"finite-grid, finite-spin diagnostic only"}
    state={"multipliers":lam,"names":np.array(names),"U":U,"V":V,"s":s,"ell":ell,
           "W":W,"p":p,"r":rv.reshape(shape),"b":bv.reshape(shape),"direction":direction}
    return result,state


def solve_dual(s,ell,W,p,u0,v0,rows,direction,warm=None):
    shape=p.shape
    Ubase=direction[0]*u0
    Vbase=direction[0]*v0+direction[1]*p
    ur=np.array([r[1] for r in rows]); vr=np.array([r[2] for r in rows])
    def fg(lam,hessian=False):
        U=Ubase+np.einsum('i,ijk->jk',lam,ur)
        V=Vbase+np.einsum('i,ijk->jk',lam,vr)
        rad=np.maximum(np.hypot(U,V),1e-200)
        r=U/rad
        b=np.where(V<0,r*(U/np.maximum(rad-V,1e-200)),1+V/rad)
        h=np.where(V<0,U*(U/np.maximum(rad-V,1e-200)),rad+V)
        val=float(np.sum(W*h))
        grad=np.einsum('jk,ijk->i',W*r,ur)+np.einsum('jk,ijk->i',W*b,vr)
        if hessian:
            vec=(ur*(V/rad)[None,:,:]-vr*r[None,:,:])*np.sqrt(W/rad)[None,:,:]
            return val,grad,vec.reshape(len(rows),-1)@vec.reshape(len(rows),-1).T
        return val,grad
    lam=np.zeros(len(rows))
    if warm:
        state=np.load(warm)
        lookup=dict(zip(state['names'],state['multipliers']))
        lam=np.array([lookup.get(r[0],0.) for r in rows])
    if fg(lam)[0]>1e8:
        # Find a threshold-feasible starting witness. This only seeds the search.
        mask=p>1e8
        mat=(vr/p[None,:,:])[:,mask].T
        rhs=-(Vbase/p)[mask]-1e-5
        scales=np.maximum(np.linalg.norm(mat,axis=0),1e-8)
        lp=linprog(np.r_[np.zeros(len(rows)),np.ones(len(rows))],
            A_ub=np.vstack([np.column_stack([mat/scales,np.zeros_like(mat)]),
                            np.column_stack([np.eye(len(rows)),-np.eye(len(rows))]),
                            np.column_stack([-np.eye(len(rows)),-np.eye(len(rows))])]),
            b_ub=np.r_[rhs,np.zeros(2*len(rows))],
            bounds=[(None,None)]*len(rows)+[(0,None)]*len(rows),method='highs')
        if lp.success: lam=lp.x[:len(rows)]/scales
    start=time.time()
    history=[]
    for stage in range(5):
        val,grad,H=fg(lam,True)
        eig,vec=np.linalg.eigh(H)
        floor=max(np.max(eig)*1e-9,1e-12)
        trans=vec/np.sqrt(np.maximum(eig,floor))
        base=lam.copy()
        def fun(x):
            v,g=fg(base+trans@x)
            return v,trans.T@g
        opt=minimize(fun,np.zeros(len(rows)),jac=True,method='L-BFGS-B',
                     options={'maxiter':1200,'maxls':60,'ftol':1e-13,'gtol':1e-8,'maxcor':40})
        lam=base+trans@opt.x
        value,g=fg(lam)
        history.append({'value':value,'nit':int(opt.nit),'message':str(opt.message),
                        'gradient_norm':float(np.linalg.norm(g))})
        print(json.dumps({'stage':stage,**history[-1]}),flush=True)
        if abs(val-value)<1e-10 and np.linalg.norm(trans.T@g)<1e-6: break
    U=Ubase+np.einsum('i,ijk->jk',lam,ur)
    V=Vbase+np.einsum('i,ijk->jk',lam,vr)
    rad=np.maximum(np.hypot(U,V),1e-200)
    r=U/rad; b=np.where(V<0,r*U/np.maximum(rad-V,1e-200),1+V/rad)
    result={'status':'direct_dual_search','support':value,'dual_support_grid':value,
            'g0':float(np.sum(W*(u0*r+v0*b))),'g2':float(np.sum(W*p*b)),
            'row_count':len(rows),'points':p.size,'seconds':time.time()-start,
            'history':history,'claim':'finite-grid, finite-spin diagnostic only'}
    state={'multipliers':lam,'names':np.array([r[0] for r in rows]),'U':U,'V':V,
           's':s,'ell':ell,'W':W,'p':p,'r':r,'b':b,'direction':direction}
    return result,state


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=16)
    ap.add_argument("--lmax",type=int,default=8)
    ap.add_argument("--hmax",type=int,default=8)
    ap.add_argument("--B",type=float,default=11.9)
    ap.add_argument("--roy-order",type=int,default=-1)
    ap.add_argument("--roy-spin",type=int,default=4)
    ap.add_argument("--direction",choices=["maxg0","ming0","maxg2"],default="maxg0")
    ap.add_argument("--radii",type=float,nargs="*",default=[])
    ap.add_argument("--out",default="pilot")
    ap.add_argument("--verbose",action="store_true")
    ap.add_argument("--engine",choices=['cone','dual','newton'],default='newton')
    ap.add_argument("--warm")
    ap.add_argument("--smax",type=float)
    ap.add_argument("--transfer-order",type=int,default=-1)
    ap.add_argument("--transfer-taper",type=int,default=0)
    ap.add_argument("--cross-degree",type=int,default=4)
    ap.add_argument("--threshold-margin",type=float,default=1e-6)
    args=ap.parse_args()
    s,ell,W,p=spectral_grid(args.n,args.lmax,args.B,args.smax)
    if args.engine=='newton':
        # A zero-weight feasibility node constrains the leading threshold sign.
        # Its separation is safely above binary64 spacing at 4. The positive
        # margin also protects against extrapolation/cancellation at s=4.
        probe=np.array([4.+1e-10])
        pp=sp.eval_legendre(ell[:,None],(probe-4/3)/(probe-4))/(probe-4/3)**3
        s=np.r_[probe,s]; W=np.column_stack([np.zeros(len(ell)),W]); p=np.column_stack([pp,p])
    u0,v0,_=coupling_kernels(s,ell)
    rows=forward_rows(s,ell,args.hmax,args.radii)
    if args.cross_degree==4:
        rows.append(("X4",np.zeros_like(p),first_crossing(s,ell)))
    else:
        from crossing_kernels import crossing_rows
        rows.extend(crossing_rows(s,ell,args.cross_degree))
    if args.roy_order>=0:
        from roy_kernels import build_roy_bank
        specs=[(args.B,J,n) for J in range(0,args.roy_spin+1,2) for n in range(args.roy_order+1)]
        bank=build_roy_bank(s,ell,specs,normalized=False)
        rows.extend((name,u+c*u0,v+c*v0) for name,u,v,c in zip(bank['labels'],bank['u'],bank['v'],bank['g0_coeff']))
    if args.transfer_order>=0:
        from transfer_kernels import build_transfer_bank,build_tapered_transfer_bank
        specs=[(A,B,4/3,j,n) for A,B in [(12,24),(24,48)] for j in range(3) for n in range(args.transfer_order+1)]
        bank=(build_transfer_bank(s,ell,specs,normalized=False) if args.transfer_taper==0
              else build_tapered_transfer_bank(s,ell,specs,power=args.transfer_taper,normalized=False))
        rows.extend((name,u+c*u0,v+c*v0) for name,u,v,c in zip(bank['labels'],bank['u'],bank['v'],bank['g0_coeff']))
    direction={"maxg0":(1.,0.),"ming0":(-1.,0.),"maxg2":(0.,1.)}[args.direction]
    if args.engine=='newton':
        from dual_newton import solve_dual as newton_solve
        result,state=newton_solve(s,ell,W,p,u0,v0,rows,direction,args.warm,eps=args.threshold_margin)
    else:
        result,state=(solve_grid(s,ell,W,p,u0,v0,rows,direction,args.verbose) if args.engine=='cone'
                      else solve_dual(s,ell,W,p,u0,v0,rows,direction,args.warm))
    result["config"]=vars(args)
    out=ROOT/"output"/args.out
    out.parent.mkdir(exist_ok=True,parents=True)
    out.with_suffix(".json").write_text(json.dumps(result,indent=2))
    if state is not None: np.savez_compressed(out.with_suffix(".npz"),**state)
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__": main()
