"""Constrained, damped Newton search for finite-grid disk support witnesses.

Threshold sign inequalities restrict dual search only. They do not impose an
extra assumption on physical amplitudes. This module DOES NOT certify spectral
quadrature, omitted spins, or continuum constraints.
"""
from __future__ import annotations
import json
import time
import warnings
from pathlib import Path
import numpy as np
import cvxpy as cp


class DiskSupport:
    def __init__(self, Ubase,Vbase,ur,vr,W,p):
        self.shape=np.shape(p)
        self.p=np.asarray(p,dtype=float).ravel()
        self.w=np.broadcast_to(W,self.shape).astype(float).ravel()
        self.ub=np.asarray(Ubase,dtype=float).ravel()
        self.vb=np.asarray(Vbase,dtype=float).ravel()/self.p
        self.ur=np.asarray(ur,dtype=float).reshape(len(ur),-1)
        self.vr=np.asarray(vr,dtype=float).reshape(len(vr),-1)
        self.vnorm=self.vr/self.p[None,:]
        if not all(np.isfinite(a).all() for a in [self.p,self.w,self.ub,self.vb,self.ur,self.vr]):
            raise ValueError('Input contains nonfinite coefficients')
        if np.any(self.p<=0) or np.any(self.w<0):
            raise ValueError('p must be positive and W nonnegative')

    def fields(self,lam):
        # Computing V/p first avoids cancellation of astronomically large p.
        return self.ub+lam@self.ur, self.p*(self.vb+lam@self.vnorm)

    def evaluate(self,lam,hessian=False):
        U,V=self.fields(lam)
        rad=np.maximum(np.hypot(U,V),1e-300)
        r=U/rad
        neg=V<0
        b=np.empty_like(V);cost=np.empty_like(V)
        b[neg]=r[neg]*(U[neg]/(rad[neg]-V[neg]))
        cost[neg]=U[neg]*(U[neg]/(rad[neg]-V[neg]))
        b[~neg]=1+V[~neg]/rad[~neg]
        cost[~neg]=rad[~neg]+V[~neg]
        value=float(self.w@cost)
        grad=self.ur@(self.w*r)+self.vr@(self.w*b)
        if not hessian:
            return value,grad
        tangent=(self.ur*(V/rad)[None,:]-self.vr*r[None,:])*np.sqrt(self.w/rad)[None,:]
        H=tangent@tangent.T
        return value,grad,H


def _solve_qp(linear,A,rhs,trust_radius=None,verbose=False):
    z=cp.Variable(len(linear))
    constraints=[]
    if len(rhs):
        # Row scaling is essential: finite-spin threshold coefficients vary
        # wildly even after division by the positive coupling kernel p.
        row_scale=np.maximum(np.linalg.norm(A,axis=1),1e-30)
        constraints.append((A/row_scale[:,None])@z<=rhs/row_scale)
    if trust_radius is not None:
        constraints.append(cp.norm(z,2)<=trust_radius)
    problem=cp.Problem(cp.Minimize(.5*cp.sum_squares(z)+linear@z),constraints)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',UserWarning)
        try:
            problem.solve(solver='CLARABEL',max_iter=250,tol_gap_abs=1e-10,
                          tol_gap_rel=1e-10,tol_feas=1e-11,verbose=verbose)
        except cp.error.SolverError as exc:
            return None, {'status':'solver_error','error':str(exc)}
    return (None if z.value is None else np.asarray(z.value)), {
        'status':problem.status,'iterations':problem.solver_stats.num_iters,
        'predicted_model_change':None if problem.value is None else float(problem.value)}


def solve_constrained_dual(Ubase,Vbase,ur,vr,W,p,start=None,*,threshold=1e4,
                           eps=1e-10,max_iter=100,trust_radius=10.,
                           eigen_floor=1e-8,ftol=1e-11,gtol=1e-7,
                           verbose=False,log_path=None):
    """Return (multipliers, diagnostics), retaining every accepted decrease.

    ur and vr have shape (rows, *p.shape). Constraints apply wherever p>threshold.
    eps is a strictly negative margin for V/p. trust_radius bounds the norm of
    each Hessian-whitened step, not the final multipliers; None disables it.
    A failed Newton or line-search step cannot overwrite a better accepted point.
    """
    tic=time.time(); fun=DiskSupport(Ubase,Vbase,ur,vr,W,p)
    nr=len(ur); lam=np.zeros(nr) if start is None else np.asarray(start,dtype=float).copy()
    if lam.shape!=(nr,): raise ValueError('Wrong number of initial multipliers')
    mask=fun.p>threshold
    C=fun.vnorm[:,mask].T; target=-eps-fun.vb[mask]
    history=[]; value_before_projection=fun.evaluate(lam)[0]
    initial_violation=float(np.max(C@lam-target)) if mask.any() else 0.
    projection=None
    if initial_violation>0:
        # Projection in the coefficient row scales used by the primal cone
        # formulation, not unscaled lambda distance.
        scale=np.maximum(np.linalg.norm(fun.ur*np.sqrt(fun.w/fun.p)[None,:],axis=1),
                         np.linalg.norm(fun.vnorm,axis=1))
        scale=np.maximum(scale,1e-20)
        # Strict additional margin protects against projection solver tolerance.
        z,projection=_solve_qp(np.zeros(nr),C/scale[None,:],
                               target-C@lam-max(eps*.01,1e-12))
        if z is None:
            return lam,{'status':'initial_projection_failed','projection':projection,
                        'support':value_before_projection,'seconds':time.time()-tic}
        candidate=lam+z/scale
        violation=float(np.max(C@candidate-target))
        projection['threshold_violation']=violation
        if violation>max(eps*.1,1e-10):
            return lam,{'status':'initial_projection_inaccurate','projection':projection,
                        'support':value_before_projection,'seconds':time.time()-tic}
        lam=candidate
    value,grad,H=fun.evaluate(lam,True)
    initial_value=value
    if not np.isfinite(value):
        raise FloatingPointError('Nonfinite support at feasible initialization')
    status='iteration_limit'
    radius=trust_radius
    for iteration in range(max_iter):
        # First diagonal-equilibrate H, then whiten its correlated row space.
        diag=np.sqrt(np.maximum(np.diag(H),1e-40))
        Hn=H/diag[:,None]/diag[None,:]
        eig,vec=np.linalg.eigh((Hn+Hn.T)*.5)
        floor=max(float(eig.max())*eigen_floor,1e-12)
        trans=(vec/np.sqrt(np.maximum(eig,floor))[None,:])/diag[:,None]
        linear=trans.T@grad
        slack=target-C@lam
        z,qp=_solve_qp(linear,C@trans,slack,radius)
        if z is None:
            status='newton_qp_failed';history.append({'iteration':iteration,'qp':qp});break
        step=trans@z
        descent=float(grad@step)
        # Upper bound on feasible alpha, retaining an interior safety factor
        # only for constraints approached by this particular direction.
        ds=C@step
        towards=ds>0
        alpha=1.
        if towards.any():
            alpha=min(alpha,float(np.min(np.maximum(slack[towards],0)/ds[towards]))*.995)
        # QP near-feasibility can generate a small outward direction at an
        # active row. Treat it as convergence/failure, never cross the boundary.
        if descent>=0 or alpha<=1e-15:
            status='stationary_or_active_constraint'
            history.append({'iteration':iteration,'value':value,'descent':descent,
                            'alpha':alpha,'qp':qp});break
        accepted=False
        for backtracks in range(70):
            trial=lam+alpha*step
            violation=float(np.max(C@trial-target)) if mask.any() else 0.
            trial_value,trial_grad=fun.evaluate(trial)
            if np.isfinite(trial_value) and violation<=max(eps*.01,1e-13) and trial_value<=value+1e-4*alpha*descent:
                accepted=True;break
            alpha*=.5
        if not accepted:
            status='line_search_failed'
            history.append({'iteration':iteration,'value':value,'descent':descent,'alpha':alpha,'qp':qp});break
        old_value=value;lam=trial
        value,grad,H=fun.evaluate(lam,True)
        entry={'iteration':iteration,'value':value,'decrease':old_value-value,
               'alpha':alpha,'descent':descent,'newton_decrement':float(np.sqrt(max(-descent,0))),
               'threshold_violation':violation,'qp':qp}
        history.append(entry)
        if verbose: print(json.dumps(entry),flush=True)
        if log_path:
            Path(log_path).write_text(json.dumps(history,indent=2),encoding='utf8')
        if old_value-value<ftol*max(1.,abs(value)) and -descent<gtol**2:
            status='converged';break
        if radius is not None:
            if alpha>.9: radius=min(radius*1.4,100.)
            elif alpha<.1: radius=max(radius*.5,.05)
    U,V=fun.fields(lam)
    final_violation=float(np.max(C@lam-target)) if mask.any() else 0.
    diagnostics={'status':status,'support':float(value),'dual_support_grid':float(value),
        'initial_support':float(initial_value),'support_before_projection':float(value_before_projection),
        'accepted_improvement':float(initial_value-value),'threshold':threshold,'eps':eps,
        'threshold_point_count':int(mask.sum()),'initial_threshold_violation':initial_violation,
        'final_threshold_violation':final_violation,'threshold_positive_cost_points':int(np.sum(V[mask]>0)),
        'row_count':nr,'points':fun.p.size,'projection':projection,'history':history,
        'seconds':time.time()-tic,'claim':'finite-grid, finite-spin dual witness only'}
    return lam,diagnostics


def solve_dual(s,ell,W,p,u0,v0,rows,direction,warm=None,**kwargs):
    """Optional replacement for massive_search.solve_dual, same basic API."""
    lam=np.zeros(len(rows))
    if warm:
        old=np.load(warm);lookup=dict(zip(old['names'],old['multipliers']))
        lam=np.array([lookup.get(row[0],0.) for row in rows])
    ur=np.array([row[1] for row in rows]);vr=np.array([row[2] for row in rows])
    ub=direction[0]*u0;vb=direction[0]*v0+direction[1]*p
    lam,result=solve_constrained_dual(ub,vb,ur,vr,W,p,lam,**kwargs)
    fun=DiskSupport(ub,vb,ur,vr,W,p);U,V=fun.fields(lam)
    U=U.reshape(p.shape);V=V.reshape(p.shape)
    rad=np.maximum(np.hypot(U,V),1e-300)
    r=U/rad;b=np.empty_like(r);neg=V<0
    b[neg]=r[neg]*U[neg]/(rad[neg]-V[neg]);b[~neg]=1+V[~neg]/rad[~neg]
    result['g0']=float(np.sum(W*(u0*r+v0*b)))
    result['g2']=float(np.sum(W*p*b))
    state={'multipliers':lam,'names':np.array([row[0] for row in rows]),
           'U':U,'V':V,'s':s,'ell':ell,'W':W,'p':p,'r':r,'b':b,'direction':direction}
    return result,state


def _self_test():
    """Independent derivatives, conic optimum, and infeasible-seed projection."""
    rng=np.random.default_rng(18)
    ub=rng.normal(size=(2,5));vb=-2+rng.normal(size=(2,5))
    ur=rng.normal(size=(3,2,5));vr=rng.normal(size=(3,2,5))
    w=np.ones((2,5));p=np.ones((2,5))
    fun=DiskSupport(ub,vb,ur,vr,w,p)
    x=np.array([.1,-.2,.3]);value,grad,H=fun.evaluate(x,True)
    h=1e-5;basis=np.eye(3)
    fdg=np.array([(fun.evaluate(x+h*e)[0]-fun.evaluate(x-h*e)[0])/(2*h) for e in basis])
    fdH=np.column_stack([(fun.evaluate(x+h*e)[1]-fun.evaluate(x-h*e)[1])/(2*h) for e in basis])
    np.testing.assert_allclose(grad,fdg,atol=1e-8)
    np.testing.assert_allclose(H,fdH,atol=1e-8)
    lam,diag=solve_constrained_dual(ub,vb,ur,vr,w,p,x,max_iter=80)
    z=cp.Variable(3)
    uu=ub.ravel()+ur.reshape(3,-1).T@z
    vv=vb.ravel()+vr.reshape(3,-1).T@z
    cone=cp.Problem(cp.Minimize(cp.sum(cp.norm(cp.vstack([uu,vv]),2,axis=0)+vv)))
    cone.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-11)
    np.testing.assert_allclose(diag['support'],cone.value,atol=1e-8)
    la,projected=solve_constrained_dual(np.array([1.]),np.array([-1e5]),
        np.array([[1.]]),np.array([[1e5]]),np.array([1.]),np.array([1e5]),
        np.array([2.]),max_iter=100)
    assert projected['support']<1e-9
    assert projected['threshold_positive_cost_points']==0
    assert projected['projection'] is not None
    print(json.dumps({'gradient_max_error':float(max(abs(grad-fdg))),
        'hessian_max_error':float(max(abs(H-fdH).ravel())),
        'independent_conic_support':float(cone.value),'newton_support':diag['support'],
        'projected_analytic_zero_support':projected['support']}))


if __name__=='__main__':
    _self_test()
