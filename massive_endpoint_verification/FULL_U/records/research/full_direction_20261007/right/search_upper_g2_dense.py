"""Direct unit-normal searches using the saved 220-row lower-sector bank.

Numerical objective only. Every run records load/build/projection/optimization
time; source and crossing components must be checked on the saved witness.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import sys,time,json,argparse,hashlib,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
for p in (ROOT/'research',ROOT/'research/directional_certification',ROOT/'research/lower_g0_upper_g2_20261006'):
    sys.path.insert(0,str(p))
import numpy as np
from angle_audit import normalized_inventory
from numerical_dual_newton import DiskSupport,solve_constrained_dual
import numerical_dual_newton as newton_module
from numerical_direction_search import harmonic_tail_transform
from lower_search import feasibility_project
from matrix_crossing_series import optimizer_guard_matrix
from scipy.optimize import linprog

BASE=ROOT/'output/lower_g0_upper_g2_20261006'

def cut_project(ub,vb,U,V,W,p,lam,eps=1e-8):
    fun=DiskSupport(ub,vb,U,V,W,p);mask=p>1e4
    C=fun.vnorm[:,mask].T;target=-eps-fun.vb[mask]
    scale=np.maximum(np.max(abs(C),axis=1),1.)
    C=C/scale[:,None];target=target/scale-1e-8
    column=np.maximum(np.max(abs(C),axis=0),1e-12)
    A=C/column[None,:];rhs=target-C@lam;nr=len(lam)
    active=set(np.argsort(rhs)[:1000].tolist());started=time.perf_counter()
    for iteration in range(20):
        idx=np.array(sorted(active));AA=A[idx]
        sol=linprog(np.r_[np.zeros(nr),np.ones(nr)],
          A_ub=np.vstack([np.c_[AA,np.zeros_like(AA)],np.c_[np.eye(nr),-np.eye(nr)],np.c_[-np.eye(nr),-np.eye(nr)]]),
          b_ub=np.r_[rhs[idx],np.zeros(2*nr)],bounds=[(None,None)]*nr+[(0,None)]*nr,
          method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9,'time_limit':45})
        if not sol.success:raise RuntimeError('Cut projection failed: '+sol.message)
        candidate=lam+sol.x[:nr]/column;residual=C@candidate-target
        print('CUT_PROJECT',iteration,len(active),float(residual.max()),'seconds',time.perf_counter()-started,flush=True)
        if residual.max()<=2e-9:return candidate
        bad=np.flatnonzero(residual>2e-9);add=bad[np.argsort(residual[bad])[-1000:]]
        old=len(active);active.update(add.tolist())
        if old==len(active):break
    raise RuntimeError('Cut projection failed full-grid separation')

def enable_constraint_generation():
    original=newton_module._solve_qp
    previous_active=set()
    def generated(linear,A,rhs,trust_radius=None,verbose=False):
        nonlocal previous_active
        if len(rhs)<1800:return original(linear,A,rhs,trust_radius,verbose)
        scale=np.maximum(np.linalg.norm(A,axis=1),1e-30)
        # A strict inward cushion prevents solver rounding at an active row
        # from forcing the full-grid feasible line-search step to zero.
        inner_rhs=rhs-1e-9*scale
        normalized_rhs=inner_rhs/scale
        active=set(np.argsort(normalized_rhs)[:800].tolist()) | previous_active
        rounds=[];started=time.perf_counter()
        for iteration in range(12):
            subset=np.array(sorted(active))
            z,diag=original(linear,A[subset],inner_rhs[subset],trust_radius,verbose)
            if z is None:return z,diag
            residual=(A@z-inner_rhs)/scale
            maximum=float(residual.max())
            rounds.append({'rows':len(active),'max_normalized_violation':maximum,'qp':dict(diag)})
            print('CUT_QP',iteration,len(active),maximum,'seconds',time.perf_counter()-started,flush=True)
            if maximum<=5e-10:
                # Constraint row identities are stable across Newton steps.
                # Retain all rows of the successful reduced QP for the next step.
                previous_active=active.copy()
                diag['constraint_generation']=rounds
                return z,diag
            bad=np.flatnonzero(residual>5e-10)
            add=bad[np.argsort(residual[bad])[-800:]]
            previous=len(active);active.update(add.tolist())
            if len(active)==previous:break
        z,diag=original(linear,A,rhs,trust_radius,verbose)
        diag['constraint_generation']=rounds
        return z,diag
    newton_module._solve_qp=generated

def grids():
    out=[(np.array([4+1e-6,4+.001,4+.01,4.1,4.5,5.,6.,8.,12.1,20.,30.,50.,75.,100.,150.,250.,500.,1000.,1e4,1e6]),np.arange(26,42,2),0.),
         (np.r_[np.arange(13.,105.,2.),150.,250.,500.,1000.,1e4,1e6],np.arange(42,74,2),1e-4),
         (np.array([8.,10.,12.1,16.,20.,30.,50.,75.,100.,150.,250.,500.,1000.,1e4,1e6]),np.array([80,100,140,200]),1e-4)]
    tau=np.unique(np.r_[np.geomspace(.001,1,12),np.arange(1,7.001,.125),np.arange(7,40.001,.25),np.geomspace(40,150,14)])
    out += [(4+ell*(ell+1)/tau,np.array([ell]),1e-4) for ell in np.r_[np.arange(42,102,2),[120,160,200,300,500,1000]]]
    return out

def general_guard_base(names,guard,gb,a,b,tail_margin):
    result=gb.copy();offset=862
    cross=[i for i,n in enumerate(names) if n.startswith('X')]
    crossnames=[names[i] for i in cross]
    worst=0.
    for s,l,oldmargin in grids():
        _,v0,_,_,v=normalized_inventory(crossnames,s,l)
        count=v0.size;old=guard[cross,offset:offset+count]
        raw=v.reshape(len(cross),-1)
        index=np.argmax(abs(old),axis=0);cols=np.arange(count)
        invscale=old[index,cols]/raw[index,cols]
        assert np.isfinite(invscale).all() and np.all(invscale>0)
        worst=max(worst,float(np.max(abs(raw*invscale-old))))
        # Recover the exact old builder's scaling from a nonzero crossing row.
        assert np.max(abs(raw*invscale-old))<1e-11
        np.testing.assert_allclose(gb[offset:offset+count],(-v0.ravel()+oldmargin)*invscale,rtol=2e-10,atol=2e-12)
        margin=tail_margin if oldmargin else 0.
        result[offset:offset+count]=(a*v0.ravel()+b+margin)*invscale
        offset+=count
    assert offset==len(gb)
    return result,worst

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--angle',type=float,default=180)
    parser.add_argument('--cache',type=Path,default=BASE/'lower_enriched_cache.npz')
    parser.add_argument('--warm',type=Path,default=BASE/'lower_bernstein_partial_frozen.npz')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--max-iter',type=int,default=8)
    parser.add_argument('--tail-margin',type=float,default=.002)
    parser.add_argument('--source-margin',type=float,default=.0001)
    parser.add_argument('--eigen-floor',type=float,default=1e-7)
    parser.add_argument('--project-only',action='store_true')
    parser.add_argument('--constraint-generation',action='store_true')
    args=parser.parse_args();started=time.perf_counter();args.out.parent.mkdir(parents=True,exist_ok=True)
    if args.constraint_generation:enable_constraint_generation()
    cache=args.cache;data=np.load(cache)
    names=list(map(str,data['names']));nr=len(names)
    ub,vb,U,V,W,p,guard,gb=[data[k].copy() for k in ('ub','vb','U','V','W','p','guard','gb')]
    cache_seconds=time.perf_counter()-started
    theta=math.radians(args.angle)
    a,b=math.cos(theta),math.sin(theta)
    if args.angle%360==180:a,b=-1.,0.
    elif args.angle%360==90:a,b=0.,1.
    elif args.angle%360==270:a,b=0.,-1.
    elif args.angle%360==0:a,b=1.,0.
    old=np.load(args.warm);lookup=dict(zip(map(str,old['names']),old['multipliers']))
    lam=np.array([lookup.get(n,0.) for n in names])
    if 'direction' in old:
        source=np.asarray(old['direction'],float)
        # Scaling the old left witness is a useful initial approximation.
        if abs(source[1])<1e-12 and source[0]<0 and a<0:lam*=a/source[0]
    ub=-a*ub;vb=-a*vb+b*p
    tic=time.perf_counter()
    gb,reconstruction_error=general_guard_base(names,guard,gb,a,b,args.tail_margin)
    gb[:861]+=args.source_margin
    exact_A,powers=optimizer_guard_matrix(names,cutoff=10,degree=8)
    A=np.array(exact_A,float);scale=np.maximum(np.max(abs(A),axis=1),1.)
    guard=np.column_stack([guard,(A/scale[:,None]).T])
    gb=np.r_[gb,np.where(np.array(powers)==0,1./scale,0.)-1e-8]
    bank=np.load(BASE/'matrix_roy_guard16x4.npz');assert list(bank['names'])==names
    A=bank['A'];scale=np.maximum(np.max(abs(A),axis=1),1.)
    guard=np.column_stack([guard,(-A/scale[:,None]).T]);gb=np.r_[gb,np.full(len(A),args.source_margin)/scale]
    guard_seconds=time.perf_counter()-tic
    physical=(ub.copy(),vb.copy(),U.copy(),V.copy(),W.copy(),p.copy())
    initial=DiskSupport(*physical).evaluate(lam)[0]
    ng=gb.size;artificial_p=1e5
    ub=np.r_[ub,np.zeros(ng)];vb=np.r_[vb,gb*artificial_p]
    U=np.column_stack([U,np.zeros((nr,ng))]);V=np.column_stack([V,guard*artificial_p])
    W=np.r_[W,np.zeros(ng)];p=np.r_[p,np.full(ng,artificial_p)]
    T=harmonic_tail_transform(names);U=T@U;V=T@V;lam=np.linalg.solve(T.T,lam)
    tic=time.perf_counter()
    lam=(cut_project if args.constraint_generation else feasibility_project)(ub,vb,U,V,W,p,lam)
    projection_seconds=time.perf_counter()-tic
    def checkpoint(x,entry):
        np.savez_compressed(str(args.out)+'_checkpoint.npz',names=np.array(names),multipliers=T.T@x,direction=np.array([a,b]))
        Path(str(args.out)+'_checkpoint.json').write_text(json.dumps({'angle':args.angle,'last_entry':entry,'elapsed_seconds':time.perf_counter()-started,'cache_seconds':cache_seconds,'guard_seconds':guard_seconds,'projection_seconds':projection_seconds},indent=2))
    checkpoint(lam,{'stage':'projected'})
    print('READY',json.dumps({'angle':args.angle,'initial':initial,'projection_seconds':projection_seconds,'guard_seconds':guard_seconds}),flush=True)
    tic=time.perf_counter()
    if args.project_only:diag={'status':'projected_only'}
    else:lam,diag=solve_constrained_dual(ub,vb,U,V,W,p,lam,eps=1e-8,max_iter=args.max_iter,eigen_floor=args.eigen_floor,verbose=True,log_path=str(args.out)+'.history.json',checkpoint=checkpoint)
    optimizer_seconds=time.perf_counter()-tic;lam=T.T@lam
    np.savez_compressed(str(args.out)+'.npz',names=np.array(names),multipliers=lam,direction=np.array([a,b]))
    support=DiskSupport(*physical).evaluate(lam)[0]
    result={'angle':args.angle,'direction':[a,b],'rows':nr,'initial_support_same_grid':initial,'support':support,'optimizer':diag,
            'timing_seconds':{'cache_load':cache_seconds,'guard_build':guard_seconds,'projection':projection_seconds,'optimizer':optimizer_seconds,'total':time.perf_counter()-started},
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'warm_sha256':hashlib.sha256(args.warm.read_bytes()).hexdigest(),'cache_sha256':hashlib.sha256(cache.read_bytes()).hexdigest(),
            'guard_reconstruction_max_error':reconstruction_error,'args':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
            'claim':'Direct finite-spin numerical support, even spins0..40. No continuum or omitted-spin certificate.'}
    Path(str(args.out)+'.json').write_text(json.dumps(result,indent=2));print('FINAL',json.dumps(result),flush=True)

if __name__=='__main__':main()
