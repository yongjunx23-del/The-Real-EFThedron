"""Separate directional numerical refinement, normal (1,b); no full-bound claim is made.

Rebuild saved physical rows on endpoint-split sine-squared velocity grids,
optionally enrich rational resummation, and constrain sampled source signs.
Existing source and final witness files are never modified.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import argparse,hashlib,json,re,time
from pathlib import Path
import numpy as np
from scipy.special import eval_legendre
from angle_audit import normalized_inventory,velocity_grid,fields
from numerical_dual_newton import solve_constrained_dual,DiskSupport


def screen(names,lam):
    s=np.linspace(20/3,11.9,151)[:,None]
    z=np.linspace(0,1,151)[None,:]*(1-(8/3)/(s-4))
    w=np.zeros_like(s+z)
    for name,x in zip(names,lam):
        m=re.fullmatch(r'Roy_B11.9_J(\d+)_n(\d+)',name)
        if m:w+=x*eval_legendre(int(m[2]),2*(s-4)/7.9-1)*eval_legendre(int(m[1]),z)
    return {'source_grid':[151,151],'source_minimum':float(w.min()),
            'sampled_roy_sign_pass':bool(np.all(w>=0)),
            'lambda_X4':float(lam[names.index('X4')]),
            'necessary_crossing_sign_pass':bool(lam[names.index('X4')]<=0)}


def harmonic_tail_transform(names):
    """Invertible row reparametrization, subtracting known low harmonics.

    Saved coefficients always return to the original physical row catalogue.
    This only improves conditioning; it adds no physical constraint.
    """
    from scipy.special import eval_chebyt,eval_chebyu
    T=np.eye(len(names));indices={name:i for i,name in enumerate(names)}
    for i,name in enumerate(names):
        m=re.fullmatch(r'Rc_([^_]+)_([^_]+)',name)
        if m:r,c,kind=float(m[1]),float(m[2]),'cos'
        else:
            m=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',name)
            if not m:continue
            kind,r,c=m[1],float(m[2]),float(m[3])
        denominator=1+2*r*c+r*r
        minusone=(1+r*c)/denominator if kind=='cos' else -r/denominator
        for h,j in indices.items():
            m=re.fullmatch(r'H(\d+)',h)
            if not m:continue
            k=int(m[1])-2
            coefficient=(r**k*eval_chebyt(k,c) if kind=='cos' else
                         (r**k*eval_chebyu(k-1,c) if k else 0.))
            if k==0:coefficient-=minusone
            T[i,j]-=coefficient
    return T


def run(args):
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    if out.with_suffix('.npz').exists():raise FileExistsError(out.with_suffix('.npz'))
    source=Path(args.warm);old=np.load(source);lookup=dict(zip(map(str,old['names']),old['multipliers']))
    # Along a positive change in b, X4/p tends to -9/32 for spin two at
    # threshold. Compensate that constant before feasibility projection.
    # Every coefficient remains free in the actual Newton optimization.
    warm_b=float(old['direction'][1])
    warm_x4_shift=max(args.b-warm_b,0.)*(32/9)
    if 'X4' in lookup:lookup['X4']+=warm_x4_shift
    names=list(lookup)
    for j in range(2,args.hmax+1):
        if f'H{j}' not in names:names.append(f'H{j}')
    if args.roy_order>=0:
        for J in range(0,args.roy_spin+1,2):
            for n in range(args.roy_order+1):
                tag=f'Roy_B11.9_J{J}_n{n}'
                if tag not in names:names.append(tag)
    if args.strip_order>=0:
        for A,B in [(12,24),(24,48)]:
            for j in range(args.strip_j+1):
                for n in range(args.strip_order+1):
                    tag=f'Taper1_Strip_A{A}_B{B}_T{4/3:g}_j{j}_n{n}'
                    if tag not in names:names.append(tag)
    if args.enrich:
        for r,cs in [(.95,[-.95,-.9,-.25,.25,.9,.95]),
                     (.98,[-1,-.95,-.9,-.75,-.5,-.25,0,.25,.5,.75,.9,.95,1])]:
            for c in cs:
                tag=f'Rc_{r:g}_{c:g}'
                if tag not in names:names.append(tag)
        for r in [.95,.98]:
            for c in [-.9,-.75,-.5,-.25,0,.25,.5,.75,.9]:
                tag=f'rational_sin_r{r:g}_c{c:g}'
                if tag not in names:names.append(tag)
    ell=np.arange(0,args.lmax+1,2)
    if args.peak_panels:
        from improvement_peak_audit import peak_velocity_grid
        s,W=peak_velocity_grid(names,args.n,ell)
    else:s,W=velocity_grid(names,args.n,ell)
    # The divided-difference recurrence overflows at delta1e-12 for ell24.
    # Audit delta1e-12 separately at ell<=20; here keep delta1e-10 as before.
    s=np.r_[4+1e-10,s];W=np.column_stack([np.zeros(len(ell)),W])
    t=time.time();u0,v0,p,ur,vn=normalized_inventory(names,s,ell,args.kernel_order,args.angular_order)
    print(json.dumps({'built_seconds':time.time()-t,'rows':len(names),'physical_cells':p.size}),flush=True)
    lam=np.array([lookup.get(name,0.) for name in names]);baseu=u0.ravel();basev=(p*(v0+args.b)).ravel()
    U=ur.reshape(len(names),-1);V=(vn*p).reshape(len(names),-1);ww=W.ravel();pp=p.ravel()
    if getattr(args,'train_bridge',False):
        from improvement_peak_audit import peak_velocity_grid
        bridge_spins=np.arange(args.lmax+2,42,2)
        bs,bw=peak_velocity_grid(names,args.n,bridge_spins,first_order=20)
        bu,bv,bp,bur,bvr=normalized_inventory(names,bs,bridge_spins,args.kernel_order,args.angular_order)
        baseu=np.r_[baseu,bu.ravel()];basev=np.r_[basev,(bp*(bv+args.b)).ravel()]
        U=np.column_stack([U,bur.reshape(len(names),-1)])
        V=np.column_stack([V,(bvr*bp).reshape(len(names),-1)])
        ww=np.r_[ww,bw.ravel()];pp=np.r_[pp,bp.ravel()]
    physical_baseu=baseu.copy();physical_basev=basev.copy()
    physical_U=U.copy();physical_V=V.copy();physical_w=ww.copy();physical_p=pp.copy()
    seed=DiskSupport(baseu,basev,U,V,ww,pp).evaluate(lam)[0]
    # Additional zero-cost cells restrict coefficients; they are not nulls.
    vectors=[];bases=[]
    for ss in np.linspace(20/3,11.9,25):
        for zz in np.linspace(0,1-(8/3)/(ss-4),13):
            vector=[]
            for name in names:
                m=re.fullmatch(r'Roy_B11.9_J(\d+)_n(\d+)',name)
                vector.append(-eval_legendre(int(m[2]),2*(ss-4)/7.9-1)*eval_legendre(int(m[1]),zz) if m else 0.)
            vectors.append(vector);bases.append(0.)
    vectors.append([float(x=='X4') for x in names]);bases.append(0.)
    gs=np.array([4+1e-6,4+.01,4.1,5.,8.,12.1,20.,50.,100.,1000.,1e5])
    ge=np.arange(args.lmax+2,42,2)
    if len(ge):
        _,gv0,_,_,gv=normalized_inventory(names,gs,ge,args.kernel_order,args.angular_order)
        vectors.extend(gv.reshape(len(names),-1).T.tolist());bases.extend((gv0+args.b).ravel().tolist())
    guards=np.array(vectors).T;gb=np.array(bases);scale=np.maximum(np.max(np.abs(guards),axis=0),1.)
    guards/=scale;gb/=scale;ng=len(gb);gp=1e5
    U=np.column_stack([U,np.zeros((len(names),ng))]);V=np.column_stack([V,guards*gp])
    baseu=np.r_[baseu,np.zeros(ng)];basev=np.r_[basev,gb*gp];ww=np.r_[ww,np.zeros(ng)];pp=np.r_[pp,np.full(ng,gp)]
    print(json.dumps({'same_grid_warm_support':seed,'artificial_guards':ng,'warm_screen':screen(names,lam)}),flush=True)
    transform=harmonic_tail_transform(names) if args.strip_rational_harmonics else np.eye(len(names))
    if args.strip_rational_harmonics:
        U=transform@U;V=transform@V;lam=np.linalg.solve(transform.T,lam)
    def checkpoint(coefficients,entry):
        np.savez_compressed(str(out)+'_checkpoint.npz',names=np.array(names),
                            multipliers=transform.T@coefficients,direction=np.array([1.,args.b]))
    lam,result=solve_constrained_dual(baseu,basev,U,V,ww,pp,lam,eps=args.eps,
        max_iter=args.max_iter,eigen_floor=args.eigen_floor,ftol=args.ftol,gtol=args.gtol,
        log_path=out.with_suffix('.history.json'),checkpoint=checkpoint)
    full=DiskSupport(baseu,basev,U,V,ww,pp);cost=full.evaluate(lam)[0]
    lam=transform.T@lam
    physical_cost=DiskSupport(physical_baseu,physical_basev,physical_U,physical_V,physical_w,physical_p).evaluate(lam)[0]
    np.savez_compressed(out.with_suffix('.npz'),names=np.array(names),multipliers=lam,direction=np.array([1.,args.b]))
    config=dict(vars(args))
    if getattr(args,'train_bridge',False):
        config['head_lmax']=args.lmax
        config['lmax']=40
    record={'support':physical_cost,'transformed_objective':cost,
        'parameter_change_roundoff':physical_cost-cost,
        'same_grid_warm_support':seed,'status':result['status'],
        'warm_X4_threshold_compensation':warm_x4_shift,
        'config':config,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'physical_cells':physical_p.size,'artificial_guard_cells':ng,'screen':screen(names,lam),
        'optimizer':result,'claim':'finite-spin quadrature only; independent replay and all-spin certification required'}
    out.with_suffix('.json').write_text(json.dumps(record,indent=2))
    print(json.dumps({k:record[k] for k in ['support','same_grid_warm_support','status','screen']}),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--b',type=float,default=0.)
    ap.add_argument('--warm',default='output/angles_resummed/angle_000p000.npz')
    ap.add_argument('--out',default='output/improvement_20261006/right_enriched')
    ap.add_argument('--n',type=int,default=80);ap.add_argument('--lmax',type=int,default=24)
    ap.add_argument('--hmax',type=int,default=40);ap.add_argument('--enrich',action='store_true')
    ap.add_argument('--roy-order',type=int,default=-1);ap.add_argument('--roy-spin',type=int,default=6)
    ap.add_argument('--strip-order',type=int,default=-1);ap.add_argument('--strip-j',type=int,default=2)
    ap.add_argument('--peak-panels',action='store_true')
    ap.add_argument('--train-bridge',action='store_true',help='Include every even spin26..40 via a separate overflow-safe quadrature')
    ap.add_argument('--strip-rational-harmonics',action='store_true')
    ap.add_argument('--eigen-floor',type=float,default=1e-8)
    ap.add_argument('--ftol',type=float,default=3e-11);ap.add_argument('--gtol',type=float,default=1e-6)
    ap.add_argument('--kernel-order',type=int,default=64);ap.add_argument('--angular-order',type=int,default=56)
    ap.add_argument('--eps',type=float,default=2e-6);ap.add_argument('--max-iter',type=int,default=100)
    run(ap.parse_args())


