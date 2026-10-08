"""New lower-g0 search. Grid objectives and sign screens are not certificates."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
sys.path.insert(0,str(HERE.parent/'directional_certification'))
import argparse,json,re,time,hashlib
import numpy as np
from scipy.special import eval_legendre
from angle_audit import normalized_inventory,velocity_grid
from numerical_dual_newton import DiskSupport,solve_constrained_dual
from numerical_direction_search import harmonic_tail_transform
from crossing_kernels import crossing_catalogue
from scipy.optimize import linprog

def feasibility_project(ub,vb,U,V,W,p,lam,eps=1e-8):
    fun=DiskSupport(ub,vb,U,V,W,p)
    mask=p>1e4;C=fun.vnorm[:,mask].T;target=-eps-fun.vb[mask]
    scale=np.maximum(np.max(np.abs(C),axis=1),1.)
    C=C/scale[:,None];target=target/scale-1e-10
    column=np.maximum(np.max(np.abs(C),axis=0),1e-12)
    A=C/column[None,:];rhs=target-C@lam;nr=lam.size
    sol=linprog(np.r_[np.zeros(nr),np.ones(nr)],
       A_ub=np.vstack([np.c_[A,np.zeros_like(A)],np.c_[np.eye(nr),-np.eye(nr)],np.c_[-np.eye(nr),-np.eye(nr)]]),
       b_ub=np.r_[rhs,np.zeros(2*nr)],bounds=[(None,None)]*nr+[(0,None)]*nr,
       method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    print('LP_PROJECT',sol.success,sol.message,flush=True)
    if not sol.success:return lam
    result=lam+sol.x[:nr]/column
    print('PROJECT_SUPPORT',fun.evaluate(result)[0],'SCALED_VIOLATION',np.max(C@result-target),flush=True)
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--warm',default='output/angle_constraints/ming0_scaled_threshold_repaired.npz')
    ap.add_argument('--out',default='output/lower_g0_upper_g2_20261006/lower_enriched')
    ap.add_argument('--cache',default='output/lower_g0_upper_g2_20261006/lower_enriched_cache.npz')
    ap.add_argument('--n',type=int,default=48)
    ap.add_argument('--max-iter',type=int,default=80)
    ap.add_argument('--tail-margin',type=float,default=1e-4)
    ap.add_argument('--source-margin',type=float,default=1e-4)
    ap.add_argument('--crossing-matrix',action='store_true')
    ap.add_argument('--source-bernstein',type=Path)
    ap.add_argument('--project-only',action='store_true')
    ap.add_argument('--eigen-floor',type=float,default=1e-7)
    args=ap.parse_args()
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    old=np.load(args.warm);lookup=dict(zip(map(str,old['names']),old['multipliers']))
    names=list(lookup)
    extra=[f'H{j}' for j in range(2,41)]
    extra += [f'Roy_B11.9_J{J}_n{n}' for J in range(0,7,2) for n in range(15)]
    extra += [f'Taper1_Strip_A{A}_B{B}_T{4/3:g}_j{j}_n{n}' for A,B in [(12,24),(24,48)] for j in range(4) for n in range(8)]
    extra += [f'rational_{kind}_r0.8_c{c:g}' for kind in ['cos','sin'] for c in [-.9,-.75,-.5,0,.5,.75,.9]]
    for name in extra:
        if name not in names:names.append(name)
    nr=len(names);lam=np.array([lookup.get(name,0.) for name in names])
    cache=Path(args.cache)
    if cache.exists():
        data=np.load(cache)
        assert list(data['names'])==names
        ub,vb,U,V,W,p,guard,gb=[data[k] for k in ['ub','vb','U','V','W','p','guard','gb']]
        print('LOADED_CACHE',cache,flush=True)
    else:
        tic=time.time();blocks=[]
        for spins,n,first in [(np.arange(0,26,2),args.n,None),(np.arange(26,42,2),args.n,20)]:
            s,ww=velocity_grid(names,n,spins,first_order=first)
            if spins[0]==0:s=np.r_[4+1e-10,s];ww=np.column_stack([np.zeros(len(spins)),ww])
            u0,v0,pp,uu,vv=normalized_inventory(names,s,spins,64,64)
            blocks.append((-u0.ravel(),-(v0*pp).ravel(),uu.reshape(nr,-1),(vv*pp).reshape(nr,-1),ww.ravel(),pp.ravel()))
            print('PHYSICAL_BLOCK',spins.tolist(),s.size,'seconds',time.time()-tic,flush=True)
        ub,vb,U,V,W,p=[np.concatenate([x[k] for x in blocks],axis=1 if k in (2,3) else 0) for k in range(6)]
        vectors=[];bases=[]
        # Positive Roy source weight prevents exponential bad-sign growth.
        for ss in np.linspace(20/3,11.9,41):
            for zz in np.linspace(0,1-(8/3)/(ss-4),21):
                row=[]
                for name in names:
                    m=re.fullmatch(r'Roy_B11.9_J(\d+)_n(\d+)',name)
                    row.append(-eval_legendre(int(m[2]),2*(ss-4)/7.9-1)*eval_legendre(int(m[1]),zz) if m else 0.)
                vectors.append(row);bases.append(0.)
        cat=crossing_catalogue(8);lead={r['name']:float(r['beta'][6]) if len(r['beta'])>6 else 0. for r in cat}
        vectors.append([lead.get(name,0.) for name in names]);bases.append(0.)
        grids=[(np.array([4+1e-6,4+.001,4+.01,4.1,4.5,5.,6.,8.,12.1,20.,30.,50.,75.,100.,150.,250.,500.,1000.,1e4,1e6]),np.arange(26,42,2),0.),
               (np.r_[np.arange(13.,105.,2.),150.,250.,500.,1000.,1e4,1e6],np.arange(42,74,2),args.tail_margin),
               (np.array([8.,10.,12.1,16.,20.,30.,50.,75.,100.,150.,250.,500.,1000.,1e4,1e6]),np.array([80,100,140,200]),args.tail_margin)]
        tau=np.unique(np.r_[np.geomspace(.001,1,12),np.arange(1,7.001,.125),np.arange(7,40.001,.25),np.geomspace(40,150,14)])
        for ell in np.r_[np.arange(42,102,2),[120,160,200,300,500,1000]]:
            grids.append((4+ell*(ell+1)/tau,np.array([ell]),args.tail_margin))
        for k,(ss,ll,margin) in enumerate(grids):
            _,v0,_,_,vv=normalized_inventory(names,ss,ll,64,64)
            assert np.isfinite(vv).all()
            vectors.extend(vv.reshape(nr,-1).T.tolist());bases.extend((-v0+margin).ravel().tolist())
            if k%5==0:print('GUARD_GRID',k,'of',len(grids),'seconds',time.time()-tic,flush=True)
        guard=np.array(vectors).T;gb=np.array(bases)
        scale=np.maximum(np.max(np.abs(guard),axis=0),1.)
        guard/=scale;gb/=scale
        np.savez_compressed(cache,names=np.array(names),ub=ub,vb=vb,U=U,V=V,W=W,p=p,guard=guard,gb=gb)
    gb[:41*21]+=args.source_margin
    if args.crossing_matrix:
        from matrix_crossing_series import optimizer_guard_matrix
        exact_A,powers=optimizer_guard_matrix(names,cutoff=10,degree=8)
        A=np.array(exact_A,float);scale=np.maximum(np.max(np.abs(A),axis=1),1.)
        # Global eps is cancelled on nonconstant powers; q0<=-1.
        base=np.where(np.array(powers)==0,1./scale,0.)-1e-8
        guard=np.column_stack([guard,(A/scale[:,None]).T]);gb=np.r_[gb,base]
    if args.source_bernstein:
        bank=np.load(args.source_bernstein);assert list(bank['names'])==names
        A=bank['A'];scale=np.maximum(np.max(np.abs(A),axis=1),1.)
        guard=np.column_stack([guard,(-A/scale[:,None]).T])
        gb=np.r_[gb,np.full(len(A),args.source_margin)/scale]
    initial=DiskSupport(ub,vb,U,V,W,p).evaluate(lam)[0]
    print('WARM_SUPPORT',initial,'ROWS',nr,'GUARDS',gb.size,flush=True)
    physical=(ub.copy(),vb.copy(),U.copy(),V.copy(),W.copy(),p.copy())
    artificial_p=1e5;ng=gb.size
    ub=np.r_[ub,np.zeros(ng)];vb=np.r_[vb,gb*artificial_p]
    U=np.column_stack([U,np.zeros((nr,ng))]);V=np.column_stack([V,guard*artificial_p])
    W=np.r_[W,np.zeros(ng)];p=np.r_[p,np.full(ng,artificial_p)]
    T=harmonic_tail_transform(names)
    U=T@U;V=T@V;lam=np.linalg.solve(T.T,lam)
    lam=feasibility_project(ub,vb,U,V,W,p,lam)
    def checkpoint(x,entry):
        np.savez_compressed(str(out)+'_checkpoint.npz',names=np.array(names),multipliers=T.T@x,direction=np.array([-1.,0.]))
    checkpoint(lam,{'stage':'projected'})
    if args.project_only:diag={'status':'projected_only','support':DiskSupport(ub,vb,U,V,W,p).evaluate(lam)[0]}
    else:
        lam,diag=solve_constrained_dual(ub,vb,U,V,W,p,lam,eps=1e-8,max_iter=args.max_iter,eigen_floor=args.eigen_floor,
                 verbose=True,log_path=str(out)+'.history.json',checkpoint=checkpoint)
    lam=T.T@lam
    np.savez_compressed(str(out)+'.npz',names=np.array(names),multipliers=lam,direction=np.array([-1.,0.]))
    objective=DiskSupport(*physical).evaluate(lam)[0]
    record={'warm_support_same_grid':initial,'head_support':objective,'status':diag['status'],
        'source':args.warm,'source_sha256':hashlib.sha256(Path(args.warm).read_bytes()).hexdigest(),
        'configuration':{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},'rows':nr,'guard_count':ng,'optimizer':diag,
        'claim':'Numerical all-energy quadrature for even spins0..40; no omitted-spin or interval certificate.'}
    Path(str(out)+'.json').write_text(json.dumps(record,indent=2))
    print('FINAL',json.dumps({k:record[k] for k in ['warm_support_same_grid','head_support','status']}),flush=True)
if __name__=='__main__':main()
