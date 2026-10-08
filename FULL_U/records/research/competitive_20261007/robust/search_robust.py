"""Improve numerical conditioning margins of the 220-row upper-g2 witness.

Numerical search only. All saved coefficients are in the original physical
row basis; stronger pointwise guards do not themselves certify a continuum.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2');os.environ.setdefault('OMP_NUM_THREADS','2')
from pathlib import Path
import sys,argparse,json,time,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'research/full_direction_20261007/right'),str(ROOT/'research/full_direction_20261007/left'),
             str(ROOT/'research'),str(ROOT/'research/directional_certification'),str(ROOT/'research/lower_g0_upper_g2_20261006')]
import numpy as np
from search_upper_g2_guarded import enable_constraint_generation,cut_project,general_guard_base
from angle_audit import normalized_inventory
from numerical_dual_newton import DiskSupport,solve_constrained_dual
from numerical_direction_search import harmonic_tail_transform
from matrix_crossing_series import optimizer_guard_matrix

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--threshold-margin',type=float,default=.001)
    ap.add_argument('--tail-margin',type=float,default=.02);ap.add_argument('--finite-margin',type=float,default=.02)
    ap.add_argument('--source-margin',type=float,default=.0001);ap.add_argument('--max-iter',type=int,default=40)
    ap.add_argument('--warm',type=Path,default=ROOT/'output/full_direction_20261007/right/g2_guarded192_refined.npz')
    ap.add_argument('--cache',type=Path,default=ROOT/'output/full_direction_20261007/right/dense220_n192k96_cache.npz')
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--project-only',action='store_true')
    args=ap.parse_args();args.out.parent.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    enable_constraint_generation();data=np.load(args.cache);names=list(map(str,data['names']));nr=len(names)
    _,_,U,V,W,p,guard,gb=[data[k].copy()for k in ['ub','vb','U','V','W','p','guard','gb']]
    ub=np.zeros_like(p);vb=p.copy();physical=(ub.copy(),vb.copy(),U.copy(),V.copy(),W.copy(),p.copy())
    old=np.load(args.warm);assert list(map(str,old['names']))==names;lam=old['multipliers'].copy()
    load_seconds=time.perf_counter()-start;tic=time.perf_counter()
    gb,reconstruction=general_guard_base(names,guard,gb,0.,1.,args.tail_margin);gb[:861]+=args.source_margin
    metadata=[]
    def add(raw,base,margin,label):
        nonlocal guard,gb
        raw=np.asarray(raw).reshape(nr,-1);base=np.asarray(base).ravel()
        scale=np.maximum(np.max(abs(raw),axis=0),1.)
        guard=np.column_stack([guard,raw/scale]);gb=np.r_[gb,(base+margin)/scale]
        metadata.append(dict(label=label,points=len(base),requested_margin=margin,max_row_scale=float(scale.max())))
    threshold_mu=4+np.array([1e-8,1e-6,1e-4]);threshold_spins=np.array([2,4])
    _,tv0,_,_,tv=normalized_inventory(names,threshold_mu,threshold_spins,96,96)
    add(tv,np.ones(tv0.shape),args.threshold_margin,'threshold_l2_l4')
    densepath=ROOT/'output/full_direction_20261007/left/dense_v_guards_l10_40_mu6p005_20p005_k96.npz'
    dense=np.load(densepath);assert list(map(str,dense['names']))==names
    add(dense['rowVnorm'],np.ones(dense['v0'].shape),args.finite_margin,'dense_all_even_l10_40_mu6p005_20p005')
    AA,powers=optimizer_guard_matrix(names,cutoff=10,degree=8);A=np.asarray(AA,float)
    scale=np.maximum(np.max(abs(A),axis=1),1.)
    guard=np.column_stack([guard,(A/scale[:,None]).T]);gb=np.r_[gb,np.where(np.array(powers)==0,1./scale,0.)-1e-8]
    roybank=ROOT/'output/lower_g0_upper_g2_20261006/matrix_roy_guard16x4.npz';bank=np.load(roybank);assert list(bank['names'])==names
    A=bank['A'];scale=np.maximum(np.max(abs(A),axis=1),1.)
    guard=np.column_stack([guard,(-A/scale[:,None]).T]);gb=np.r_[gb,np.full(len(A),args.source_margin)/scale]
    guard_seconds=time.perf_counter()-tic;initial=DiskSupport(*physical).evaluate(lam)[0]
    ng=len(gb);gp=1e5
    ub=np.r_[ub,np.zeros(ng)];vb=np.r_[vb,gb*gp]
    U=np.column_stack([U,np.zeros((nr,ng))]);V=np.column_stack([V,guard*gp])
    W=np.r_[W,np.zeros(ng)];p=np.r_[p,np.full(ng,gp)]
    T=harmonic_tail_transform(names);U=T@U;V=T@V;lam=np.linalg.solve(T.T,lam)
    tic=time.perf_counter();lam=cut_project(ub,vb,U,V,W,p,lam);projection_seconds=time.perf_counter()-tic
    def save(x,entry):
        np.savez_compressed(str(args.out)+'_checkpoint.npz',names=np.array(names),multipliers=T.T@x,direction=np.array([0.,1.]))
        Path(str(args.out)+'_checkpoint.json').write_text(json.dumps(dict(last_entry=entry,seconds=time.perf_counter()-start),indent=2))
    save(lam,dict(stage='projected'));print('PROJECTED',DiskSupport(*physical).evaluate(T.T@lam)[0],flush=True)
    tic=time.perf_counter()
    if args.project_only:diag={'status':'projection_only'}
    else:lam,diag=solve_constrained_dual(ub,vb,U,V,W,p,lam,eps=1e-8,max_iter=args.max_iter,eigen_floor=1e-7,
                                       verbose=True,log_path=str(args.out)+'.history.json',checkpoint=save)
    optimizer_seconds=time.perf_counter()-tic;lam=T.T@lam
    path=Path(str(args.out)+'.npz');np.savez_compressed(path,names=np.array(names),multipliers=lam,direction=np.array([0.,1.]))
    report=dict(status='NUMERICAL_ONLY',support=DiskSupport(*physical).evaluate(lam)[0],initial_support=initial,rows=nr,optimizer=diag,
        args={k:str(v)if isinstance(v,Path)else v for k,v in vars(args).items()},added_guards=metadata,
        guard_reconstruction_error=reconstruction,artificial_guards=ng,
        timings=dict(load=load_seconds,guard_build=guard_seconds,projection=projection_seconds,optimization=optimizer_seconds,total=time.perf_counter()-start),
        witness_sha256=sha(path),input_sha256={str(p):sha(p)for p in [args.warm,args.cache,densepath,roybank,Path(__file__)]},
        threshold_probe_V_over_p=(1+np.einsum('i,ijk->jk',lam,tv)).tolist(),
        dense_finite_max_V_over_p=float(np.max(1+np.einsum('i,ijk->jk',lam,dense['rowVnorm']))),
        claim='Finite-grid numerical candidate only. Exact component checks and denser head/tail audit required; no complete continuum certificate.')
    Path(str(args.out)+'.json').write_text(json.dumps(report,indent=2));print('FINAL',json.dumps({k:report[k]for k in ['support','initial_support','timings','threshold_probe_V_over_p','dense_finite_max_V_over_p']}),flush=True)

if __name__=='__main__':main()
