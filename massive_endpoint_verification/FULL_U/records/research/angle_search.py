"""Checkpointed directional search. Results are candidate supports, not bounds."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import argparse,json,time
from pathlib import Path
import numpy as np
from search_problem import build_problem
from dual_newton import solve_constrained_dual,DiskSupport

ROOT=Path(__file__).resolve().parent.parent

def run(config,angles,out,warm,max_iter=35):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    tic=time.time(); p=build_problem(**config)
    ur=np.array([r[1] for r in p['rows']]);vr=np.array([r[2] for r in p['rows']])
    names=np.array([r[0] for r in p['rows']]); nrows=len(names)
    st=np.load(warm); lookup=dict(zip(st['names'],st['multipliers']))
    lam=np.array([lookup.get(x,0.) for x in names])
    records=[]
    print(json.dumps({'event':'built','seconds':time.time()-tic,'rows':nrows,'points':p['p'].size}),flush=True)
    for theta in angles:
        direction=np.array([np.cos(np.deg2rad(theta)),np.sin(np.deg2rad(theta))])
        direction[np.abs(direction)<1e-14]=0.
        label=f'angle_{theta%360:07.3f}'.replace('.','p')
        npz=out/(label+'.npz'); js=out/(label+'.json')
        if js.exists() and npz.exists():
            record=json.loads(js.read_text());lam=np.load(npz)['multipliers'];records.append(record)
            continue
        Ubase=direction[0]*p['u0'];Vbase=direction[0]*p['v0']+direction[1]*p['p']
        start=time.time()
        if direction[0]==0 and direction[1]<0:
            lam=np.zeros(nrows); result={'support':0.,'status':'analytic_g2_nonnegative',
                'claim':'analytic positivity bound','history':[]}
        else:
            newlam,result=solve_constrained_dual(Ubase,Vbase,ur,vr,p['W'],p['p'],lam,
                eps=2e-6,max_iter=max_iter,ftol=3e-11,gtol=1e-6)
            if not result['status'].startswith('initial_projection'):
                lam=newlam
            else:
                result['unusable']=True
        fun=DiskSupport(Ubase,Vbase,ur,vr,p['W'],p['p']);U,V=fun.fields(lam)
        r=U/np.maximum(np.hypot(U,V),1e-300)
        rad=np.maximum(np.hypot(U,V),1e-300);b=np.empty_like(r);neg=V<0
        b[neg]=r[neg]*U[neg]/(rad[neg]-V[neg]);b[~neg]=1+V[~neg]/rad[~neg]
        g0=float(np.sum(p['W'].ravel()*(p['u0'].ravel()*r+p['v0'].ravel()*b)))
        g2=float(np.sum(p['W'].ravel()*p['p'].ravel()*b))
        record={'theta_deg':float(theta%360),'direction':direction.tolist(),
                'support':float(fun.evaluate(lam)[0]),'g0_phase':g0,'g2_phase':g2,
                'status':result['status'],'seconds':time.time()-start,'config':config,
                'usable_grid_result':not result.get('unusable',False),
                'claim':'finite-grid, finite-spin support; continuum validation pending',
                'optimizer':result}
        np.savez_compressed(npz,multipliers=lam,names=names,direction=direction,
                            s=p['s'],ell=p['ell'],W=p['W'],p=p['p'],U=U.reshape(p['p'].shape),V=V.reshape(p['p'].shape))
        js.write_text(json.dumps(record,indent=2))
        records.append(record)
        (out/'summary.json').write_text(json.dumps(records,indent=2))
        print(json.dumps({k:record[k] for k in ['theta_deg','support','status','seconds','usable_grid_result']}),flush=True)
    return records

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--branch',choices=['positive','negative'],default='positive')
    ap.add_argument('--angles',type=float,nargs='*')
    ap.add_argument('--n',type=int,default=24);ap.add_argument('--lmax',type=int,default=20)
    ap.add_argument('--hmax',type=int,default=20);ap.add_argument('--roy-order',type=int,default=10)
    ap.add_argument('--roy-spin',type=int,default=4);ap.add_argument('--transfer-order',type=int,default=5)
    ap.add_argument('--cross-degree',type=int,default=4)
    ap.add_argument('--radii',type=float,nargs='*',default=[])
    ap.add_argument('--fractional',action='store_true')
    ap.add_argument('--max-iter',type=int,default=35)
    ap.add_argument('--warm',default='output/maxg0_guarded.npz')
    ap.add_argument('--out',default='output/angles')
    args=ap.parse_args()
    config=dict(n=args.n,lmax=args.lmax,hmax=args.hmax,roy_order=args.roy_order,
                roy_spin=args.roy_spin,transfer_order=args.transfer_order,transfer_taper=1,
                cross_degree=args.cross_degree,radii=args.radii,fractional=args.fractional)
    angles=args.angles
    if not angles:
        angles=([0,15,30,45,60,70,75,80,82,84,85,86,87,88,89,90,95,105,120,135,150,165,180]
                if args.branch=='positive' else [0,-15,-30,-45,-60,-75,-90,-105,-120,-135,-150,-165,-180])
    run(config,angles,ROOT/args.out,ROOT/args.warm,args.max_iter)
