"""Build an independent dense physical kernel cache, retaining the 220-row guard bank."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import sys,time,json,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research'))
import numpy as np
from angle_audit import normalized_inventory,velocity_grid

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--n',type=int,default=96)
    ap.add_argument('--kernel-order',type=int,default=96)
    args=ap.parse_args();started=time.perf_counter()
    source=ROOT/'output/lower_g0_upper_g2_20261006/lower_enriched_cache.npz'
    old=np.load(source);names=list(map(str,old['names']));nr=len(names)
    records=[];blocks=[]
    for spins,first in [(np.arange(0,26,2),None),(np.arange(26,42,2),20)]:
        tic=time.perf_counter();s,ww=velocity_grid(names,args.n,spins,first_order=first)
        if spins[0]==0:
            s=np.r_[4+1e-10,s];ww=np.column_stack([np.zeros(len(spins)),ww])
        u0,v0,pp,uu,vv=normalized_inventory(names,s,spins,args.kernel_order,args.kernel_order)
        blocks.append((-u0.ravel(),-(v0*pp).ravel(),uu.reshape(nr,-1),(vv*pp).reshape(nr,-1),ww.ravel(),pp.ravel()))
        record={'spins':spins.tolist(),'energy_nodes':len(s),'first_order':first,'seconds':time.perf_counter()-tic}
        records.append(record);print(json.dumps(record),flush=True)
    keys=('ub','vb','U','V','W','p')
    arrays={key:np.concatenate([x[k] for x in blocks],axis=1 if k in (2,3) else 0) for k,key in enumerate(keys)}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.out,names=np.array(names),guard=old['guard'],gb=old['gb'],**arrays)
    meta={'n':args.n,'source_order':args.kernel_order,'angular_order':args.kernel_order,'blocks':records,
          'source_cache':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
          'cache_sha256':hashlib.sha256(args.out.read_bytes()).hexdigest(),'seconds':time.perf_counter()-started,
          'claim':'Numerical physical kernels and inherited sampled guards, not a certificate.'}
    args.out.with_suffix('.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)

if __name__=='__main__':main()
