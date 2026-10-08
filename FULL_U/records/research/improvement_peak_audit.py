"""Fixed-witness replay with explicit rational-peak panel boundaries.

This refines numerical quadrature only, not continuum/all-spin certification.
The first beta panel stays [0,.12] to retain the existing overflow-safe
head/bridge behavior and to avoid weakening its threshold coverage.
"""
from __future__ import annotations
import argparse,json,math,re
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
import angle_audit as aa
HEAD_FIRST_ORDER=80


def peak_edges(names):
    edges=set([0.,.12,.35,.95,.995,1.]+[
        math.sqrt(1-4/e) for e in aa.inventory_breakpoints(names)])
    for name in names:
        m=re.fullmatch(r'Rc_([^_]+)_([^_]+)',str(name))
        if m:r,c=float(m[1]),float(m[2])
        else:
            m=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',str(name))
            if not m:continue
            r,c=float(m[2]),float(m[3])
        if r<.94 or abs(c)>=1:continue
        center=math.tan(math.acos(c)/4)
        width=(1+center*center)*(-math.log(r))/4
        for value in [center-2*width,center,center+2*width]:
            if .12<value<.995:edges.add(value)
    return sorted(edges)


def peak_velocity_grid(names,n,spins,first_order=None):
    edges=peak_edges(names);ss=[];ws=[]
    for a,b in zip(edges[:-1],edges[1:]):
        # Keep head first-panel refinement comparable to the original n80/160.
        order=(first_order or max(HEAD_FIRST_ORDER,n)) if a==0 else n
        x,w=roots_legendre(order)
        if a==0:beta=(a+b)/2+(b-a)/2*x;wb=(b-a)/2*w
        else:
            angle=np.pi*(x+1)/4;beta=a+(b-a)*np.sin(angle)**2
            wb=w*(b-a)*np.pi/4*np.sin(2*angle)
        ss.append(4/(1-beta**2));ws.append(4/np.pi*wb/(1-beta**2)**2)
    return np.concatenate(ss),(2*np.asarray(spins)[:,None]+1)*np.concatenate(ws)


def main():
    global HEAD_FIRST_ORDER
    ap=argparse.ArgumentParser();ap.add_argument('witness',type=Path)
    ap.add_argument('--n',type=int,default=24);ap.add_argument('--lmax',type=int,default=40)
    ap.add_argument('--head-lmax',type=int,default=24)
    ap.add_argument('--threshold-order',type=int,default=80)
    ap.add_argument('--kernel-order',type=int,default=80);ap.add_argument('--angular-order',type=int,default=80)
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    HEAD_FIRST_ORDER=args.threshold_order
    names=np.load(args.witness)['names'];edges=peak_edges(names)
    aa.velocity_grid=peak_velocity_grid
    record=aa.audit(args.witness,args.n,args.lmax,args.head_lmax,args.kernel_order,args.angular_order)
    record['quadrature']='Rational-peak split sine-squared panels; first panel linear, bridge20'
    record['head_first_panel_order']=max(args.threshold_order,args.n)
    record['velocity_edges']=edges;record['panel_count']=len(edges)-1
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(record,indent=2))
    print(json.dumps({k:record[k] for k in ['support','status','threshold_pass','large_p_guard_pass','panel_count','elapsed_seconds']}),flush=True)


if __name__=='__main__':main()
