"""Independent fixed-multiplier audits for arbitrary massive coupling normals.

Supports H_j, Rc, X4..X8 (and higher generated crossing rows), Roy, Strip,
and TaperN_Strip inventories. No optimization occurs. The spin head and
bridge use endpoint-split velocity quadrature; omitted spins remain unproved.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import time
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from roy_kernels import build_roy_bank
from transfer_kernels import build_transfer_bank
from crossing_kernels import build_crossing_bank
from existing_solver.massive_forward_kernels import coupling_kernels,forward_geometry,rational_generator,fractional_difference
from numpy.polynomial.legendre import legmul

_CACHE={}


def parse_strip(name):
    m=re.fullmatch(r'(?:Taper(\d+)_)?Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
    if not m:return None
    power=int(m[1] or 0);A,B,T=float(m[2]),float(m[3]),float(m[4])
    if abs(T-4/3)<1e-5:T=4/3  # The stored :g label truncates canonical 4/3.
    return power,(A,B,T,int(m[5]),int(m[6]))


def taper_expansion(spec,power):
    A,B,T,j,n=spec;coeff=np.zeros(n+1);coeff[n]=1
    for _ in range(power):coeff=legmul(coeff,[2/3,0,-2/3])
    return [((A,B,T,j,k),float(c)) for k,c in enumerate(coeff) if c]


def inventory_breakpoints(names):
    edges=set()
    for name in names:
        m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
        if m:edges.add(float(m[1]))
        st=parse_strip(name)
        if st:edges.update(st[1][:2])
    return sorted(edges)


def velocity_grid(names,n,spins,first_order=None):
    edges=sorted(set([0.,.12,.35,.95,.995,1.]+[
        math.sqrt(1-4/e) for e in inventory_breakpoints(names)]))
    ss=[];ws=[]
    for a,b in zip(edges[:-1],edges[1:]):
        x,w=roots_legendre(first_order if a==0 and first_order else n)
        if a==0:
            beta=(a+b)/2+(b-a)/2*x;wb=(b-a)/2*w
        else:
            angle=np.pi*(x+1)/4
            beta=a+(b-a)*np.sin(angle)**2
            wb=w*(b-a)*np.pi/4*np.sin(2*angle)
        ss.append(4/(1-beta*beta))
        ws.append(4/np.pi*wb/(1-beta*beta)**2)
    return np.concatenate(ss),(2*np.asarray(spins)[:,None]+1)*np.concatenate(ws)


def normalized_inventory(names,s,spins,kernel_order=64,angular_order=56):
    """Return u0, v0/p, p, row U kernels, and row V/p kernels.

    Independent assembly sums normalized V terms before multiplying by p;
    this avoids cancellation between separate enormous threshold kernels.
    """
    names=tuple(map(str,names));spins=np.asarray(spins,int);s=np.asarray(s,float)
    key=(names,tuple(spins),hashlib.sha256(s.tobytes()).hexdigest(),kernel_order,angular_order)
    if key in _CACHE:return _CACHE[key]
    cu,cv=coupling_kernels(s,spins);p=cv[1];sqrtp=np.sqrt(p)
    u0=cu[0];v0=cv[0]/p
    u=np.zeros((len(names),len(spins),len(s)));v=np.zeros_like(u)
    F,omega,_=forward_geometry(s)
    roy_specs=[];roy_inds=[];strip_specs=[];strip_map={};strip_rows=[];cross_inds=[];deriv_inds=[]
    for i,name in enumerate(names):
        m=re.fullmatch(r'H(\d+)',name)
        if m:
            q=F*omega**int(m[1]);u[i]=q.imag;v[i]=q.real/p;continue
        m=re.fullmatch(r'Rc_([^_]+)_([^_]+)',name)
        if m:
            q=rational_generator(s,float(m[1]),float(m[2]),'cos');u[i]=q.imag;v[i]=q.real/p;continue
        m=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',name)
        if m:
            q=rational_generator(s,float(m[2]),float(m[3]),m[1]);u[i]=q.imag;v[i]=q.real/p;continue
        m=re.fullmatch(r'Delta(\d+)_C(\d+)_a([^_]+)_h([^_]+)',name)
        if m:
            if int(m[1])!=int(m[2])+1:raise ValueError('Unexpected fractional difference order')
            q=fractional_difference(s,int(m[2]),float(m[3]),float(m[4]));u[i]=q.imag;v[i]=q.real/p;continue
        if re.fullmatch(r'Dt[123](?:H_j\d+|G_regular)',name):deriv_inds.append(i);continue
        if re.fullmatch(r'X\d+(?:_\d+)?',name):cross_inds.append(i);continue
        m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
        if m:
            roy_specs.append((float(m[1]),int(m[2]),int(m[3])));roy_inds.append(i);continue
        st=parse_strip(name)
        if st:
            terms=[]
            for spec,c in taper_expansion(st[1],st[0]):
                if spec not in strip_map:strip_map[spec]=len(strip_specs);strip_specs.append(spec)
                terms.append((strip_map[spec],c))
            strip_rows.append((i,terms));continue
        raise ValueError(f'Unknown physical row {name!r}')
    if cross_inds:
        degree=max(int(re.match(r'X(\d+)',names[i])[1]) for i in cross_inds)
        bank=build_crossing_bank(s,spins,degree,normalized=True)
        lookup={str(label):k for k,label in enumerate(bank['labels'])}
        for i in cross_inds:v[i]=bank['v'][lookup[names[i]]]
    if deriv_inds:
        from forward_transfer_derivatives import build_derivative_bank
        order=max(int(names[i][2]) for i in deriv_inds)
        hmax=max([int(names[i].split('_j')[1]) for i in deriv_inds if '_j' in names[i]]+[2])
        bank={name:(uu,vv) for name,uu,vv in build_derivative_bank(s,spins,order,hmax)}
        for i in deriv_inds:
            u[i]=bank[names[i]][0];v[i]=bank[names[i]][1]/p
    if roy_specs:
        bank=build_roy_bank(s,spins,roy_specs,kernel_order,angular_order,normalized=True)
        for k,i in enumerate(roy_inds):
            u[i]=bank['u'][k]*sqrtp+bank['g0_coeff'][k]*u0
            v[i]=bank['v'][k]+bank['g0_coeff'][k]*v0
    if strip_specs:
        bank=build_transfer_bank(s,spins,strip_specs,kernel_order,angular_order,normalized=True)
        for i,terms in strip_rows:
            for k,c in terms:
                u[i]+=c*(bank['u'][k]*sqrtp+bank['g0_coeff'][k]*u0)
                v[i]+=c*(bank['v'][k]+bank['g0_coeff'][k]*v0)
    result=(u0,v0,p,u,v)
    # Bound cache growth in long campaigns with changing inventories.
    if len(_CACHE)>12:_CACHE.clear()
    _CACHE[key]=result
    return result


def fields(names,lam,direction,s,spins,**kwargs):
    u0,v0,p,u,v=normalized_inventory(names,s,spins,**kwargs)
    U=direction[0]*u0+np.einsum('i,ijk->jk',lam,u)
    vnorm=direction[0]*v0+direction[1]+np.einsum('i,ijk->jk',lam,v)
    return U,vnorm,p


def panel_support(names,lam,direction,n,spins,first_order=None,**kwargs):
    s,W=velocity_grid(names,n,spins,first_order)
    U,v,p=fields(names,lam,direction,s,spins,**kwargs)
    V=p*v;radius=np.hypot(U,V)
    h=radius+V
    np.multiply(U,np.divide(U,radius-V,out=np.zeros_like(U),where=V<0),out=h,where=V<0)
    costs=np.sum(W*h,axis=1)
    positive_large_p=[]
    for i,l in enumerate(spins):
        bad=(p[i]>1e4)&(v[i]>0)
        if l and np.any(bad):
            k=int(np.argmax(np.where(bad,W[i]*h[i],-1)))
            positive_large_p.append(dict(spin=int(l),energy=float(s[k]),
                                         V_over_p=float(v[i,k]),p=float(p[i,k]),
                                         node_cost=float(W[i,k]*h[i,k])))
    return costs,dict(max_V_over_p={str(l):float(np.max(v[i])) for i,l in enumerate(spins)},
                      positive_large_p_samples=positive_large_p,
                      energy_nodes=len(s),s_min=float(s.min()),s_max=float(s.max()),
                      first_panel_order=first_order or n,
                      finite=bool(np.isfinite(costs).all()))


def audit(path,n=80,lmax=40,head_lmax=24,kernel_order=64,angular_order=56):
    tic=time.time();path=Path(path);data=np.load(path)
    names=tuple(map(str,data['names']));lam=np.asarray(data['multipliers']);normal=np.asarray(data['direction'],float)
    opts=dict(kernel_order=kernel_order,angular_order=angular_order)
    head_spins=np.arange(0,min(head_lmax,lmax)+1,2)
    head,head_info=panel_support(names,lam,normal,n,head_spins,**opts)
    bridge_spins=np.arange(head_lmax+2,lmax+1,2)
    if len(bridge_spins):
        bridge,bridge_info=panel_support(names,lam,normal,n,bridge_spins,first_order=min(n,20),**opts)
    else:bridge=np.array([]);bridge_info={}
    # Zero-weight threshold probes are never integrated or fitted away.
    deltas=np.array([1e-4,1e-6,1e-8,1e-10,1e-12]);probe_spins=np.arange(0,min(lmax,20)+1,2)
    _,vp,_=fields(names,lam,normal,4+deltas,probe_spins,**opts)
    limiting=vp[1:,-1]
    threshold_pass=bool(np.isfinite(limiting).all() and np.all(limiting<=0.0))
    # Raw P_l at delta=1e-12 would overflow beyond the low-spin probe.
    # This separately named finite-offset test reaches all audited spins.
    extended_spins=np.arange(2,lmax+1,2)
    extended_deltas=np.array([1e-2,1e-4,1e-6])
    _,extended_vp,_=fields(names,lam,normal,4+extended_deltas,extended_spins,**opts)
    extended_pass=bool(np.isfinite(extended_vp[:,-1]).all() and np.all(extended_vp[:,-1]<=0.0))
    total=float(head.sum()+bridge.sum())
    source_grid=None;source_angle=None;training_lmax=None
    adjacent=path.with_suffix('.json')
    if adjacent.exists():
        old=json.loads(adjacent.read_text(encoding='utf8'))
        source_grid=old.get('dual_support_grid',old.get('support'))
        source_angle=old.get('theta_deg')
        training_lmax=old.get('config',{}).get('lmax')
    angle=float(np.degrees(np.arctan2(normal[1],normal[0])))
    all_spins=np.r_[head_spins,bridge_spins];all_costs=np.r_[head,bridge]
    training_replay=float(all_costs[all_spins<=training_lmax].sum()) if training_lmax is not None else None
    result=dict(witness=str(path.resolve()),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                direction=normal.tolist(),angle_degrees=angle,angle_0_360=angle%360,
                source_angle_degrees=source_angle,
                support=total if np.isfinite(total) else None,
                source_grid_support=source_grid,
                difference_from_source_grid=total-source_grid if source_grid is not None and np.isfinite(total) else None,
                training_lmax=training_lmax,support_through_training_lmax=training_replay,
                same_spin_difference_from_source_grid=training_replay-source_grid if source_grid is not None and training_replay is not None else None,
                per_spin={str(l):float(c) for l,c in zip(np.r_[head_spins,bridge_spins],np.r_[head,bridge])},
                n=n,lmax=lmax,head_lmax=min(head_lmax,lmax),head=head_info,bridge=bridge_info,
                kernel_order=kernel_order,angular_order=angular_order,
                threshold_probe_deltas=deltas.tolist(),threshold_probe_spins=probe_spins.tolist(),
                threshold_V_over_p=vp.tolist(),threshold_sign_tolerance=0.,
                threshold_pass=threshold_pass,threshold_max_nonzero_spin=float(np.max(limiting)),
                extended_threshold_probe_spins=extended_spins.tolist(),
                extended_threshold_probe_deltas=extended_deltas.tolist(),
                extended_threshold_V_over_p=extended_vp.tolist(),
                extended_threshold_probe_pass=extended_pass,
                extended_threshold_max=float(np.max(extended_vp[:,-1])),
                all_spin_certified=False,continuum_certified=False,
                status='finite_head_audited_tail_unproved' if threshold_pass else 'failed_threshold_sign',
                row_count=len(names),elapsed_seconds=time.time()-tic)
    result['large_p_guard_pass']=not (head_info['positive_large_p_samples'] or bridge_info.get('positive_large_p_samples',[]))
    fractional_alphas=[float(re.fullmatch(r'Delta\d+_C\d+_a([^_]+)_h[^_]+',name)[1])
                       for name,coefficient in zip(names,lam) if name.startswith('Delta') and coefficient!=0]
    singular_fractional=bool(fractional_alphas and min(fractional_alphas)<0)
    result['fractional_threshold_treatment']=('ordinary_velocity_Gauss_unresolved' if singular_fractional else
                  'bounded_fractional_rows_velocity_Gauss' if fractional_alphas else 'not_applicable')
    if singular_fractional:result['status']='unresolved_fractional_threshold'
    if not result['large_p_guard_pass']:result['status']='failed_dense_large_p_guard'
    if not extended_pass:result['status']='failed_extended_threshold_probe'
    if not head_info['finite'] or (bridge_info and not bridge_info['finite']):result['status']='failed_nonfinite_kernel_or_support'
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('witnesses',nargs='*',type=Path)
    ap.add_argument('--directory',type=Path)
    ap.add_argument('--pattern',default='*.npz')
    ap.add_argument('--n',type=int,default=80);ap.add_argument('--lmax',type=int,default=40)
    ap.add_argument('--head-lmax',type=int,default=24)
    ap.add_argument('--kernel-order',type=int,default=64);ap.add_argument('--angular-order',type=int,default=56)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();paths=list(args.witnesses)
    if args.directory:paths.extend(sorted(args.directory.glob(args.pattern)))
    args.output.mkdir(parents=True,exist_ok=True);results=[]
    duplicate_stems={path.stem for path in paths if sum(p.stem==path.stem for p in paths)>1}
    for path in paths:
        result=audit(path,args.n,args.lmax,args.head_lmax,args.kernel_order,args.angular_order)
        label=f'{path.parent.name}_{path.stem}' if path.stem in duplicate_stems else path.stem
        (args.output/f'{label}_audit.json').write_text(json.dumps(result,indent=2),encoding='utf8')
        results.append(result)
        print(json.dumps({k:result[k] for k in ['witness','angle_degrees','support','threshold_pass','status','elapsed_seconds']}),flush=True)
        keys=['witness','angle_degrees','angle_0_360','source_angle_degrees','support','source_grid_support','difference_from_source_grid','training_lmax','support_through_training_lmax','same_spin_difference_from_source_grid','threshold_pass','threshold_max_nonzero_spin','extended_threshold_probe_pass','extended_threshold_max','large_p_guard_pass','status','all_spin_certified']
        with (args.output/'summary.csv').open('w',newline='',encoding='utf8') as stream:
            writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader()
            writer.writerows({k:r[k] for k in keys} for r in results)


if __name__=='__main__':main()
