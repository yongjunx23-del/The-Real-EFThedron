"""Bernstein coefficient guards for the complete dangerous Roy source region.

Search matrix is a floating image of exact rational polynomial coefficients.
After optimization, --verify recomputes the same cell coefficients exactly.
Successful fixed-cell verification proves W>=margin on the whole region.
"""
from __future__ import annotations
import argparse, hashlib, json, math, re, sys, time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sy


def bernstein_poly(expr, variable, degree):
    p=sy.Poly(sy.expand(expr),variable)
    values=[F(p.nth(j)) for j in range(degree+1)]
    return [sum((values[j]*F(math.comb(i,j),math.comb(degree,j))
                 for j in range(i+1)),F(0)) for i in range(degree+1)]


def exact_global_bank(names):
    specs=[]
    for col,name in enumerate(map(str,names)):
        m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
        if m:specs.append((col,F(m[1]),int(m[2]),int(m[3])))
    endpoints={row[1] for row in specs}
    if len(endpoints)!=1:raise ValueError('One common Roy endpoint is required')
    endpoint=endpoints.pop();ymin=F(8,3);ymax=endpoint-4
    Jmax=max(row[2] for row in specs);nmax=max(row[3] for row in specs)
    da=nmax+Jmax;dq=Jmax
    y,a=sy.symbols('y a');yy=sy.Rational(ymin)+sy.Rational(ymax-ymin)*a
    univariate={}
    for n in range(nmax+1):
        pn=sy.legendre(n,2*y/sy.Rational(ymax)-1)
        for j in range(Jmax+1):
            if j%2:continue
            poly=sy.Poly(sy.expand(pn*(y-sy.Rational(ymin))**j*y**(Jmax-j)/sy.Rational(ymax)**Jmax),y)
            # Polynomial composition preserves all rational coefficients.
            expr=sum(coef*yy**power[0] for power,coef in poly.terms())
            univariate[n,j]=bernstein_poly(expr,a,da)
    bank=np.empty((da+1,dq+1,len(specs)),dtype=object)
    bank.fill(F(0))
    for idx,(_,_,J,n) in enumerate(specs):
        pj=sy.Poly(sy.legendre(J,y),y)
        for j in range(J+1):
            coefficient=F(pj.nth(j))
            if not coefficient:continue
            bb=univariate[n,j]
            for h in range(j,dq+1):
                factor=coefficient*F(math.comb(h,j),math.comb(dq,j))
                for i in range(da+1):bank[i,h,idx]+=factor*bb[i]
    reference=np.empty((da+1,dq+1),dtype=object)
    ref=bernstein_poly((yy/sy.Rational(ymax))**Jmax,a,da)
    for i in range(da+1):reference[i,:]=ref[i]
    return bank,reference,[row[0] for row in specs],dict(endpoint=str(endpoint),
        ymin=str(ymin),ymax=str(ymax),clearing_power=Jmax,degree=[da,dq])


def split_axis(array,axis):
    a=np.moveaxis(array,axis,0)
    degree=a.shape[0]-1
    left=np.empty_like(a);right=np.empty_like(a);row=a.copy()
    left[0]=row[0];right[degree]=row[-1]
    for step in range(1,degree+1):
        row=(row[:-1]+row[1:])/2
        left[step]=row[0];right[degree-step]=row[-1]
    return np.moveaxis(left,0,axis),np.moveaxis(right,0,axis)


def fixed_leaves(array,energy_cells,q_cells):
    if energy_cells&(energy_cells-1) or q_cells&(q_cells-1):
        raise ValueError('Numbers of cells must be powers of two')
    leaves=[array]
    for axis,count in [(0,energy_cells),(1,q_cells)]:
        for _ in range(count.bit_length()-1):
            leaves=[half for leaf in leaves for half in split_axis(leaf,axis)]
    return leaves


def build_guard_cache(witness,out,energy_cells=16,q_cells=4):
    st=np.load(witness);names=list(map(str,st['names']));tic=time.time()
    bank,reference,columns,metadata=exact_global_bank(names)
    print('EXACT_GLOBAL_BANK_SECONDS',time.time()-tic,flush=True)
    # Use the platform's long-double type for the search transform when it
    # offers extra precision. Final feasibility needs the exact checker below.
    bank_float=np.array(bank,dtype=np.longdouble)
    ref_float=np.array(reference,dtype=np.longdouble)
    matrices=fixed_leaves(bank_float,energy_cells,q_cells)
    references=fixed_leaves(ref_float,energy_cells,q_cells)
    reduced=np.concatenate([(m/r[:,:,None]).reshape(-1,len(columns)) for m,r in zip(matrices,references)])
    A=np.zeros((reduced.shape[0],len(names)))
    A[:,columns]=np.array(reduced,dtype=float)
    np.savez_compressed(out,names=np.array(names),A=A,energy_cells=energy_cells,q_cells=q_cells)
    record={'scope':'Numerical search matrix; exact postcheck required','construction':'A@lambda>=margin',
        'domain':'s in [20/3,B], z in [0,1-(8/3)/(s-4)]','metadata':metadata,
        'energy_cells':energy_cells,'q_cells':q_cells,'rows':len(A),
        'source_witness':str(Path(witness).resolve()),'elapsed_seconds':time.time()-tic}
    Path(str(out)+'.json').write_text(json.dumps(record,indent=2),encoding='utf8')
    print(json.dumps(record,indent=2),flush=True)


def verify_fixed(witness,out,energy_cells=16,q_cells=4,margin=F(0)):
    st=np.load(witness);names=list(map(str,st['names']));tic=time.time()
    bank,reference,columns,metadata=exact_global_bank(names)
    lam=[F(float(st['multipliers'][i])) for i in columns]
    total=np.empty(bank.shape[:2],dtype=object)
    for i in range(total.shape[0]):
        for j in range(total.shape[1]):
            total[i,j]=sum((bank[i,j,k]*v for k,v in enumerate(lam)),F(0))
    leaves=fixed_leaves(total,energy_cells,q_cells)
    refs=fixed_leaves(reference,energy_cells,q_cells)
    lower=[];positive=True
    for index,(leaf,ref) in enumerate(zip(leaves,refs)):
        # Ratio bound proves P >= lower*R because R has positive coefficients.
        cell_lower=min(leaf[i,j]/ref[i,j] for i in range(leaf.shape[0]) for j in range(leaf.shape[1]))
        positive=positive and cell_lower>=margin
        ai,qi=divmod(index,q_cells)
        lower.append({'cell_index':index,
                      'a_interval':[str(F(ai,energy_cells)),str(F(ai+1,energy_cells))],
                      'q_interval':[str(F(qi,q_cells)),str(F(qi+1,q_cells))],
                      'W_lower_rational':str(cell_lower),
                      'W_lower_approx':float(cell_lower)})
    minimum=min(F(item['W_lower_rational']) for item in lower)
    result={'scope':'Exact rational whole-domain Roy source sign component certificate; not a complete dual bound',
        'witness':str(Path(witness).resolve()),'sha256':hashlib.sha256(Path(witness).read_bytes()).hexdigest(),
        'metadata':metadata,'energy_cells':energy_cells,'q_cells':q_cells,
        'margin_requested':str(margin),'passed':positive,
        'W_lower_rational':str(minimum),'W_lower_approx':float(minimum),
        'cells':lower,'elapsed_seconds':time.time()-tic}
    Path(out).write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in result.items() if k!='cells'},indent=2),flush=True)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('witness');p.add_argument('--out',required=True)
    p.add_argument('--energy-cells',type=int,default=16);p.add_argument('--q-cells',type=int,default=4)
    p.add_argument('--verify',action='store_true');p.add_argument('--margin',default='0')
    a=p.parse_args()
    if a.verify:verify_fixed(a.witness,a.out,a.energy_cells,a.q_cells,F(a.margin))
    else:build_guard_cache(a.witness,a.out,a.energy_cells,a.q_cells)


if __name__=='__main__':main()
