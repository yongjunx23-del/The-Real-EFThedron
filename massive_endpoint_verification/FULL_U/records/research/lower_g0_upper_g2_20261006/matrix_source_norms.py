"""Exact rational Bernstein integral bounds for |W K1| and |W K2|."""
import argparse, hashlib, json, math, re, time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sy
from matrix_roy_guards import fixed_leaves


def tensor_bernstein(expression,x,z):
    poly=sy.Poly(sy.expand(expression),x,z);n,m=poly.degree(x),poly.degree(z)
    a=[[F(poly.coeff_monomial(x**i*z**j)) for j in range(m+1)] for i in range(n+1)]
    b=np.empty((n+1,m+1),dtype=object)
    for i in range(n+1):
        for j in range(m+1):
            b[i,j]=sum((a[k][l]*F(math.comb(i,k),math.comb(n,k))*F(math.comb(j,l),math.comb(m,l))
                         for k in range(i+1) for l in range(j+1)),F(0))
    return b


def source_norms(witness,energy_cells=16,angular_cells=16):
    st=np.load(witness);x,z=sy.symbols('x z');c=sy.Rational(4,3);roy={};strip={}
    for name,lam in zip(map(str,st['names']),st['multipliers']):
        v=sy.Rational(float(lam));m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
        if m:
            B=sy.Rational(m[1]);J,n=int(m[2]),int(m[3])
            roy[B]=roy.get(B,0)+2*v*sy.legendre(n,2*x-1)*sy.legendre(J,z)
        m=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
        if m:
            A,B=sy.Rational(m[1]),sy.Rational(m[2]);j,n=int(m[4]),int(m[5])
            assert abs(float(m[3])-4/3)<1e-5
            strip[A,B]=strip.get((A,B),0)+v*(1-(2*x-1)**2)*sy.legendre(n,2*x-1)*sy.legendre(j,2*z-1)
    specs=[]
    for B,W in roy.items():
        s=4+(B-4)*x;t=-(B-4)*x*(1-z)/2
        specs.append((f'Roy_B{B}',W,s,t))
    for (A,B),W in strip.items():
        s=A+(B-A)*x;t=-c*(1-z)
        specs.append((f'Strip_{A}_{B}',W,s,t))
    totals={1:F(0),2:F(0)};signed={1:F(0),2:F(0)};records=[]
    for name,W,s,t in specs:
        for kind,K in [(1,2*(s-c)*(c-(4-s-t))),(2,2*(t-c)**2)]:
            b=tensor_bernstein(W*K,x,z);leaves=fixed_leaves(b,energy_cells,angular_cells)
            signed_value=sum(b.ravel(),F(0))/b.size
            absolute=sum((sum((abs(v) for v in leaf.ravel()),F(0))/leaf.size
                          for leaf in leaves),F(0))/(energy_cells*angular_cells)
            totals[kind]+=absolute;signed[kind]+=signed_value
            records.append({'source':name,'kernel':kind,'degree':[b.shape[0]-1,b.shape[1]-1],
                'absolute_integral_upper':str(absolute),'absolute_integral_upper_approx':float(absolute),
                'signed_integral':str(signed_value),'signed_integral_approx':float(signed_value)})
    return {'witness':str(Path(witness).resolve()),'sha256':hashlib.sha256(Path(witness).read_bytes()).hexdigest(),
        'scope':'Exact rational upper bounds for source integrals of |W K1| and |W K2| on the complete source unit squares',
        'K1':'2(s-4/3)(4/3-u), u=4-s-t','K2':'2(t-4/3)^2',
        'energy_cells':energy_cells,'angular_cells':angular_cells,
        'norm1_rational':str(totals[1]),'norm2_rational':str(totals[2]),
        'norm1_approx':float(totals[1]),'norm2_approx':float(totals[2]),
        'signed1_rational':str(signed[1]),'signed2_rational':str(signed[2]),'records':records}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('witness');p.add_argument('--out',required=True)
    p.add_argument('--energy-cells',type=int,default=16);p.add_argument('--angular-cells',type=int,default=16)
    a=p.parse_args();tic=time.time();r=source_norms(a.witness,a.energy_cells,a.angular_cells)
    r['seconds']=time.time()-tic;Path(a.out).write_text(json.dumps(r,indent=2),encoding='utf8')
    print(json.dumps({k:r[k] for k in ['norm1_approx','norm2_approx','seconds']},indent=2))
