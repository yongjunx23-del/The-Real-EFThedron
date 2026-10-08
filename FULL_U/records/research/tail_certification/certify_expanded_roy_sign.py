"""Exact rational Bernstein sign check for this witness's dangerous Roy region.

This certifies one component of the all-spin argument, NOT the whole witness.
Witness binary64 multipliers are interpreted as their exact rational values.
"""
import json
import hashlib
import argparse
import math
import re
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as sy


def split_axis(b,axis):
    if axis==1:
        l,r=split_axis(list(map(list,zip(*b))),0)
        return list(map(list,zip(*l))),list(map(list,zip(*r)))
    n=len(b)-1;m=len(b[0])
    left=[[None]*m for _ in range(n+1)]
    right=[[None]*m for _ in range(n+1)]
    for j in range(m):
        row=[b[i][j] for i in range(n+1)]
        left[0][j]=row[0];right[n][j]=row[-1]
        for k in range(1,n+1):
            row=[(a+c)/2 for a,c in zip(row[:-1],row[1:])]
            left[k][j]=row[0];right[n-k][j]=row[-1]
    return left,right


def main():
    root=Path(__file__).resolve().parent.parent
    parser=argparse.ArgumentParser()
    parser.add_argument('witness',nargs='?',type=Path,default=root/'output/maxg0_transfer.npz')
    parser.add_argument('--output',type=Path,default=root/'research/roy_dangerous_sign_certificate.json')
    args=parser.parse_args()
    st=np.load(args.witness)
    y,q,a=sy.symbols('y q a')
    endpoints=set();angular_degrees=[]
    for name in st['names']:
        m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',str(name))
        if m:endpoints.add(m[1]);angular_degrees.append(int(m[2]))
    if len(endpoints)!=1:raise ValueError('Current certificate requires one common Roy endpoint.')
    endpoint=next(iter(endpoints))
    ymax=sy.Rational(endpoint)-4;ymin=sy.Rational(6,5)
    clearing_power=max(angular_degrees)
    if ymax<=ymin:raise ValueError('No dangerous Roy angular region for this endpoint.')
    poly=0
    for name,lam in zip(st['names'],st['multipliers']):
        m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',str(name))
        if not m:continue
        J,n=int(m[2]),int(m[3])
        rational=sy.Rational(float(lam))
        poly+=rational*sy.legendre(n,2*y/ymax-1)*sy.legendre(J,q*(y-ymin)/y)*y**clearing_power
    poly=sy.Poly(sy.cancel(poly).expand(),y,q)
    poly=sy.Poly(poly.as_expr().subs(y,ymin+(ymax-ymin)*a).expand(),a,q)
    da,db=poly.degree(a),poly.degree(q)
    power=[[Fraction(poly.coeff_monomial(a**i*q**j)) for j in range(db+1)]
           for i in range(da+1)]
    bern=[[sum((power[k][l]*Fraction(math.comb(i,k),math.comb(da,k))
                *Fraction(math.comb(j,l),math.comb(db,l))
                for k in range(i+1) for l in range(j+1)),Fraction(0))
           for j in range(db+1)] for i in range(da+1)]
    stack=[(bern,0)];leaves=[];maxdepth=0
    while stack:
        b,depth=stack.pop()
        lower=min(min(row) for row in b)
        if lower>0:
            leaves.append(lower);maxdepth=max(maxdepth,depth);continue
        if depth>=20:raise RuntimeError('Could not certify positive Bernstein coefficients.')
        left,right=split_axis(b,depth%2)
        stack.extend([(left,depth+1),(right,depth+1)])
    lower=min(leaves)/Fraction(ymax**clearing_power)
    result=dict(witness=str(args.witness),witness_sha256=hashlib.sha256(args.witness.read_bytes()).hexdigest(),arithmetic='exact rational (binary64 witness values)',
                domain=f'26/5 <= s <= {endpoint}, 0 <= z <= 1-(6/5)/(s-4)',
                polynomial=f'(s-4)^{clearing_power} sum_Jn lambda_Jn P_n(2(s-4)/{ymax}-1) P_J(z)',
                denominator_clearing_power=clearing_power,
                degree=[da,db],positive_bernstein_leaves=len(leaves),max_subdivision_depth=maxdepth,
                W_certified_lower_bound=float(lower),
                W_lower_bound_rational=f'{lower.numerator}/{lower.denominator}',
                conclusion='W>0 on t<=-3/5; this is a source-sign component only.',
                limitation='This is NOT a full all-spin V/p upper bound or an exclusion certificate.')
    args.output.write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()

