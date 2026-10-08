"""Certified polynomial source contractions for an even-Roy witness.

Gauss rules are only used on polynomials, with exactness checked by degree.
The node/weight enclosures come from validated_quadrature.py. No approximate
quadrature assumption enters these moments. Rational geometric truncations
built from the moments still require their own analytic remainder bounds.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,re,json,time,hashlib,argparse
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import arb,arb_mat,ctx
from validated_quadrature import gauss_legendre

def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator

def powers(x,n):
    out=[arb(1)]
    for _ in range(n):out.append(out[-1]*x)
    return out

def legendre(x,n):
    out=[arb(1)]
    if n:out.append(x)
    for k in range(1,n):out.append(((2*k+1)*x*out[-1]-k*out[-2])/(k+1))
    return out

def encode(x):
    return [[str(v) for v in x.lower().man_exp()], [str(v) for v in x.upper().man_exp()]]

def decode(p):
    ends=[arb(int(m))*arb(2)**int(e) for m,e in p]
    return ends[0].union(ends[1])

def build(path,jmax=40,kmax=32,skmax=40,order=None):
    path=Path(path);data=np.load(path);rows=[];ends=set()
    for name,lam in zip(data['names'],data['multipliers']):
        m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',str(name))
        if m:
            B,J,n=F(m[1]),int(m[2]),int(m[3]);ends.add(B)
            assert J%2==0
            rows.append((J,n,F(float(lam))))
    assert len(ends)==1 and rows
    B=next(iter(ends));nmax=max(n for J,n,v in rows);Jmax=max(J for J,n,v in rows)
    degree_s=nmax+jmax+max(kmax-1,skmax)
    degree_z=Jmax+jmax+max(kmax-1,0)
    minorder=(max(degree_s,degree_z)+2)//2
    order=max(48,minorder) if order is None else order
    assert 2*order-1>=max(degree_s,degree_z)
    gx,gw=gauss_legendre(order)
    energy=[legendre(x,nmax) for x in gx]
    angular=[legendre((x+1)/2,Jmax) for x in gx]
    # Contract the energy basis first, retaining all coefficient rounding balls.
    source=[];tscale=(B-4)/2;uscale=F(4);sscale=B-4;tcscale=F(2)
    for i in range(order):
        byJ=[sum((ar(v)*energy[i][n] for J,n,v in rows if J==j),arb(0)) for j in range(Jmax+1)]
        for h in range(order):
            ds=(ar(B)-4)*(gx[i]+1)/2;z=(gx[h]+1)/2
            t=-ds*(1-z)/2;u=-ds-t
            # dx dz times W=2 sum lambda Pn PJ: GL weights give /2.
            weight=gw[i]*gw[h]/2*sum((byJ[j]*angular[h][j] for j in range(Jmax+1)),arb(0))
            source.append((weight,powers(t/ar(tscale),jmax),powers((u+4)/ar(uscale),kmax-1),powers(ds/ar(sscale),skmax),powers((t+2)/ar(tcscale),kmax-1)))
    count=len(source)
    A=arb_mat([[w*tp[j] for w,tp,up,sp,t2p in source] for j in range(jmax+1)])
    U=arb_mat([up for w,tp,up,sp,t2p in source])
    S=arb_mat([sp for w,tp,up,sp,t2p in source])
    mu=A*U;ms=A*S
    # Moment of t^j is column0 of either matrix. t-denominator needs no P_l.
    mt=[sum((w*t2p[k] for w,tp,up,sp,t2p in source),arb(0))*ar(tcscale**k) for k in range(kmax)]
    norm=2*sum((abs(v) for J,n,v in rows),F(0))
    return {'metadata':{'witness_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'precision_bits':ctx.prec,'order':order,'B':str(B),'Jmax':Jmax,'nmax':nmax,
                'jmax':jmax,'kmax':kmax,'skmax':skmax,'energy_degree_bound':degree_s,
                'angular_degree_bound':degree_z,'absolute_W_mass_bound':str(norm),
                'convention':'integrate dx dz with W=2*sum lambda P_n(2x-1)P_J(z), s=4+(B-4)x'},
            'M_u':[[mu[j,k]*ar(tscale**j*uscale**k) for k in range(kmax)] for j in range(jmax+1)],
            'M_s':[[ms[j,k]*ar(tscale**j*sscale**k) for k in range(skmax+1)] for j in range(jmax+1)],
            'T':[mu[j,0]*ar(tscale**j) for j in range(jmax+1)],'M_t':mt}

def save(d,path):
    out={'metadata':d['metadata']}
    for key in ('M_u','M_s'):out[key]=[[encode(x) for x in row] for row in d[key]]
    for key in ('T','M_t'):out[key]=[encode(x) for x in d[key]]
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(out,separators=(',',':')))

def load(path):
    d=json.loads(Path(path).read_text());out={'metadata':d['metadata']}
    for key in ('M_u','M_s'):out[key]=[[decode(x) for x in row] for row in d[key]]
    for key in ('T','M_t'):out[key]=[decode(x) for x in d[key]]
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness',type=Path);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--jmax',type=int,default=40);ap.add_argument('--kmax',type=int,default=32);ap.add_argument('--skmax',type=int,default=40)
    a=ap.parse_args();ctx.prec=220;tic=time.time();d=build(a.witness,a.jmax,a.kmax,a.skmax);save(d,a.out)
    print(json.dumps(d['metadata'],indent=2));print('T0',d['T'][0]);print('seconds',time.time()-tic)
