"""Arb-enclosed polynomial source moments for the power-one strip rows.

Every integral is polynomial and the Gauss degree criterion is asserted.
Uses canonical T=4/3, exact decimal A,B, and exact binary64 multipliers.
"""
from certified_roy_moments import ROOT,ar,powers,legendre,save
from validated_quadrature import gauss_legendre
from flint import arb,arb_mat,ctx
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,re,time
import numpy as np

def build(path,jmax=40,kmax=32,order=None):
    path=Path(path);data=np.load(path);groups={}
    for name,lam in zip(data['names'],data['multipliers']):
        m=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',str(name))
        if m:
            A,B=F(m[1]),F(m[2]);assert abs(float(m[3])-4/3)<1e-5
            groups.setdefault((A,B),[]).append((int(m[4]),int(m[5]),F(float(lam))))
    results=[];T=F(4,3)
    for (A,B),rows in groups.items():
        nmax=max(n for j,n,v in rows);Jmax=max(j for j,n,v in rows)
        de=nmax+2+kmax-1;da=Jmax+jmax+kmax-1
        size=max(40,(max(de,da)+2)//2) if order is None else order
        assert 2*size-1>=max(de,da)
        gx,gw=gauss_legendre(size);energy=[legendre(x,nmax) for x in gx];angular=[legendre(x,Jmax) for x in gx]
        h=(A+B)/2-4-T/2;radius=(B-A)/2+T/2;smid=(A+B)/2;sradius=(B-A)/2;source=[]
        for i in range(size):
            byJ=[sum((ar(v)*energy[i][n] for j,n,v in rows if j==jj),arb(0)) for jj in range(Jmax+1)]
            s=ar(A)+ar(B-A)*(gx[i]+1)/2;taper=1-gx[i]*gx[i]
            for k in range(size):
                t=ar(T)*(gx[k]-1)/2;u=4-s-t
                w=gw[i]*gw[k]/4*taper*sum((byJ[j]*angular[k][j] for j in range(Jmax+1)),arb(0))
                source.append((w,powers(t/ar(T),jmax),powers((u+ar(h))/ar(radius),kmax-1),powers((t+ar(T/2))/ar(T/2),kmax-1),powers((s-ar(smid))/ar(sradius),kmax-1)))
        aa=arb_mat([[w*tp[j] for w,tp,up,t2p,sp in source] for j in range(jmax+1)])
        uu=arb_mat([up for w,tp,up,t2p,sp in source]);mu=aa*uu
        ss=arb_mat([sp for w,tp,up,t2p,sp in source]);ms=aa*ss
        mt=[sum((w*t2p[k] for w,tp,up,t2p,sp in source),arb(0))*ar((T/2)**k) for k in range(kmax)]
        meta={'witness_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'precision_bits':ctx.prec,
              'order':size,'A':str(A),'B':str(B),'T':str(T),'u_center_minus':str(h),'u_radius':str(radius),
              't_center_minus':str(T/2),'t_radius':str(T/2),'s_center':str(smid),'s_radius':str(sradius),'jmax':jmax,'kmax':kmax,
              'energy_degree_bound':de,'angular_degree_bound':da,
              'absolute_W_mass_bound':str(sum((abs(v) for j,n,v in rows),F(0))),
              'convention':'dx dy, W=(1-(2x-1)^2)*sum lambda P_n(2x-1)P_j(2y-1)'}
        results.append({'metadata':meta,'M_u':[[mu[j,k]*ar(T**j*radius**k) for k in range(kmax)] for j in range(jmax+1)],
                        'M_s':[[ms[j,k]*ar(T**j*sradius**k) for k in range(kmax)] for j in range(jmax+1)],
                        'T':[mu[j,0]*ar(T**j) for j in range(jmax+1)],'M_t':mt})
    return results

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('witness',type=Path);ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--jmax',type=int,default=40);ap.add_argument('--kmax',type=int,default=32)
    a=ap.parse_args();ctx.prec=240;tic=time.time()
    for d in build(a.witness,a.jmax,a.kmax):
        name='strip_'+d['metadata']['A']+'_'+d['metadata']['B']+'_moments.json'
        save(d,a.outdir/name);print(json.dumps(d['metadata'],indent=2));print('T0',d['T'][0])
    print('seconds',time.time()-tic)
