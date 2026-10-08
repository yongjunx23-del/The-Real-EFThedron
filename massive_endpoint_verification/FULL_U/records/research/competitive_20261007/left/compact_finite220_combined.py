"""Finite compact-energy V/p for the frozen left180 witness.

Reuse the reviewed removable-pole/source geometric evaluator, with new
source moments and exact witness-specific norms. Construct the polynomial
divided difference by long division to avoid cubic coefficient bookkeeping.
No old certificate inputs are modified.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,math,hashlib,pickle,copy,re
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'research/tail_certification'),str(ROOT/'research/all_direction_certification_20261007/head'),str(ROOT/'research')]
import certified_finite_v as old
from certified_finite_v import arb,arb_poly,ctx,ar,coeff,rev,tight,geometric_num,taper_q,load
from forward_combined import ForwardCombined
from head_certification_lowspin import pe
X=arb_poly([0,1]);C=F(4,3)
WITNESS=ROOT/'output/full_direction_20261007/left/angle180_final.npz'
EXPECTED='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'

class Source(old.Source):
    def __init__(self,mom,spins,norm,roy=False):
        assert mom['metadata']['witness_sha256']==EXPECTED
        self.m=mom;self.roy=roy;self.uv=False;self.N=int(mom['metadata']['kmax']);N=self.N
        assert N>=3 and mom['metadata']['jmax']>=max(spins)
        md=mom['metadata']
        if roy:
            assert F(md['B'])==F(119,10) and md['Jmax']==6
            assert md['skmax']>=max(spins)-1
            self.h=F(4);self.rt=F(2);self.ru=F(4);self.tm=F(2)
        else:
            assert F(md['T'])==C
            self.h=F(md['u_center_minus']);self.ru=F(md['u_radius']);self.tm=F(md['t_center_minus']);self.rt=F(md['t_radius'])
        self.norm=F(norm);assert self.norm>0
        self.u_den=X+ar(4+self.h);self.t_den=X+ar(4+self.tm);self.q=X+ar(2*C)
        Du=self.u_den**N;Dt=self.t_den**N;self.polys={}
        for ell in spins:
            cs=[coeff(ell,j)for j in range(ell+1)]
            rs=arb_poly([ar(cs[ell-i]*C**(ell-i))for i in range(ell+1)])
            ir=arb_poly([cs[ell-i]*mom['T'][ell-i]for i in range(ell+1)])
            pus=[arb_poly([cs[ell-i]*mom['M_u'][ell-i][k]for i in range(ell+1)])for k in range(N)]
            Nu=geometric_num(pus,self.u_den,N)
            Nt=geometric_num([arb_poly([m])for m in mom['M_t']],self.t_den,N)
            # Quotient of (Rstar-R(d,t))/(t+d+c); even ell makes remainder0.
            qr=arb_poly([-cs[ell]]);dd=mom['T'][ell-1]*qr
            for r in range(ell-1,0,-1):
                qr=-cs[r]*(X**(ell-r))-(X+ar(C))*qr
                dd+=mom['T'][r-1]*qr
            pv=arb_poly([])
            if roy:
                pv=arb_poly([sum((cs[j]*mom['M_s'][j][ell-j-1-p]for j in range(ell-p)),arb(0))for p in range(ell)])
            base=(pv+dd)*self.q*Du*Dt-(ir+2*mom['T'][0]*rs)*Du*Dt+self.q*Nu*Dt+self.q*rs*Nt*Du
            if roy:
                target=ell+2*N-2
                assert all(base[k].contains(0)for k in range(target+1,base.degree()+1))
                base=arb_poly([base[k]for k in range(target+1)])
            self.polys[ell]=(base,rs,None)
    def eval(self,mu,ell):
        # The geometric remainder uses a scalar lower bound for Rstar;
        # a wide-cell enclosure crossing zero must be rejected, not divided.
        rs=tight(self.polys[ell][1],mu-4)
        if not rs.lower()>0:raise ValueError('Rstar interval is not positive')
        return super().eval(mu,ell)

class CompactFinite220:
    def __init__(self,moment_dir,spins=range(42,100,2)):
        self.spins=list(spins);self.path=WITNESS
        assert hashlib.sha256(WITNESS.read_bytes()).hexdigest()==EXPECTED
        self.forward=ForwardCombined(WITNESS,self.spins[0]);self.f={}
        assert self.forward.direction==[F(-1),F(0)]
        for ell in self.spins:
            f=copy.copy(self.forward);f.ell=ell;f.coeff=[coeff(ell,j)for j in range(ell+1)];self.f[ell]=f
        cd=pickle.loads((ROOT/'output/all_direction_certification_20261007/tail/left/refined/compact_norms.pkl').read_bytes())
        assert cd['witness_sha256']==EXPECTED
        self.roy=Source(load(Path(moment_dir)/'roy_moments.json'),self.spins,cd['roy_abs'],True)
        self.strips=[];self.hp={}
        for A,B in [(12,24),(24,48)]:
            src=Source(load(Path(moment_dir)/f'strip_{A}_{B}_moments.json'),self.spins,cd['strip_abs'])
            self.strips.append((A,B,src))
        self.groups={}
        for name,value in self.forward.rows.items():
            mm=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
            if mm:self.groups.setdefault((int(mm[1]),int(mm[2])),[]).append((int(mm[4]),int(mm[5]),ar(value)))
        self.jmax=max(j for group in self.groups.values() for j,n,v in group)
        for ell in self.spins:
            for j in range(self.jmax+1):
                self.hp[ell,j]=arb_poly([ar(F((-1)**(k+j)*math.factorial(ell+k),math.factorial(ell-k)*math.factorial(k-j)*math.factorial(k+j+1))*C**k)if k>=j else 0 for k in reversed(range(ell+1))])
        self.cross={}
        for ell in self.spins:
            cs=[coeff(ell,j)for j in range(ell+1)]
            rc=arb_poly([ar(cs[ell-i]*(ell-i)*C**(ell-i-1))if ell-i>=1 else 0 for i in range(ell+1)])
            rcc=arb_poly([ar(cs[ell-i]*(ell-i)*(ell-i-1)*C**(ell-i-2))if ell-i>=2 else 0 for i in range(ell+1)])
            self.cross[ell]=((X+ar(2*C))*rcc-3*rc)*ar(F(1,2))
        self.lx=ar(self.forward.rows['X4'])
        self.cross_all={};qp=X+ar(2*C)
        for ell in self.spins:
            rs=self.roy.polys[ell][1];cs=self.f[ell].coeff
            num=self.lx*self.cross[ell]*qp**5+ar(self.forward.direction[1])*rs*qp**6
            for row,lam in self.forward.cross:
                for j,bj in enumerate(row['beta']):
                    if not bj or j>ell:continue
                    q=row['degree']-2-j;assert 0<=q<=6
                    pj=arb_poly([ar(F(cs[ell-i]*math.comb(ell-i,j))*C**(ell-i-j)) if ell-i>=j else 0 for i in range(ell+1)])
                    num+=lam*ar(F(bj))*pj*qp**(6-q)
            self.cross_all[ell]=num
    def fields(self,mu,spins=None):
        # The caller's exact interval is inside[4,49]. Arb's radius
        # representation may enlarge an endpoint slightly; polynomial
        # source factors remain analytic in this harmless enlargement.
        assert mu.lower()>ar(F(399,100)) and mu.upper()<50
        delta=mu-4;Q=mu-ar(C);e=4/mu
        ml,mh=arb(mu.lower()),arb(mu.upper())
        blo=arb(0) if ml<=4 else (1-4/ml).sqrt().lower()
        bhi=(1-4/mh).sqrt().upper();beta=arb(blo).union(arb(bhi))
        strip_pv=[]
        for A,B,src in self.strips:
            zz=2*(mu-A)/(B-A)-1;vals=[taper_q(n,zz)for n in range(max(n for j,n,v in self.groups[A,B])+1)]
            strip_pv.append({j:sum((v*vals[n]*2/(B-A)for jj,n,v in self.groups[A,B]if jj==j),arb(0))for j in range(self.jmax+1)})
        out={}
        for ell in (self.spins if spins is None else spins):
            f=self.f[ell];v,rs=self.roy.eval(mu,ell);_,vf=f.core(beta)
            # core returns V/e^2 including the symmetric subtraction.
            v+=vf*e*e*delta**ell*Q**3/rs-f.alpha*32*Q**2/(9*(delta+ar(C))*mu)*(1-delta**ell/rs)
            v+=tight(self.cross_all[ell],delta)/(Q**6*rs)
            v+=f.derivative(beta,normalized=True)[1]
            for k,(A,B,src)in enumerate(self.strips):
                vs,_=src.eval(mu,ell)
                pv=sum((strip_pv[k][j]*tight(self.hp[ell,j],delta)for j in range(self.jmax+1)),arb(0))
                v+=vs-Q**3/rs*pv
            out[ell]=v
        return out
