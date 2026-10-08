"""General 220-row finite-spin Arb head enclosure.
Constructor source-polynomial identities adapted from head_certification_lowspin.
This file does not alter or write to the earlier fixed-150 certificate paths.
"""
from forward220 import *
import argparse,hashlib,json,time
from functools import lru_cache
from head_certification_lowspin import LowSpin,C,qdata,val
from head_certification_spin0 import SpinZero,integrate
from certified_roy_moments import load
from certified_uv_moments import load_uv
from certify_real_tail import power
from flint import ctx
@lru_cache(None)
def angular_moment(j,J):
    z=sy.Symbol('z');p=sy.Poly(sy.expand(((1-z)/2)**j*sy.legendre(J,z)),z)
    return sum((p.nth(k)/sy.Integer(k+1) for k in range(p.degree()+1)),sy.Integer(0))
@lru_cache(None)
def strip_moment(k,j):
    y=sy.Symbol('y');p=sy.Poly(sy.expand(sy.legendre(j,y)*(sy.Rational(4,3)*(y-1)/2)**k),y)
    return sum((p.nth(m)/sy.Integer(m+1) for m in range(p.degree()+1) if m%2==0),sy.Integer(0))
class FiniteSpin(LowSpin):
    def __init__(self,path,ell,moment_dir,N=64,uv_start=128):
        self.path=Path(path);self.sha=hashlib.sha256(self.path.read_bytes()).hexdigest();self.base=Forward220(path,ell);self.cap_status={}
        self.N=N;self.ell=ell;self.uv_start=F(uv_start);self.alpha=self.base.alpha
        self.endpoint_cache={};self.factor=4*(2*ell+1)/arb.pi()
        st=np.load(path);self.rows=[(str(n),F(float(v))) for n,v in zip(st['names'],st['multipliers'])]
        self.lamx=exact(dict(self.rows)['X4']);self.coeff=[math.comb(ell,j)*math.comb(ell+j,j) for j in range(ell+1)]
        d,t,y,z=sy.symbols('d t y z');cl=sy.Rational(C)
        rp=sum(self.coeff[j]*d**(ell-j)*t**j for j in range(ell+1));rs=rp.subs(t,cl)
        brem=sy.div(rs-rp,d+cl+t,t)
        assert sy.expand(brem[1])==0
        self.Rstar=apoly(rs,d)
        lp=sy.Poly(sy.legendre(ell,z),z)
        self.crossZ=apoly(sum(lp.nth(m)*m*(2*m-5)*z**m for m in range(ell+1)),z)
        self.crossD=apoly(sum(lp.nth(m)*m*(2*m-5)*(d+4-cl)**m*d**(ell-m) for m in range(ell+1)),d)
        self.sources=[]
        for label in ('roy','strip_12_24','strip_24_48'):
            mompath=Path(moment_dir)/f'{label}_moments.json'
            uvpath=Path(moment_dir)/f'{label}_uv_moments.json'
            m=load(mompath);uv=load_uv(uvpath);meta=m['metadata'];um=uv['metadata']
            assert meta['witness_sha256']==self.sha==um['witness_sha256']
            roy=label=='roy';A=F(4) if roy else F(meta['A']);B=F(meta['B']);width=B-A
            s=dict(kind='roy' if roy else 'strip',A=A,B=B,width=width,T=m['T'],mass=m['T'][0],
                norm=exact(F(meta['absolute_W_mass_bound'])),Tmax=exact(F(um['Tmax'])),
                hu=exact(4 if roy else F(meta['u_center_minus'])),ru=exact(4 if roy else F(meta['u_radius'])),
                ht=exact(2 if roy else F(meta['t_center_minus'])),rt=exact(2 if roy else F(meta['t_radius'])),
                U=[arb_poly(row[:N]) for row in m['M_u'][:ell+1]],Mt=arb_poly(m['M_t'][:N]),
                UVa=[arb_poly(row) for row in uv['A'][:ell+1]],UVb=arb_poly(uv['B']),UVn=int(um['N']),
                moment_path=str(mompath),moment_sha256=hashlib.sha256(mompath.read_bytes()).hexdigest(),
                uv_path=str(uvpath),uv_sha256=hashlib.sha256(uvpath.read_bytes()).hexdigest())
            rb=sy.Poly(brem[0],d,t)
            bc=[arb(0) for _ in range(ell+1)]
            for (di,ti),q in rb.terms():bc[di]+=exact(F(q))*m['T'][ti]
            s['Bpoly']=arb_poly(bc)
            s['Ravg']=arb_poly([exact(self.coeff[ell-i])*m['T'][ell-i] for i in range(ell+1)])
            if roy:
                selected=[];allrows=[]
                for name,v in self.rows:
                    mm=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
                    if mm:
                        J,n=int(mm[2]),int(mm[3]);allrows.append((J,n,v))
                        if J==ell:selected.append((n,v))
                expr=sum((sy.Rational(v)*sy.legendre(n,y) for n,v in selected),sy.Integer(0))
                s['q']=[qdata(expr,y,N,sum((abs(v) for n,v in selected),F(0)))];s['H']=[arb_poly([1])]
                tp=sy.Integer(0)
                for j in range(1,ell+1):
                    for a in range(j):
                        moment=F(0)
                        for J,n,v in allrows:
                            if a<n:continue
                            energy=F(math.factorial(a)**2,math.factorial(a-n)*math.factorial(a+n+1))
                            angular=angular_moment(j,J)
                            moment+=2*v*width**a*energy*F(angular)
                        tp+=(-1)**(j+1)*self.coeff[j]*sy.Rational(moment)*d**(ell-1-a)
                s['Tphys']=apoly(tp,d)
            else:
                byj={}
                for name,v in self.rows:
                    mm=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
                    if mm and F(mm[1])==A and F(mm[2])==B:
                        j,n=int(mm[4]),int(mm[5]);byj.setdefault(j,[]).append((n,v))
                s['q']=[];s['H']=[]
                for j,rows in sorted(byj.items()):
                    expr=(1-y*y)*sum((sy.Rational(v)*sy.legendre(n,y) for n,v in rows),sy.Integer(0))
                    s['q'].append(qdata(expr,y,N,sum((abs(v) for n,v in rows),F(0))))
                    hp=sum(self.coeff[k]*d**(ell-k)*strip_moment(k,j) for k in range(ell+1))
                    s['H'].append(apoly(hp,d))
            self.sources.append(s)

    def forward_series(self,beta):return self.base.core(beta)

    def evaluate(self,a,b,series=False,omit_sources=()):
        U,V,error,logs=super().evaluate(a,b,series,omit_sources)
        if a==b==1 and not series:return U,V,error,logs
        ball=arb(exact((a+b)/2),exact((b-a)/2));beta=arb_series([ball,1],prec=3) if series else ball
        du,dv=self.base.derivative(beta,u_beta=False)
        return U+du,V+dv+self.base.higher(beta),error,logs

    def normalized_interval(self,a,b):
        U,Vn,pinv,e=super().normalized_interval(a,b)
        beta=arb(exact((a+b)/2),exact((b-a)/2))
        du,dv=self.base.derivative(beta,normalized=True,u_beta=False)
        return U+du,Vn+dv+self.base.higher(beta,normalized=True),pinv,e

    def threshold_cap(self,b):
        # Copy the old exact source cap algebra, but bound beta*U rather
        # than U.  Derivative rows have a true beta^-1 factor at ell>=2.
        assert self.ell>=2
        beta=arb(0).union(exact(b));e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C)
        rs=pe(self.Rstar,d);dl=power(d,self.ell)
        U,V0=self.base.core(beta);Ub=beta*U
        vn=V0*e*e*dl*Q**3/rs-self.alpha*32*Q**2/(9*(d+exact(C))*mu)*(1-dl/rs)
        vn+=self.lamx*pe(self.crossD,d)/(Q*Q*rs)
        du,dv=self.base.derivative(beta,normalized=True)
        Ub+=du;vn+=dv+self.base.higher(beta,normalized=True)
        for s in self.sources:
            den=mu+s['hu'];ud=sum(self.coeff[j]*d**(self.ell-j)*pe(s['U'][j],1/den)/den for j in range(self.ell+1))
            td=pe(s['Mt'],1/(mu+s['ht']))/(mu+s['ht'])
            term=pe(s['Bpoly'],d)-pe(s['Ravg'],d)/Q+ud+rs*(td-2*s['mass']/Q)
            if s['kind']=='roy':term+=pe(s['Tphys'],d)
            pmax=sum(self.coeff[j]*d.abs_upper()**(self.ell-j)*s['Tmax']**j for j in range(self.ell+1))
            er=s['norm']*pmax*power(s['ru']/den.abs_lower(),self.N)/(den.abs_lower()-s['ru'])
            dt=mu+s['ht'];er+=rs.abs_upper()*s['norm']*power(s['rt']/dt.abs_lower(),self.N)/(dt.abs_lower()-s['rt'])
            vn-=Q**3/rs*(term+arb(0,er.upper()))
            y=2*(mu-exact(s['A']))/exact(s['width'])-1
            for q,hp in zip(s['q'],s['H']):
                if s['kind']=='roy':
                    eps=d.abs_upper();ww=pe(q['W'],y).abs_upper();rr=pe(q['R'],y).abs_upper()
                    assert eps<arb('0.01')
                    qbound=eps**self.ell*(ww/2*(abs(eps.log())+abs((exact(s['width'])-d).log()).upper())+rr)
                    vn+=arb(0,(4*Q**3/(exact((2*self.ell+1)*s['width'])*rs)*qbound).abs_upper())
                    Ub+=beta*2*arb.pi()/exact((2*self.ell+1)*s['width'])*pe(q['W'],y)/(e*e)
                else:
                    qq,qe,_=self.qfun(q,y)
                    vn-=2*Q**3*pe(hp,d)/(exact(s['width'])*rs)*(qq+arb(0,qe.upper()))
        if not vn.upper()<0:raise ValueError('Threshold sign failed: '+str(vn))
        # pinv/beta^2 = 4^ell beta^(2ell-2) Q^3/(e^ell Rstar).
        # Integrate the explicit beta power rather than bounding a 0/0.
        coef=4**self.ell*Q.abs_upper()**3/(e.abs_lower()**self.ell*rs.abs_lower())
        bound=self.factor*Ub.abs_upper()**2*coef*exact(b)**(2*self.ell-1)/(2*(2*self.ell-1)*(-vn.upper()))
        self.cap_status={'V_over_p':str(vn),'beta_end':str(b),'beta_U_scaled':str(Ub),
            'upper':str(bound.upper()),'upper_dyadic':list(map(str,bound.upper().man_exp()))}
        return bound.upper()

    def gaps(self,bits=30):
        gaps=SpinZero.gaps(self,bits)
        if self.ell:
            # Cap width is a proof parameter.  Find a passing sign enclosure.
            for exponent in (14,18,22,26,30,34):
                end=F(1,2**exponent)
                try:self.threshold_cap(end);break
                except ValueError:continue
            else:raise ValueError('No finite-width threshold sign cap')
            gaps[0]=(F(0),end)
        return gaps

    def gap_bound(self,a,b):
        if a==0 and self.ell==0:return SpinZero.gap_bound(self,a,b)
        return LowSpin.gap_bound(self,a,b)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--moments',type=Path,required=True)
    ap.add_argument('--spin',type=int,required=True);ap.add_argument('--target',default='2');ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--bits',type=int,default=256);ap.add_argument('--max-cells',type=int,default=5000);ap.add_argument('--validate',action='store_true')
    a=ap.parse_args();ctx.prec=a.bits;tic=time.time();model=FiniteSpin(a.witness,a.spin,a.moments)
    print(json.dumps(dict(event='built',spin=a.spin,seconds=time.time()-tic)),flush=True)
    if a.validate:
        from angle_audit import fields
        energies=np.array([4.0001,4.02,4.5,5.,7.,10.,11.8,11.99,12.1,18.,23.9,24.1,36.,47.9,48.1,80.,127.,128.,160.,1000.,1e6])
        st=np.load(a.witness);uu,vv,pp=fields(st['names'],st['multipliers'],st['direction'],energies,[a.spin],kernel_order=96,angular_order=96)
        rec=[]
        for i,s in enumerate(energies):
            beta=F(float(math.sqrt(1-4/s)));u,v,_=model.fields(beta,beta);e=1-float(beta)**2
            un=float(u)*e**2;vn=float(v)*e**2/pp[0,i]
            rec.append(dict(mu=float(s),U_relative=abs(un-uu[0,i])/(1+abs(uu[0,i])),Vnorm_relative=abs(vn-vv[0,i])/(1+abs(vv[0,i]))))
        print(json.dumps(rec),flush=True)
        print('GAPS',model.gaps(),flush=True)
        result=dict(spin=a.spin,witness_sha256=model.sha,comparison=rec,cap=model.cap_status,scope='Floating-point comparison diagnostic, cap enclosure separate')
    else:
        result=integrate(model,a.target,a.max_cells,30)
        result.update(spin=a.spin,witness_sha256=model.sha,threshold=model.cap_status,claim='Rigorous fixed-witness finite single-spin all-energy upper bound only; remaining spins and tails separate')
    result['elapsed_seconds']=time.time()-tic
    result['dependencies']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('forward220.py'),ROOT/'research/head_certification_lowspin.py',ROOT/'research/head_certification_spin0.py']}
    result['source_manifest_sha256']=hashlib.sha256((a.moments/'manifest.json').read_bytes()).hexdigest()
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='cells_data'},indent=2),flush=True)
if __name__=='__main__':main()

