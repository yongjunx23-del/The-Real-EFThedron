"""Arb enclosure of fixed-witness spin 2,4,6 disk-support integrals.

The right-channel PV and removable crossed pole are polynomial identities.
All remaining source integrals use certified exact polynomial moments with
explicit geometric remainders; UV uses a separately factored source series.
The outer integral uses convex interpolation and interval Riemann bounds.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re,sys,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sy
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import arb,arb_poly,arb_series,ctx
from certify_real_tail import exact,power
from certified_roy_moments import load
from certified_uv_moments import load_uv
from head_certification_spin0 import SpinZero,poly_eval,series_poly,support_upper,support_lower,integrate
C=F(4,3)

def apoly(expr,x):
    p=sy.Poly(sy.expand(expr),x)
    return arb_poly([exact(F(p.nth(i))) for i in range(max(0,p.degree())+1)])

def pe(poly,x):
    return series_poly(poly,x) if isinstance(x,arb_series) else poly_eval(poly,x)

def val(x):return x[0] if isinstance(x,arb_series) else x

def qdata(expr,y,N,norm):
    p=sy.Poly(sy.expand(expr),y);a=[F(p.nth(i)) for i in range(max(0,p.degree())+1)]
    r=[sum((a[m]*F(1,m-i) for m in range(i+1,len(a)) if (m-1-i)%2==0),F(0)) for i in range(max(1,len(a)-1))]
    moments=[sum((q/F(i+k+1) for i,q in enumerate(a) if (i+k)%2==0),F(0)) for k in range(N)]
    return dict(W=arb_poly(list(map(exact,a))),R=arb_poly(list(map(exact,r))),moments=arb_poly(list(map(exact,moments))),norm=exact(norm))

class LowSpin(SpinZero):
    def __init__(self,path,ell,N=64,uv_start=128):
        self.base=SpinZero(path,N,uv_start);self.path=self.base.path;self.sha=self.base.sha
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
            mompath=ROOT/f'output/certification_head/N64/{label}_moments.json'
            uvpath=ROOT/f'output/certification_head/{label}_uv_moments.json'
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
                            angular=sy.integrate(((1-z)/2)**j*sy.legendre(J,z),(z,0,1))
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
                    hp=sum(self.coeff[k]*d**(ell-k)*sy.integrate(sy.legendre(j,y)*(cl*(y-1)/2)**k,(y,-1,1))/2 for k in range(ell+1))
                    s['H'].append(apoly(hp,d))
            self.sources.append(s)

    def forward_series(self,beta):
        b2=beta*beta;e=1-b2;xmu=e/4;x=(1-6*b2+b2*b2)/(1+b2)**2
        rg=series_poly(self.base.G,x);kg=series_poly(self.base.H,x)
        for r,c,v in self.base.rats:
            y2=1-x*x;dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x
            nr=1-r*c*x;ni=-r*c;mminus=(1+r*c)/(1+2*r*c+r*r)
            den=(2*r*x-c*(1+r*r))**2+(1-c*c)*(1-r*r)**2
            mr=(nr*dr+y2*ni*di)/den;mi=(ni*dr-nr*di)/den
            rg+=v*((2*x*x-1)*(mr-mminus)-2*x*y2*mi)
            kg+=v*((2*x*x-1)*mi+2*x*(mr-mminus))
        return 16*beta*e*e/(1+b2)**5*kg,4*e/(1+b2)**3*rg-self.alpha*xmu/(2*(1-2*xmu)*(3-4*xmu))

    def qfun(self,q,y,omit=False):
        yy=val(y)
        if not omit and yy.abs_lower()>arb('1.5'):
            inv=1/y;out=pe(q['moments'],inv)*inv
            err=q['norm']*power(1/yy.abs_lower(),self.N)/(yy.abs_lower()-1)
            return out,err.upper(),None
        w=pe(q['W'],y);r=pe(q['R'],y)
        if omit:return -r,arb(0),w/2
        da=y+1;db=y-1
        if val(da)<0:da=-da
        if val(db)<0:db=-db
        return w/2*(da.log()-db.log())-r,arb(0),None

    def evaluate(self,a,b,series=False,omit_sources=()):
        if a==b==1 and not series:return arb(0),arb(0),arb(0),[]
        ball=arb(exact((a+b)/2),exact((b-a)/2));beta=arb_series([ball,1],prec=3) if series else ball
        e=1-beta*beta;x=e/4
        U,V=self.forward_series(beta) if series else self.base.forward(beta)
        ell=self.ell;logs=[];error=arb(0)
        if 4/(1-exact(a)*exact(a))>=exact(self.uv_start):
            h=x/(1-4*x);pstar=sum(self.coeff[j]*(exact(C)*h)**j for j in range(ell+1))
            z=1+2*exact(C)*h
            V-=2*self.alpha*x*(pstar-1)/(9*(1-exact(C)*x)*(1-exact(2*C)*x))
            V+=self.lamx*x**3/(16*(1-exact(C)*x)**5)*pe(self.crossZ,z)
            for s in self.sources:
                vs=sum(self.coeff[j]*h**j*pe(s['UVa'][j],x) for j in range(ell+1))+pstar*pe(s['UVb'],x)
                V-=x/16*vs
                xx=val(x);hh=val(h);pp=val(pstar);q=exact(s['B'])*xx
                pmax=sum(self.coeff[j]*(s['Tmax']*hh).abs_upper()**j for j in range(ell+1))
                er=xx.abs_upper()/16*4*s['norm']*exact(s['B']**2)*(pmax+pp.abs_upper())*power(q.abs_upper(),s['UVn'])/(1-q.abs_upper())
                error+=er.abs_upper()
            return U,V,error,logs
        mu=4/e;d=mu-4;Q=mu-exact(C);rs=pe(self.Rstar,d);dl=d**ell
        V-=self.alpha*32*(rs/dl-1)/(9*Q*(d+exact(C))*mu*e*e)
        V+=self.lamx*pe(self.crossD,d)/(dl*Q**5*e*e)
        for si,s in enumerate(self.sources):
            den=mu+s['hu'];du=val(den);ud=0
            for j in range(ell+1):ud+=self.coeff[j]*d**(ell-j)*pe(s['U'][j],1/den)/den
            td=pe(s['Mt'],1/(mu+s['ht']))/(mu+s['ht'])
            term=pe(s['Bpoly'],d)-pe(s['Ravg'],d)/Q+ud+rs*(td-2*s['mass']/Q)
            if s['kind']=='roy':term+=pe(s['Tphys'],d)
            V-=term/(dl*e*e)
            pmax=sum(self.coeff[j]*val(d).abs_upper()**(ell-j)*s['Tmax']**j for j in range(ell+1))
            er_u=s['norm']*pmax*power(s['ru']/du.abs_lower(),self.N)/(du.abs_lower()-s['ru'])
            dt=val(mu+s['ht'])
            er_t=s['norm']*power(s['rt']/dt.abs_lower(),self.N)/(dt.abs_lower()-s['rt'])
            error+=(er_u+val(rs).abs_upper()*er_t)/(val(dl)*val(e)*val(e)).abs_lower()
            y=2*(mu-exact(s['A']))/exact(s['width'])-1
            active=not(val(mu).upper()<=exact(s['A']) or val(mu).lower()>=exact(s['B']))
            interior=val(mu).lower()>exact(s['A']) and val(mu).upper()<exact(s['B'])
            if series:
                # The analytic branch is determined by exact beta endpoints,
                # not an inflated interval image of mu. Bulk cells exclude
                # actual source endpoints via the separate logarithmic gaps.
                ma=4/(1-a*a);mb=4/(1-b*b)
                active=ma>s['A'] and mb<s['B']
                if not active and not (mb<=s['A'] or ma>=s['B']):
                    raise ValueError('Bulk series cell crosses a source endpoint')
                interior=active
            for q,hp in zip(s['q'],s['H']):
                hpv=pe(hp,d)
                fac=4/(exact((2*ell+1)*s['width'])*e*e) if s['kind']=='roy' else 2*hpv/(exact(s['width'])*dl*e*e)
                qq,qe,logc=self.qfun(q,y,si in omit_sources)
                V-=fac*qq;error+=val(fac).abs_upper()*qe
                if logc is not None:logs.append((si,val(fac*logc)))
                if active:
                    ul=arb.pi()/2*fac*pe(q['W'],y)
                    if not interior and not series:ul=ul.union(arb(0))
                    U+=ul
        return U,V,error,logs

    def fields(self,a,b,omit_sources=()):
        U,V,err,logs=self.evaluate(a,b,False,omit_sources)
        return U,V+arb(0,err.upper()),logs

    def fields_series(self,a,b):
        U,V,err,_=self.evaluate(a,b,True)
        return U,V,err

    def trapezoid(self,a,b):
        # Both formulas are separately rigorous; use the sharper one.
        width=exact(b-a)
        u,v,_=self.fields(a,b)
        riemann=self.factor*width*support_upper(u,v)
        if b<F(1,2):
            nu,nv,pinv,e=self.normalized_interval(a,b)
            if nv.upper()<0:
                small=self.factor*width*nu.abs_upper()**2*pinv.abs_upper()*e.abs_upper()**2/(2*(-nv.upper()))
                if small.is_finite() and (not riemann.is_finite() or small<riemann):riemann=small
        try:
            U,V,err=self.fields_series(a,b)
            for q in (a,b):
                if q not in self.endpoint_cache:
                    qu,qv,_=self.fields(q,q);self.endpoint_cache[q]=support_upper(qu,qv)
            rem=width**3/8*(2*U[2].abs_upper()+4*V[2].abs_upper())+4*width*err
            trap=self.factor*(width/2*(self.endpoint_cache[a]+self.endpoint_cache[b])+rem)
            up=min((q for q in (riemann,trap) if q.is_finite()),key=float)
        except (ValueError,ZeroDivisionError):up=riemann
        mid=(a+b)/2;um,vm,_=self.fields(mid,mid)
        return up.upper(),(self.factor*width*support_lower(um,vm)).lower()

    def normalized_interval(self,a,b):
        beta=arb(exact((a+b)/2),exact((b-a)/2));e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C)
        rs=pe(self.Rstar,d);dl=d**self.ell;pinv=dl*Q**3/rs
        U,V0=self.base.forward(beta)
        vn=V0*e*e*pinv-self.alpha*32*Q**2/(9*(d+exact(C))*mu)*(1-dl/rs)
        vn+=self.lamx*pe(self.crossD,d)/(Q*Q*rs)
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
                qq,qe,_=self.qfun(q,y)
                fac=4*dl*Q**3/(exact((2*self.ell+1)*s['width'])*rs) if s['kind']=='roy' else 2*Q**3*pe(hp,d)/(exact(s['width'])*rs)
                vn-=fac*(qq+arb(0,qe.upper()))
                if s['kind']=='roy':U+=2*arb.pi()/exact((2*self.ell+1)*s['width'])*pe(q['W'],y)/(e*e)
        return U,vn,pinv,e

    def gaps(self,bits=30):
        gaps=super().gaps(bits)
        # Coarse analytic threshold caps avoid enormous removable intermediate
        # interval values. Their sign and support bounds are checked below.
        gaps[0]=(F(0),F(1,2**14) if self.ell==2 else F(1,2**9))
        return gaps

    def threshold_cap(self,b):
        """Negative V/p and elementary beta cap, including beta=0 exactly."""
        beta=arb(0).union(exact(b));e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C)
        rs=pe(self.Rstar,d);dl=power(d,self.ell);pinv=dl*Q**3/rs
        U,V0=self.base.forward(beta)
        vn=V0*e*e*pinv-self.alpha*32*Q**2/(9*(d+exact(C))*mu)*(1-dl/rs)
        vn+=self.lamx*pe(self.crossD,d)/(Q*Q*rs)
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
                    # d^ell log d is monotone in this very small cap.
                    assert eps<arb('0.01')
                    qbound=eps**self.ell*(ww/2*(abs(eps.log())+abs((exact(s['width'])-d).log()).upper())+rr)
                    vn+=arb(0,(4*Q**3/(exact((2*self.ell+1)*s['width'])*rs)*qbound).abs_upper())
                    U+=2*arb.pi()/exact((2*self.ell+1)*s['width'])*pe(q['W'],y)/(e*e)
                else:
                    qq,qe,_=self.qfun(q,y)
                    vn-=2*Q**3*pe(hp,d)/(exact(s['width'])*rs)*(qq+arb(0,qe.upper()))
        if not vn.upper()<0:raise ValueError('Threshold sign failed: '+str(vn))
        # H(U,V)<=U²/(-2V), scaled measure gives 4(2l+1)/pi.
        # pinv<=dmax^ell Qmax³/Rstar_min and e<=1.
        bound=self.factor*exact(b)*U.abs_upper()**2*pinv.abs_upper()/(2*(-vn.upper()))
        self.cap_status={'V_over_p':str(vn),'V_over_p_lower':str(vn.lower()),'V_over_p_upper':str(vn.upper()),
            'V_over_p_upper_dyadic':list(map(str,vn.upper().man_exp())),'beta_end':str(b),'upper':str(bound.upper())}
        return bound.upper()

    def gap_bound(self,a,b):
        if a==0:return self.threshold_cap(b)
        omit=[]
        for i,s in enumerate(self.sources):
            if any(not (r<exact(a) or r>exact(b)) for r in [(1-exact(4)/exact(s[k])).sqrt() for k in ('A','B')]):omit.append(i)
        U,V,logs=self.fields(a,b,omit);total=exact(b-a)*(U.abs_upper()+2*V.abs_upper())
        for si,coef in logs:
            s=self.sources[si];ra=(1-exact(4)/exact(s['A'])).sqrt();rb=(1-exact(4)/exact(s['B'])).sqrt()
            lm=exact(b-a)*abs((exact(s['A'])/exact(s['B'])).log()).upper()
            for r in (ra,rb,-ra,-rb):lm+=self.log_mass(a,b,r)
            total+=2*coef.abs_upper()*lm
        return (self.factor*total).upper()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,default=ROOT/'output/angles_resummed/angle_000p000.npz')
    ap.add_argument('--spin',type=int,choices=[2,4,6],default=2);ap.add_argument('--target',default='.015')
    ap.add_argument('--precision',type=int,default=224);ap.add_argument('--max-cells',type=int,default=20000)
    ap.add_argument('--gap-bits',type=int,default=30);ap.add_argument('--validate',action='store_true')
    ap.add_argument('--out',type=Path);a=ap.parse_args();ctx.prec=a.precision;t=time.time();m=LowSpin(a.witness,a.spin)
    print(json.dumps({'event':'built','seconds':time.time()-t,'spin':a.spin}),flush=True)
    if a.validate:
        from angle_audit import fields
        ss=np.array([4.0001,4.02,4.5,5.,7.,10.,11.8,11.99,12.1,18.,23.9,24.1,36.,47.9,48.1,80.,127.,128.,160.,1000.,1e6])
        st=np.load(a.witness);uu,vv,pp=fields(st['names'],st['multipliers'],st['direction'],ss,[a.spin],kernel_order=96,angular_order=96)
        for i,s in enumerate(ss):
            beta=F(float(math.sqrt(1-4/s)));u,v,_=m.fields(beta,beta);e=1-float(beta)**2
            print(s,float(u)-uu[0,i]/e**2,float(v)-vv[0,i]*pp[0,i]/e**2,float(v),flush=True)
        print('cap',m.threshold_cap(F(1,2**30)),m.cap_status,flush=True)
        return
    r=integrate(m,a.target,a.max_cells,a.gap_bits)
    r.update(spin=a.spin,witness_sha256=m.sha,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        arithmetic=f'Arb{a.precision}',threshold=m.cap_status,method='Exact polynomial PV division; rigorous polynomial source moments plus geometric remainders; convex support interpolation and interval Riemann bounds',
        dependency_sha256={q:hashlib.sha256((ROOT/'research'/q).read_bytes()).hexdigest() for q in
            ('head_certification_spin0.py','certify_real_tail.py','certified_roy_moments.py','certified_uv_moments.py','validated_quadrature.py')},
        source_tables=[{k:s[k] for k in ('moment_path','moment_sha256','uv_path','uv_sha256')} for s in m.sources],
        claim='Fixed-witness single-spin dual support upper bound; other spins separate')
    out=a.out or ROOT/f'output/certification_head/spin{a.spin}_component.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(r,indent=2))
    print(json.dumps({k:v for k,v in r.items() if k!='cells_data'},indent=2))

if __name__=='__main__':main()
