"""Exact contraction of source moment polynomials before interval evaluation.
This replaces ell+1 degreeN polynomial evaluations by one degree(N+ell)
polynomial. It changes neither moments nor their analytic remainder bounds.
"""
from finite_gap import *
class FiniteFast(FiniteGap):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        q=arb_poly([0,1]);ell=self.ell
        for s in self.sources:
            A=4+s['hu'];one=arb_poly([1,-A]);u=arb_poly([]);uv=arb_poly([])
            for j,c in enumerate(self.coeff):
                u+=exact(c)*q**j*one**(ell-j)*s['U'][j]
                uv+=exact(c)*q**j*arb_poly([1,-4])**(ell-j)*s['UVa'][j]
            s['combined_u']=q*u;s['combined_uv']=uv
    def full_norm(self,beta):
        e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C);rs=pe(self.Rstar,d);dl=d**self.ell
        U,V0=self.base.core(beta)
        vn=V0*e*e*dl*Q**3/rs-self.alpha*32*Q**2/(9*(d+exact(C))*mu)*(1-dl/rs)
        vn+=self.lamx*pe(self.crossD,d)/(Q*Q*rs)
        _,dv=self.base.derivative(beta,normalized=True);vn+=dv+self.base.higher(beta,normalized=True)
        err=arb(0);log_factors=[]
        for si,s in enumerate(self.sources):
            den=mu+s['hu'];qi=1/den;ud=pe(s['combined_u'],qi)/qi**self.ell
            td=pe(s['Mt'],1/(mu+s['ht']))/(mu+s['ht'])
            term=pe(s['Bpoly'],d)-pe(s['Ravg'],d)/Q+ud+rs*(td-2*s['mass']/Q)
            if s['kind']=='roy':term+=pe(s['Tphys'],d)
            pmax=sum(self.coeff[j]*val(d).abs_upper()**(self.ell-j)*s['Tmax']**j for j in range(self.ell+1))
            er=s['norm']*pmax*power(s['ru']/val(den).abs_lower(),self.N)/(val(den).abs_lower()-s['ru'])
            dt=val(mu+s['ht']);er+=val(rs).abs_upper()*s['norm']*power(s['rt']/dt.abs_lower(),self.N)/(dt.abs_lower()-s['rt'])
            vn-=Q**3/rs*term;err+=val(Q**3/rs).abs_upper()*er
            y=2*(mu-exact(s['A']))/exact(s['width'])-1
            for q,hp in zip(s['q'],s['H']):
                # Freeze the same PV approximant for midpoint and range jets.
                # The complete range chooses the branch; never re-choose at
                # the midpoint, where a different approximant could be used.
                roy=s['kind']=='roy'
                if si in self._omit_pv:
                    qq=-pe(q['R'],y);qe=arb(0)
                elif self._pv_modes[si]:
                    yy=val(y);inv=1/y;qq=pe(q['moments'],inv)*inv
                    qe=q['norm']*power(1/yy.abs_lower(),self.N)/(yy.abs_lower()-1)
                else:
                    da=y+1;db=y-1
                    if val(da)<0:da=-da
                    if val(db)<0:db=-db
                    qq=pe(q['W'],y)/2*(da.log()-db.log())-pe(q['R'],y);qe=arb(0)
                fac=4*dl*Q**3/(exact((2*self.ell+1)*s['width'])*rs) if roy else 2*Q**3*pe(hp,d)/(exact(s['width'])*rs)
                vn-=fac*qq;err+=val(fac).abs_upper()*qe
        return vn,err,log_factors

    def evaluate(self,a,b,series=False,omit_sources=()):
        if a==b==1 and not series:return arb(0),arb(0),arb(0),[]
        ball=arb(exact((a+b)/2),exact((b-a)/2));beta=arb_series([ball,1],prec=3) if series else ball
        e=1-beta*beta;x=e/4
        U,V=self.forward_series(beta) if series else self.base.forward(beta)
        du,dv=self.base.derivative(beta,u_beta=False);U+=du;V+=dv+self.base.higher(beta)
        ell=self.ell;logs=[];error=arb(0)
        if 4/(1-exact(a)*exact(a))>=exact(self.uv_start):
            h=x/(1-4*x);pstar=sum(self.coeff[j]*(exact(C)*h)**j for j in range(ell+1))
            z=1+2*exact(C)*h
            V-=2*self.alpha*x*(pstar-1)/(9*(1-exact(C)*x)*(1-exact(2*C)*x))
            V+=self.lamx*x**3/(16*(1-exact(C)*x)**5)*pe(self.crossZ,z)
            for s in self.sources:
                vs=pe(s['combined_uv'],x)/(1-4*x)**ell+pstar*pe(s['UVb'],x)
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
            den=mu+s['hu'];du=val(den);qi=1/den;ud=pe(s['combined_u'],qi)/qi**ell
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

