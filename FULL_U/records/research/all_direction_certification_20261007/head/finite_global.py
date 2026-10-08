"""Centered normalized support bounds across all finite source windows.
Source activity is decided from exact rational beta endpoints. The same
functional, moment errors and source-PV formulas are used as finite220.py.
"""
from finite_centered import *
from head_certification_spin0 import support_lower
class FiniteGlobal(FiniteCentered):
    def full_norm(self,beta):
        e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C);rs=pe(self.Rstar,d);dl=d**self.ell
        U,V0=self.base.core(beta)
        vn=V0*e*e*dl*Q**3/rs-self.alpha*32*Q**2/(9*(d+exact(C))*mu)*(1-dl/rs)
        vn+=self.lamx*pe(self.crossD,d)/(Q*Q*rs)
        _,dv=self.base.derivative(beta,normalized=True);vn+=dv+self.base.higher(beta,normalized=True)
        err=arb(0);log_factors=[]
        for si,s in enumerate(self.sources):
            den=mu+s['hu'];ud=sum(self.coeff[j]*d**(self.ell-j)*pe(s['U'][j],1/den)/den for j in range(self.ell+1))
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
                if self._pv_modes[si]:
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



    def normalized_interval(self,a,b):
        ma=F(4)/(1-a*a);mb=F(4)/(1-b*b)
        active=[]
        for i,s in enumerate(self.sources):
            if ma>=s['A'] and mb<=s['B']:active.append(i)
            elif not (mb<=s['A'] or ma>=s['B']):raise ValueError('Cell crosses a source endpoint')
        self._active_sources=active
        self.choose_pv_modes(a,b)
        return super().normalized_interval(a,b)

    def choose_pv_modes(self,a,b):
        beta=arb(exact((a+b)/2),exact((b-a)/2));mu=4/(1-beta*beta)
        self._pv_modes=[(2*(mu-exact(s['A']))/exact(s['width'])-1).abs_lower()>arb('1.5') for s in self.sources]

    def real_field(self,beta):
        e=1-beta*beta;mu=4/e;d=mu-4
        U=self.base.core(beta)[0]+self.base.derivative(beta,normalized=True,u_beta=False)[0]
        for i in self._active_sources:
            s=self.sources[i];y=2*(mu-exact(s['A']))/exact(s['width'])-1
            for q,hp in zip(s['q'],s['H']):
                fac=4/(exact((2*self.ell+1)*s['width'])*e*e) if s['kind']=='roy' else 2*pe(hp,d)/(exact(s['width'])*d**self.ell*e*e)
                U+=arb.pi()/2*fac*pe(q['W'],y)
        return U

    def trapezoid(self,a,b):
        up,low=super().trapezoid(a,b)
        if self.ell and F(1,2)<=b<F(99,100):
            try:
                U,Vn,pinv,e=self.normalized_interval(a,b)
                if Vn.upper()<0:
                    candidate=self.factor*exact(b-a)*U.abs_upper()**2*pinv.abs_upper()*e.abs_upper()**2/(2*(-Vn.upper()))
                    if candidate.is_finite() and (not up.is_finite() or candidate<up):up=candidate.upper()
            except (ValueError,ZeroDivisionError):pass
        return up,low
