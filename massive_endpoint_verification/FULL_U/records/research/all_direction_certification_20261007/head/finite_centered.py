"""Taylor-centered normalized intervals suppress removable threshold cancellation.
The underlying functional is exactly unchanged. Only enclosure evaluation differs.
"""
from finite_taylor import *
class FiniteCentered(FiniteTaylor):
    def full_norm(self,beta):
        e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C);rs=pe(self.Rstar,d);dl=d**self.ell
        U,V0=self.base.core(beta)
        vn=V0*e*e*dl*Q**3/rs-self.alpha*32*Q**2/(9*(d+exact(C))*mu)*(1-dl/rs)
        vn+=self.lamx*pe(self.crossD,d)/(Q*Q*rs)
        _,dv=self.base.derivative(beta,normalized=True);vn+=dv+self.base.higher(beta,normalized=True)
        err=arb(0);log_factors=[]
        for s in self.sources:
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
                roy=s['kind']=='roy';qq,qe,_=self.qfun(q,y,omit=False)
                fac=4*dl*Q**3/(exact((2*self.ell+1)*s['width'])*rs) if roy else 2*Q**3*pe(hp,d)/(exact(s['width'])*rs)
                vn-=fac*qq;err+=val(fac).abs_upper()*qe
                if roy:
                    # Absolute logarithmic coefficient after removing beta^(2ell).
                    c=2*(4/e)**self.ell*Q**3*pe(q['W'],y)/(exact((2*self.ell+1)*s['width'])*rs)
                    rest=(4/(e*(exact(s['width'])-d))).log()
                    log_factors.append((val(c).abs_upper(),val(rest).abs_upper()))
        return vn,err,log_factors


    def normalized_interval(self,a,b):
        if self.ell==0:return super().normalized_interval(a,b)
        mid=exact((a+b)/2);rad=exact((b-a)/2);beta=arb(mid,rad)
        bj=arb_series([beta,1],prec=3);mj=arb_series([mid,1],prec=3)
        vm,em,_=self.full_norm(mj);vb,eb,_=self.full_norm(bj)
        vn=vm[0]+arb(0,(rad*vm[1].abs_upper()+rad*rad*vb[2].abs_upper()+eb).upper())
        e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C);rs=pe(self.Rstar,d);pinv=d**self.ell*Q**3/rs
        um=self.real_field(mj);ub=self.real_field(bj)
        U=um[0]+arb(0,(rad*um[1].abs_upper()+rad*rad*ub[2].abs_upper()).upper())
        return U,vn,pinv,e

    def real_field(self,beta):
        e=1-beta*beta;mu=4/e
        U=self.base.core(beta)[0]+self.base.derivative(beta,normalized=True,u_beta=False)[0]
        for s in self.sources:
            if s['kind']!='roy':continue
            y=2*(mu-exact(s['A']))/exact(s['width'])-1
            for q in s['q']:U+=2*arb.pi()/exact((2*self.ell+1)*s['width'])*pe(q['W'],y)/(e*e)
        return U
