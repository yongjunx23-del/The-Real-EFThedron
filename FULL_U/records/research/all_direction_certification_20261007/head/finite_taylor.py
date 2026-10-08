"""Threshold cap with a rigorously negative quadratic V/p coefficient.

V/p = A(beta)+sum c(beta) log(d/(width-d))+source error.  A is an
even analytic function at beta=0.  Arb Taylor coefficients enclose A''/2
on the whole cap; the explicit beta^(2ell) log beta terms are bounded
against beta^2. This avoids tiny constant-margin interval subdivisions.
"""
from finite_combined import *

class FiniteTaylor(FiniteCombined):
    def analytic_vnorm(self,beta):
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
                roy=s['kind']=='roy';qq,qe,_=self.qfun(q,y,omit=roy)
                fac=4*dl*Q**3/(exact((2*self.ell+1)*s['width'])*rs) if roy else 2*Q**3*pe(hp,d)/(exact(s['width'])*rs)
                vn-=fac*qq;err+=val(fac).abs_upper()*qe
                if roy:
                    # Absolute logarithmic coefficient after removing beta^(2ell).
                    c=2*(4/e)**self.ell*Q**3*pe(q['W'],y)/(exact((2*self.ell+1)*s['width'])*rs)
                    rest=(4/(e*(exact(s['width'])-d))).log()
                    log_factors.append((val(c).abs_upper(),val(rest).abs_upper()))
        return vn,err,log_factors

    def threshold_cap(self,b):
        if self.ell>=6:return FiniteSpin.threshold_cap(self,b)
        if self.ell<2:return super().threshold_cap(b)
        beta=arb(0).union(exact(b));jet=arb_series([beta,1],prec=3)
        av,err,logs=self.analytic_vnorm(jet)
        a0,e0,_=self.analytic_vnorm(arb(0))
        logco=sum((c*exact(b)**(2*self.ell-2)*(2*abs(exact(b).log())+rest) for c,rest in logs),arb(0))
        assert exact(b)<arb('0.01')
        # x^n |log x| increases on this interval: n>=2 and b<.01.
        constant=a0.upper()+err.upper()
        quadratic=av[2].upper()+logco.upper()
        if not (constant<0 and quadratic<0):raise ValueError('Quadratic threshold sign not resolved: '+str((constant,quadratic)))
        e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C);rs=pe(self.Rstar,d)
        U,_=self.base.core(beta);Ub=beta*U+self.base.derivative(beta,normalized=True)[0]
        for s in self.sources:
            if s['kind']!='roy':continue
            y=2*(mu-exact(s['A']))/exact(s['width'])-1
            for q in s['q']:Ub+=beta*2*arb.pi()/exact((2*self.ell+1)*s['width'])*pe(q['W'],y)/(e*e)
        coef=4**self.ell*Q.abs_upper()**3/(e.abs_lower()**self.ell*rs.abs_lower())
        bound=self.factor*Ub.abs_upper()**2*coef*exact(b)**(2*self.ell-3)/(2*(2*self.ell-3)*(-quadratic))
        self.cap_status=dict(method='Even analytic quadratic cap plus explicit beta^(2ell)log beta bound',beta_end=str(b),
            Vnorm_constant_upper=str(constant),Vnorm_beta2_upper=str(quadratic),beta_U_scaled=str(Ub),
            upper=str(bound.upper()),upper_dyadic=list(map(str,bound.upper().man_exp())))
        return bound.upper()

    def gaps(self,bits=30):
        gaps=SpinZero.gaps(self,bits)
        if self.ell:
            for exponent in range(8,31):
                end=F(1,2**exponent)
                try:self.threshold_cap(end);break
                except ValueError:continue
            else:raise ValueError('No finite-width quadratic threshold sign cap')
            gaps[0]=(F(0),end)
        return gaps

if __name__=='__main__':
    ctx.prec=256;path=ROOT/'output/full_direction_20261007/right/g2_guarded192_refined.npz';mom=ROOT/'output/all_direction_certification_20261007/head/g2/moments128'
    for ell in (2,4,6,14,40):
        m=FiniteTaylor(path,ell,mom,N=128);m.gaps();print(ell,m.cap_status,flush=True)
