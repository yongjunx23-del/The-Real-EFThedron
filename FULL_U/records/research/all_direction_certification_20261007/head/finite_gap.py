"""Normalized finite source-endpoint caps for ell>6.
No Roy PV logarithm remains above the finite local Roy spins. Taper1 strips
vanish linearly at source endpoints, so their log products are bounded by
explicit |x log|x|| estimates. These are uniform bounds across the endpoint.
"""
from finite_quartic import *
class FiniteGap(FiniteQuartic):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self._omit_pv=set();self.gap_records=[]
        y=sy.Symbol('y')
        for s in self.sources:
            if s['kind']!='strip':continue
            byj={}
            for name,v in self.rows:
                mm=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
                if mm and F(mm[1])==s['A'] and F(mm[2])==s['B']:
                    byj.setdefault(int(mm[4]),[]).append((int(mm[5]),v))
            assert len(byj)==len(s['q'])
            for q,(j,rows) in zip(s['q'],sorted(byj.items())):
                q['untapered']=apoly(sum((sy.Rational(v)*sy.legendre(n,y) for n,v in rows),sy.Integer(0)),y)
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

    @staticmethod
    def xlog_abs(x):
        """Uniform upper for |x log|x||, with continuous value0 at0."""
        if x.contains(0):
            r=x.abs_upper();assert r<arb('0.1')
            return (-r*r.log()).upper() if r>0 else arb(0)
        return (x.abs_upper()*abs(abs(x).log()).upper()).upper()

    def gap_bound(self,a,b):
        if a==0 or self.ell<=6:return super().gap_bound(a,b)
        beta=arb(exact((a+b)/2),exact((b-a)/2));e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C)
        rs=pe(self.Rstar,d);pinv=d**self.ell*Q**3/rs
        self.choose_pv_modes(a,b)
        self._omit_pv={i for i,s in enumerate(self.sources) if any(not (r<exact(a) or r>exact(b)) for r in [(1-exact(4)/exact(s[k])).sqrt() for k in ('A','B')])}
        try:vn,err,_=self.full_norm(beta)
        finally:omitted=self._omit_pv;self._omit_pv=set()
        U=self.base.core(beta)[0]+self.base.derivative(beta,normalized=True,u_beta=False)[0]
        ma=F(4)/(1-a*a);mb=F(4)/(1-b*b);log_error=arb(0)
        for si,s in enumerate(self.sources):
            y=2*(mu-exact(s['A']))/exact(s['width'])-1
            for q,hp in zip(s['q'],s['H']):
                if s['kind']=='roy':
                    # Current inventory has J<=6, so this polynomial is
                    # identically zero for the present ell>6.
                    assert not q['W'].coeffs()
                    continue
                hpv=pe(hp,d)
                if not (mb<=s['A'] or ma>=s['B']):
                    u=arb.pi()*hpv*pe(q['W'],y)/(exact(s['width'])*d**self.ell*e*e)
                    if not (ma>=s['A'] and mb<=s['B']):u=u.union(arb(0))
                    U+=u
                if si in omitted:
                    ym=y-1;yp=y+1
                    tapered_log=yp.abs_upper()*self.xlog_abs(ym)+ym.abs_upper()*self.xlog_abs(yp)
                    # Vn log coefficient is -(Q^3 H/(width Rstar))*W.
                    fac=Q**3*hpv/(exact(s['width'])*rs)
                    log_error+=fac.abs_upper()*pe(q['untapered'],y).abs_upper()*tapered_log
        vn+=arb(0,(err+log_error).upper())
        if vn.upper()<0:
            bound=self.factor*exact(b-a)*U.abs_upper()**2*pinv.abs_upper()*e.abs_upper()**2/(2*(-vn.upper()))
            self.gap_records.append(dict(a=str(a),b=str(b),V_over_p=str(vn),upper=str(bound.upper()),log_error=str(log_error),method='Negative normalized V plus exact taper log-product bound'))
            return bound.upper()
        # A valid but possibly weaker fallback preserves completeness.
        bound=super().gap_bound(a,b)
        self.gap_records.append(dict(a=str(a),b=str(b),V_over_p=str(vn),upper=str(bound),method='Fallback absolute log-integral bound'))
        return bound



