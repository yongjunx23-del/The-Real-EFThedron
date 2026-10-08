"""Fourth-order Taylor enclosures of normalized fields on finite-energy cells.

Source PV uses one exact polynomial/log expression for both point and interval
jets. Source geometric approximants have uniform value remainders; their
derivatives are not treated as derivatives of those error balls.
"""
from finite_global import *
import head_certification_lowspin as _low

def polynomial_series_all(poly,x):
    mid=x[0].mid();shift=poly(arb_poly([mid,arb(1)]));delta=x-mid
    result=x*0
    for v in reversed(shift.coeffs()):result=result*delta+v
    return result

# This affects only this Python process, never edits the frozen helper file.
# The old helper always initialized prec=3; preserve the caller's precision.
_low.series_poly=polynomial_series_all

class FiniteQuartic(FiniteGlobal):
    def normalized_interval(self,a,b):
        if self.ell==0:return FiniteCentered.normalized_interval(self,a,b)
        ma=F(4)/(1-a*a);mb=F(4)/(1-b*b);active=[]
        for i,s in enumerate(self.sources):
            if ma>=s['A'] and mb<=s['B']:active.append(i)
            elif not (mb<=s['A'] or ma>=s['B']):raise ValueError('Cell crosses a source endpoint')
        self._active_sources=active
        self.choose_pv_modes(a,b)
        mid=exact((a+b)/2);rad=exact((b-a)/2);beta=arb(mid,rad)
        bj=arb_series([beta,1],prec=5);mj=arb_series([mid,1],prec=5)
        vm,em,_=self.full_norm(mj);vb,eb,_=self.full_norm(bj)
        vr=sum((rad**k*vm[k].abs_upper() for k in range(1,4)),arb(0))+rad**4*vb[4].abs_upper()+eb
        vn=vm[0]+arb(0,vr.upper())
        um=self.real_field(mj);ub=self.real_field(bj)
        ur=sum((rad**k*um[k].abs_upper() for k in range(1,4)),arb(0))+rad**4*ub[4].abs_upper()
        U=um[0]+arb(0,ur.upper())
        e=1-beta*beta;mu=4/e;d=mu-4;Q=mu-exact(C);rs=pe(self.Rstar,d);pinv=d**self.ell*Q**3/rs
        return U,vn,pinv,e
