"""Algebraically combined derivative evaluation for tighter interval widths.
Uses the same exact rational functional as forward220.py. No rows are removed.
"""
from forward220 import *
class ForwardCombined(Forward220):
    def derivative(self,beta,normalized=False,u_beta=True):
        """Return (beta*U/e^2, V/e^2 or V/p) for derivative rows.

        For ell>=2 the inner coordinate is t=beta^2*tau. Factoring
        (1-omega)^(2p-2) cancels all but beta^-2 before coefficient
        extraction.  For ell=0 the original t jets have no threshold pole.
        """
        z=beta*beta;e=1-z;zero=z*0;one=zero+1
        # All jets have four exact Taylor coefficients, through order three.
        scaled=self.ell>0;ts=z if scaled else one
        A=[one/2,ts/8,zero,zero];nr=[e*q for q in A]
        q0=(1+z)/2;q1=e*ts/8
        invq=[(-q1)**k/q0**(k+1) for k in range(4)]
        invq2=jp(invq,2)
        S=[exact(F(math.prod([F(1,2)-j for j in range(k)]),math.factorial(k)))*(e*ts/4)**k for k in range(4)]
        # omega=(nr+i beta*S)^2/q^2.
        nr2=jm(nr,nr);s2=jm(S,S)
        om=(jm([a-z*b for a,b in zip(nr2,s2)],invq2),jm([2*x for x in jm(nr,S)],invq2))
        # D=(1-omega); D^2/beta^2 has a regular pair representation.
        dr=jm([2*x for x in s2],invq2);di=jm([-2*x for x in jm(nr,S)],invq2)
        d2=([z*a-b for a,b in zip(jm(dr,dr),jm(di,di))],[2*x for x in jm(dr,di)])
        ff=[2*e*x for x in jm(jp(A,2),jp(invq,3))]
        pt=[exact(self.coeff[k])*(e/4)**k if k<=self.ell else zero for k in range(4)]
        if not scaled:pt=[one,zero,zero,zero]
        base=(jm(ff,pt),[zero]*4)
        # Exact Taylor coefficients about omega=1 combine cancellations
        # before interval evaluation. omega-1 uses its factored expression.
        omm=([-z*x for x in dr],[-x for x in di]);poly={}
        for p in (1,2,3):
            exact_rows=[(int(n.split('_j')[1]) if '_j' in n else 1,v) for n,v in self.rows.items() if n.startswith(f'Dt{p}')]
            coeff=[sum((v*math.comb(j,k) for j,v in exact_rows if j>=k),F(0)) for k in range(13)]
            pp=([zero]*4,[zero]*4)
            for c in reversed(coeff):
                pp=cm(pp,omm,z);pp[0][0]+=exact(c)
            poly[p]=pp
        d2pow={p:cp(d2,p-1,z) for p,_,_,_ in self.deriv}
        U=zero;V=zero;dbar=4*z;Qbar=4-exact(F(4,3))*e
        R=sum((exact(c)*dbar**(self.ell-k)*(exact(F(4,3))*e)**k for k,c in enumerate(self.coeff)),zero)
        ka=[rat(q,z) for q in self.ka];kb=[rat(q,z) for q in self.kb]
        for p in (1,2,3):
            regular=True;lam=exact(self.rows.get(f'Dt{p}G_regular',F(0)))
            vv=cm(cm(base,poly[p],z),d2pow[p],z)
            # For ell=0 there was no scaled coordinate; undo division by
            # beta^(2p-2) introduced by the common D^2 representation.
            factor=math.factorial(p)*(one if scaled else z**(p-1))
            U+=factor*vv[1][p]*((one if u_beta else 1/beta) if scaled else (z if u_beta else beta))
            if normalized:
                if scaled:
                    ratio=4**self.ell*z**(self.ell-1)*Qbar**3/(e*R)
                    # ff carries one e, so beta=1 requires a factored field;
                    # normalized evaluation is used only away from UV infinity.
                    V+=factor*vv[0][p]*ratio
                else:V+=factor*vv[0][p]*Qbar**3/e
            else:V+=factor*vv[0][p]/(z if scaled else one)
            if regular:
                if normalized:
                    if p==1:
                        term=sum((exact(self.coeff[k])*dbar**(self.ell-k)*rat(self.ka_div[p-k],z) for k in range(min(p,self.ell)+1)),zero)/R+rat(self.kb1_div,z)
                        sub=-math.factorial(p)*Qbar**3*term
                    else:
                        term=sum((exact(self.coeff[k])*dbar**(self.ell-k)*ka[p-k] for k in range(min(p,self.ell)+1)),zero)/R+kb[p]
                        sub=-math.factorial(p)*e**(p-2)*Qbar**3*term
                else:
                    term=sum((exact(self.coeff[k])*dbar**(self.ell-k)*ka[p-k] for k in range(min(p,self.ell)+1)),zero)+R*kb[p]
                    sub=-math.factorial(p)*e**(p-1)*term/dbar**self.ell
                V+=lam*sub
        return U,V

