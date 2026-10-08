"""Exact beta-factored Arb forward rows for the 220-row fixed witnesses.

All binary64 multipliers are interpreted as exact dyadic rationals.  Inner
Taylor coefficients use algebra, never finite differences.  A complex value
is stored as (real, imaginary/beta); multiplying such pairs is exact with
i^2=-1 and the explicitly supplied beta^2.
"""
from __future__ import annotations
import math,re,sys
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
import numpy as np
import sympy as sy
from flint import arb,arb_poly,arb_series
from certify_real_tail import exact
from head_certification_lowspin import pe,apoly
from crossing_kernels import crossing_catalogue

def jm(a,b):
    return [sum((a[k]*b[j-k] for k in range(j+1)),a[0]*0) for j in range(len(a))]
def jp(a,n):
    out=[a[0]*0+1]+[a[0]*0]*(len(a)-1)
    while n:
        if n%2:out=jm(out,a)
        n//=2
        if n:a=jm(a,a)
    return out
def cm(a,b,z):
    rr=jm(a[0],b[0]);ii=jm(a[1],b[1]);ri=jm(a[0],b[1]);ir=jm(a[1],b[0])
    return ([r-z*i for r,i in zip(rr,ii)],[r+i for r,i in zip(ri,ir)])
def cp(a,n,z):
    zero=a[0][0]*0;out=([zero+1]+[zero]*3,[zero]*4)
    while n:
        if n%2:out=cm(out,a,z)
        n//=2
        if n:a=cm(a,a,z)
    return out
def scale(a,x):return ([q*x for q in a[0]],[q*x for q in a[1]])
def rational(expr,z):
    num,den=sy.fraction(sy.cancel(expr));return apoly(num,z),apoly(den,z)
def rat(pair,x):return pe(pair[0],x)/pe(pair[1],x)

class Forward220:
    def __init__(self,path,ell):
        st=np.load(path);self.rows={str(n):F(float(v)) for n,v in zip(st['names'],st['multipliers'])}
        self.direction=list(map(lambda v:F(float(v)),st['direction']));self.ell=ell
        self.coeff=[math.comb(ell,j)*math.comb(ell+j,j) for j in range(ell+1)]
        alpha=self.direction[0]
        for n,v in self.rows.items():
            m=re.fullmatch(r'Roy_B[^_]+_J(\d+)_n(\d+)',n)
            if m and int(m[1])==int(m[2])==0:alpha-=2*v
            m=re.fullmatch(r'Taper1_Strip_A[^_]+_B[^_]+_T[^_]+_j(\d+)_n(\d+)',n)
            if m and int(m[1])==0:
                if int(m[2])==0:alpha-=F(2,3)*v
                if int(m[2])==2:alpha+=F(2,15)*v
        self.alpha=exact(alpha);self.alpha_exact=alpha
        x=sy.Symbol('x');g=sy.Rational(alpha)*x;h=sy.Rational(alpha)
        self.rats=[];self.deriv=[]
        for n,v in self.rows.items():
            m=re.fullmatch(r'H(\d+)',n)
            if m:
                j=int(m[1]);g+=sy.Rational(v)*sy.chebyshevt(j,x);h+=sy.Rational(v)*sy.chebyshevu(j-1,x)
            m=re.fullmatch(r'rational_(cos|sin)_r([^_]+)_c([^_]+)',n)
            if m:self.rats.append((exact(F(m[2])),exact(F(m[3])),m[1],exact(v)))
            m=re.fullmatch(r'Dt([123])(?:H_j(\d+)|G_regular)',n)
            if m:self.deriv.append((int(m[1]),int(m[2]) if m[2] else 1,not bool(m[2]),exact(v)))
        self.G=apoly(g,x);self.H=apoly(h,x)
        self.cross=[(q,exact(self.rows[q['name']])) for q in crossing_catalogue(8) if q['name'] in self.rows and q['name']!='X4']
        z=sy.Symbol('z');e=1-z;qa=4-2*e;qb=4-sy.Rational(8,3)*e;Q=4-sy.Rational(4,3)*e
        self.ka=[];self.kb=[];self.ka_div=[]
        for k in range(4):
            ka=(-1)**k*(2/(2**k*qa**(k+1))-1/qb**(k+1))-(1/Q if k==0 else 0)
            kb=1/sy.Integer(4)**(k+1)+(-1)**k/qb**(k+1)-(2/Q if k==0 else 0)
            self.ka.append(rational(ka,z));self.kb.append(rational(kb,z))
            if k<=1:self.ka_div.append(rational(ka/e,z))
        self.kb1_div=rational((1/sy.Integer(16)-1/qb**2)/e,z)

    def core(self,beta):
        z=beta*beta;e=1-z;x=(1-6*z+z*z)/(1+z)**2
        rg=pe(self.G,x);kg=pe(self.H,x)
        for r,c,kind,v in self.rats:
            y2=1-x*x;dr=1-2*r*c*x+r*r*(2*x*x-1);di=-2*r*c+2*r*r*x
            if kind=='cos':nr=1-r*c*x;ni=-r*c;mm=(1+r*c)/(1+2*r*c+r*r)
            else:nr=r*x;ni=r;mm=-r/(1+2*r*c+r*r)
            den=(2*r*x-c*(1+r*r))**2+(1-c*c)*(1-r*r)**2
            mr=(nr*dr+y2*ni*di)/den;mi=(ni*dr-nr*di)/den
            rg+=v*((2*x*x-1)*(mr-mm)-2*x*y2*mi)
            kg+=v*((2*x*x-1)*mi+2*x*(mr-mm))
        return 16*beta*e*e/(1+z)**5*kg,e*(4/(1+z)**3*rg-self.alpha/(8*(1-e/2)*(3-e)))

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
        ompow={j:cp(om,j,z) for _,j,_,_ in self.deriv}
        d2pow={p:cp(d2,p-1,z) for p,_,_,_ in self.deriv}
        U=zero;V=zero;dbar=4*z;Qbar=4-exact(F(4,3))*e
        R=sum((exact(c)*dbar**(self.ell-k)*(exact(F(4,3))*e)**k for k,c in enumerate(self.coeff)),zero)
        ka=[rat(q,z) for q in self.ka];kb=[rat(q,z) for q in self.kb]
        for p,j,regular,lam in self.deriv:
            vv=cm(cm(base,ompow[j],z),d2pow[p],z)
            # For ell=0 there was no scaled coordinate; undo division by
            # beta^(2p-2) introduced by the common D^2 representation.
            factor=math.factorial(p)*lam*(one if scaled else z**(p-1))
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

    def higher(self,beta,normalized=False):
        z=beta*beta;e=1-z;dbar=4*z;Qbar=4-exact(F(4,3))*e;zero=z*0
        R=sum((exact(c)*dbar**(self.ell-k)*(exact(F(4,3))*e)**k for k,c in enumerate(self.coeff)),zero)
        v=zero
        for row,lam in self.cross:
            D=row['degree']
            for j,bj in enumerate(row['beta']):
                if not bj or j>self.ell:continue
                num=sum((exact(F(self.coeff[k]*math.comb(k,j)))*exact(F(4,3))**(k-j)*dbar**(self.ell-k)*e**k for k in range(j,self.ell+1)),zero)
                v+=lam*exact(F(bj))*num*e**(D-2-j)/(R*Qbar**(D-2-j))
        v+=exact(self.direction[1])
        if not normalized:v*=e*R/(dbar**self.ell*Qbar**3)
        return v

    def forward(self,beta):return self.core(beta)
