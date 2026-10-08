"""Exact rational kernel construction without symbolic expansion.

First transfer derivative/forward harmonics/X4 only. The resulting four
polynomials encode the positively congruent full disk coefficient matrix.
"""
from fractions import Fraction as F
from exact_bernstein import ROOT,q,fmpq,fmpq_poly,bernstein
import json,time

class R:
    def __init__(self,n=0,d=1):
        if isinstance(n,R):self.n,self.d=n.n,n.d;return
        self.n=n if isinstance(n,fmpq_poly) else fmpq_poly([q(n)])
        self.d=d if isinstance(d,fmpq_poly) else fmpq_poly([q(d)])
        assert self.d
        if not self.n:self.d=fmpq_poly([1]);return
        g=self.n.gcd(self.d)
        self.n=divmod(self.n,g)[0];self.d=divmod(self.d,g)[0]
        c=self.d.leading_coefficient();self.n=self.n/c;self.d=self.d/c
    def __add__(self,o):
        o=R(o);return R(self.n*o.d+o.n*self.d,self.d*o.d)
    __radd__=__add__
    def __neg__(self):return R(-self.n,self.d)
    def __sub__(self,o):return self+-R(o)
    def __rsub__(self,o):return R(o)+-self
    def __mul__(self,o):
        o=R(o);return R(self.n*o.n,self.d*o.d)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=R(o);return R(self.n*o.d,self.d*o.n)
    def __rtruediv__(self,o):return R(o)/self
    def __pow__(self,n):
        assert isinstance(n,int)
        return R(self.n**n,self.d**n) if n>=0 else R(self.d**(-n),self.n**(-n))

def compose(p,x):
    v=R(0)
    for a in reversed(p.coeffs()):v=v*x+R(a)
    return v

def legendre(n):
    z=fmpq_poly([0,1]);p0=fmpq_poly([1])
    if n==0:return p0
    p1=z
    for k in range(1,n):p0,p1=p1,((2*k+1)*z*p1-k*p0)/(k+1)
    return p1

def prove_positive(p):
    stack=[(F(0),F(1),0)];cover=[]
    while stack:
        a,b,depth=stack.pop();low=min(bernstein(p,a,b))
        if low>0:cover.append((str(a),str(b),str(low)));continue
        assert depth<20
        m=(a+b)/2;stack.extend([(m,b,depth+1),(a,m,depth+1)])
    return cover

def build(ell,names,weights,a=0,b=1):
    beta=R(fmpq_poly([0,1]));r=beta**2;one=1-r
    E=4/one;delta=4*r/one;Q=E-F(4,3);e2=one**2;x=one/(1+r)
    f=32/(E-2)**3;c=2*x*x-1
    fp=f*(F(1,2)-F(3,2)/(E-2));cp=4*delta/(E-2)**3
    ka=2/(E-2)-1/Q-1/(E-F(8,3));kap=-1/(E-2)**2+1/(E-F(8,3))**2
    kbp=1/E**2-1/(E-F(8,3))**2
    U0=16*beta*one**4/(1+r)**5
    top=max([int(n[6:]) if n.startswith('Dt1H_j') else int(n[1:]) for n in names if n.startswith('H') or n.startswith('Dt1H_j')]+[1])
    ts=[R(1),c];us=[R(1),2*c]
    for j in range(2,top+1):
        ts.append(2*c*ts[-1]-ts[-2]);us.append(2*c*us[-1]-us[-2])
    Z=Q/delta;pl=legendre(ell);P=compose(pl,Z);p=P/Q**3
    CP=8/((3*E-8)*(3*E-4))*(1/(E-2)-4*P/E)
    U=F(a)*U0;V=F(a)*(f*c+CP)+F(b)*p;L=ell*(ell+1)
    for name,weight in zip(names,weights):
        weight=F(weight)
        if name=='X4':
            cross=(2*Z*Z*compose(pl.derivative().derivative(),Z)-3*Z*compose(pl.derivative(),Z))/Q**5
            V+=weight*cross;continue
        sub=name=='Dt1G_regular';deriv=name.startswith('Dt')
        j=1 if sub else int(name[6:] if deriv else name[1:])
        assert sub or j>=2
        if deriv:
            uu=U0*(L/delta*us[j-1]+((2-x-3*x*x)*us[j-1]-j*ts[j])/(4*(1+x)))
            vv=(fp+L/delta*f)*ts[j]+f*cp*j*us[j-1]
            if sub:vv-=kap+L/delta*ka+P*kbp
        else:uu=U0*us[j-1];vv=f*ts[j]
        U+=weight*uu;V+=weight*vv
    off=-beta**ell*U/e2;bb=-2*beta**(2*ell)*V/e2
    D=divmod(off.d*bb.d,off.d.gcd(bb.d))[0]
    if D(q(F(1,2)))<0:D=-D
    bbquot,bbrem=divmod(D,bb.d);offquot,offrem=divmod(D,off.d)
    assert not bbrem and not offrem
    pp=[fmpq_poly([0]*(2*ell)+[1])*D,bb.n*bbquot,off.n*offquot,D]
    n=max(p.degree() for p in pp)+1
    arr=[[p[k] for k in range(n)] for p in pp]
    return arr,prove_positive(D.coeffs())

def self_test():
    d=json.loads((ROOT/'output/larger_grams_20261006/certificate_m011_head.json').read_text())
    tic=time.perf_counter()
    for row in d['finite_spins']:
        arr,cover=build(row['ell'],d['names'],d['multipliers'])
        old=[fmpq_poly(list(map(fmpq,p))) for p in row['polynomials_A_B_C_D']]
        new=[fmpq_poly(p) for p in arr]
        for j in range(3):assert new[j]*old[3]==old[j]*new[3]
        print('MATCH',row['ell'],len(arr[0])-1,time.perf_counter()-tic,flush=True)

if __name__=='__main__':self_test()
