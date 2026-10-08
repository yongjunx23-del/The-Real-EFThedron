"""Same certified source integrator with exact combined derivative polynomials."""
from finite220 import *
from forward_combined import ForwardCombined
class FiniteCombined(FiniteSpin):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.base=ForwardCombined(self.path,self.ell)

def comparison():
    path=ROOT/'output/full_direction_20261007/right/g2_guarded192_refined.npz';ctx.prec=256
    for ell in (0,2,4,14,40):
        a=Forward220(path,ell);b=ForwardCombined(path,ell)
        for mu in (4.00001,4.01,5,10,100,1e6):
            beta=exact(F(float(math.sqrt(1-4/mu))))
            x=a.derivative(beta);y=b.derivative(beta)
            assert all(u.overlaps(v) for u,v in zip(x,y)),(ell,mu,x,y)
        print('OVERLAP',ell,flush=True)
    for aa,bb in ((F(1,10000),F(1,9000)),(F(1,100),F(11,1000)),(F(1,10),F(101,1000))):
        beta=arb(exact((aa+bb)/2),exact((bb-aa)/2));a=Forward220(path,2);b=ForwardCombined(path,2)
        print('WIDTHS',aa,bb,[float(v.rad()) for v in a.derivative(beta,True)],[float(v.rad()) for v in b.derivative(beta,True)],flush=True)
if __name__=='__main__':comparison()
