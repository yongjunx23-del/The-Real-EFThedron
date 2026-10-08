"""All-spin support bound allowing positive V on bounded joint-tau bands.

For eta>0, disk unitarity implies
  h(U,V) <= U^2/(2 eta p) + 2 p (V/p+eta)_+.
On mu>=M, tau=l(l+1)/(mu-4), and even l>=L,
  sum_l integral_[a,b] p dnu
  <= (b-a)b I0(2sqrt((4/3)b))/(2 pi sqrt(1-4/M) L^2).
The latter uses Q>=l(l+1)/tau, the positive Legendre Taylor majorant,
and sum_all_integer_l>=L (2l+1)/(l^2(l+1)^2)=1/L^2.
All omitted tau>=256 and compact-energy regions must obey V/p<=-eta.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, pickle, sys, time
ROOT=Path(__file__).resolve().parents[3]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/tail_certification'),str(ROOT/'tmp/certification_deps')]
from certify_uv import ar,arb,ctx,errors,leg_coeffs,poly
from compact_jensen import run as compact_run
OUT=ROOT/'output/all_direction_certification_20261007/tail/old150_extension'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def up(q):
    assert q.is_finite();m,e=q.upper().man_exp();return F(int(m))*F(2)**int(e)
def uv_bounds(d,cases,L,M,step):
    coeff=leg_coeffs(L);d2=[coeff[k]*k*(k-1) for k in range(2,len(coeff))]
    ea,en=errors(d,M)
    A=-F(32,9)*d['alpha']-d['signed2']+ea
    B=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+en
    pp=poly(coeff,ar(F(4,3)*256));nmax=B+2*sum((v for v,t in d['correction']),F(0))
    big=ar(A)+ar(max(nmax,F(0)))/pp
    am=ar(F(4,3)*256)
    mean=poly([k*v for k,v in enumerate(coeff)],am)/poly(coeff,am)
    cc=1+ar(F(4)/(F(M)-F(4,3)))
    assert mean>cc/2
    jlarge=(mean*mean-mean*cc)/ar(F(32,9))
    large={b:up(big+ar(d['lx']+delta)*jlarge+ar(b)) for b,delta in cases}
    cuts=[F(0)]
    while cuts[-1]<8:cuts.append(min(F(8),cuts[-1]+step))
    cuts.extend(F(i) for i in range(9,257));rho=ar(1-F(4,M)).sqrt();records=[]
    for lo,hi in zip(cuts,cuts[1:]):
        mcell=max(F(M),F(4)+F(L*(L+1))/hi)
        eac,enc=errors(d,mcell)
        Ac=-F(32,9)*d['alpha']-d['signed2']+eac
        Bc=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+enc
        al=ar(F(4,3)*lo);ah=ar(F(4,3)*hi)
        pl=poly(coeff,al);pu=(2*ah.sqrt()).bessel_i(0)
        cor=sum((min(2*v,hi*t) for v,t in d['correction']),F(0));num=Bc+cor
        ratio=ar(num)/(pl if num>=0 else pu)
        if lo==0 or lo>=8:xlo=arb(0)
        else:
            p2=poly(d2,al);pa=(2*ah.sqrt()).bessel_i(1)/ah.sqrt()
            neg=ar(F(3,2)*hi/(mcell-F(4,3)))*pa/pl
            raw=ar(lo*lo/F(2))*p2/pu-neg;xlo=raw if raw>0 else arb(0)
        mean=poly([k*v for k,v in enumerate(coeff)],al)/pl
        cc=1+ar(F(4)/(mcell-F(4,3)))
        if mean>cc/2:
            jlower=(mean*mean-mean*cc)/ar(F(32,9))
            if jlower.lower()>xlo.lower():xlo=arb(jlower.lower())
        bounds={str(b):str(up(ar(Ac)+ratio+ar(d['lx']+delta)*xlo+ar(b))) for b,delta in cases}
        mass=up(ar((hi-lo)*hi)*pu/(2*arb.pi()*rho*L**2))
        records.append(dict(tau_a=str(lo),tau_b=str(hi),p_mass_upper_exact=str(mass),V_over_p_upper_by_slope=bounds))
    return large,records
