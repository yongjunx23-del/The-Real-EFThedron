"""Arb component proof: V_l/p_l is strictly negative for 4<mu<=49.

Requires the exact expanded Roy-source sign certificate W>=0 at t<=-3/5.
Uses exact source norms/moments, physical Legendre bounds, and X4 coercivity.
"""
import argparse,pickle,json,hashlib
from pathlib import Path
from fractions import Fraction as F
import sys
from certify_uv import ar,arb,ctx,ROOT,leg_coeffs,poly

def P(L,z):
 # Positive Taylor series P_L(1+2a/[L(L+1)]).
 return poly(leg_coeffs(L),(z-1)*ar(F(L*(L+1),2)))

def forward_abs(d):
 # Uniform over all mu>=4, conservative complex rather than real cancellation.
 base=abs(d['alpha'])+sum((abs(v) for j,v in d['harm']),F(0))
 base+=sum((4*abs(v)/(1-r)**2 for r,co,v in d['rats']),F(0))
 return F(2048,27)*base+38*abs(d['alpha'])

def run_scoped(L,d,cd,step):
 assert L>=42 and L%2==0 and d['lx']<0
 eta=-d['lx'];T=L*(L+1);nf=forward_abs(d)
 p8=P(L,ar(F(5,3)));e0=ar(F(9,20)**L)+1/p8
 edd=ar(F(5,4)*T*F(9,20)**L)+ar(F(5,6)*T)/p8
 # Q^3 K2(t) <= 2d^2 Q0^2/(Q0^2-d^2), d<=29/15 and Q0=8/3.
 k2neg=F(2)*F(29,15)**2*F(8,3)**2/(F(8,3)**2-F(29,15)**2)
 lowmax=arb(-10**20);lowarg=None;lowcells=0;derivmax=arb(-10**20)
 step=F(str(step));a=F(4)
 while a<5:
  b=min(F(8),a+step);m=arb(ar((a+b)/2),ar((b-a)/2));q=m-ar(F(4,3));delta=m-4
  # g(delta)=max(delta,6/5-delta), exact minimum on interval.
  dl,dh=a-4,b-4
  gmin=F(3,5) if dl<=F(3,5)<=dh else min(max(dl,F(6,5)-dl),max(dh,F(6,5)-dh))
  roy=q**3*ar(cd['roy_pos'])*((1+1/p8)/ar(gmin)+3/q)
  roy+=ar(cd['roy_neg']*k2neg)+q**3*ar(cd['roy_abs'])*(edd+2*e0)
  strips=q**3*ar(cd['strip_pos'])*(1/(12-m)+3/q)
  D=(q*q-delta*delta).sqrt()
  C1=q**3*ar(2*cd['strip_neg'])/D
  strips+=q**3*ar(cd['strip_neg'])*(1/(12-m)+2/m+1/q)
  c0=ar(F(128,27)*abs(d['alpha']))+ar(nf)/p8+roy+strips
  aa=ar(F(4,3)*T/(b-4));cs=leg_coeffs(L)
  rr=(poly([k*v for k,v in enumerate(cs)],aa)/poly(cs,aa))/L
  coef=1+4/q
  assert rr*L>coef/2
  jchi=(rr*rr*L*L-rr*L*coef)/ar(F(32,9))
  jderiv=(2*rr*rr*L-rr*coef)/ar(F(32,9))
  classical=ar(L*(2*L-5))/q**2
  classical_derivative=ar(4*L-5)/q**2
  if jchi.lower()>classical.lower():
   upper=c0+C1*L-ar(eta)*jchi
   deriv=C1-ar(eta)*jderiv
  else:
   upper=c0+C1*L-ar(eta)*classical
   deriv=C1-ar(eta)*classical_derivative
  if not deriv<0:raise RuntimeError('Low-energy spin monotonicity failed '+str((a,b,deriv)))
  if not upper<0:raise RuntimeError('Low-energy bound failed '+str((a,b,upper)))
  if upper.upper()>lowmax.upper():lowmax=upper;lowarg=[float(a),float(b)]
  if deriv.upper()>derivmax.upper():derivmax=deriv
  lowcells+=1;a=b
 middlemax=arb(-10**20);middlearg=None;middlecells=0;a=F(8);N=len(cd['moments'])
 while a<49:
  b=min(F(49),a+step,F(12) if a<12 else F(49));m=arb(ar((a+b)/2),ar((b-a)/2));q=m-ar(F(4,3));delta=m-4
  # K2 = sum_{k>=0} 2(t-c)^(2k+2)/Q^(2k).
  signed=arb(0)
  for moment in reversed(cd['moments']):signed=signed/(q*q)+ar(moment)
  r=ar(F(317,60))/q
  rem=ar(d['norm2'])*r**(2*N)/(1-r*r)
  A=-ar(F(32,9)*d['alpha'])*q*q/(m*(m-ar(F(8,3))))-signed+rem
  pl=P(L,ar((b-F(4,3))/(b-4)))
  if a<12:
   nr=ar(cd['roy_abs'])*(ar(F(T,2))/delta+1/m+1/q+1/(m-ar(F(397,60))))
  else:
   nr=ar(cd['roy_abs'])*(1/(m-ar(F(119,10)))+1/m+1/q+1/(m-ar(F(397,60))))
  ns=ar(sum((2*v/(B-A) for (A,B),v in cd['strip_D'].items()),F(0)))
  ns+=ar(cd['strip_abs'])*(1/m+1/q+1/delta)
  # Positive Legendre coefficients make m_L=E[k] increasing in Z.
  # Log-convexity in spin gives m_l>=m_L for every l>=L.
  aa=ar(F(4,3)*T/(b-4));cs=leg_coeffs(L)
  mean=poly([k*v for k,v in enumerate(cs)],aa)/poly(cs,aa)
  coef=1+4/q
  assert mean>coef/2
  jensen=(mean*mean-mean*coef)/ar(F(32,9))
  crude=ar(L*(2*L-5))/q**2
  # Each is a valid lower bound; select by exact directed lower endpoint.
  chi=arb(max(jensen.lower(),crude.lower()))
  upper=A+(ar(nf)+q**3*(nr+ns))/pl-ar(eta)*chi
  if not upper<0:raise RuntimeError('Middle-energy bound failed '+str((a,b,upper)))
  if upper.upper()>middlemax.upper():middlemax=upper;middlearg=[float(a),float(b)]
  middlecells+=1;a=b
 return dict(L=L,low_energy_upper=str(lowmax.upper()),low_energy_worst_cell=lowarg,low_energy_cells=lowcells,
  low_energy_spin_derivative_upper=str(derivmax.upper()),middle_energy_upper=str(middlemax.upper()),
  middle_energy_worst_cell=middlearg,middle_energy_cells=middlecells,moment_terms=N,
  arithmetic='Directed Arb inequalities with exact rational source moments and Bernstein norms',
  domain=f'all even l>={L}, 4<mu<=5 or8<=mu<=49',
  prerequisite='Exact positivity of combined Roy W on t<=-3/5, certificate expanded_roy_sign.json',
  claim='Strict negative V_l/p_l on two stated compact-energy subdomains only',crossing_lower_bound='Jensen positive Legendre mean plus classical root-pair bound')
