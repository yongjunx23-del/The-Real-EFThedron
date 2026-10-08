"""Exact standalone Bernstein panel/maxima ledger for the95 adaptation."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[2];CODE=ROOT/'research/six_targets_20261007'
sys.path.insert(0,str(CODE/'tail_upper95'))
from angleupper95_campaign import WITNESS,SHA,OBJECTIVE,OUT
from matrix_crossing_series import coefficient_polynomial,bernstein_univariate,shifted_certificate
import sympy as sy
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
p=coefficient_polynomial(WITNESS);k,z=p.gens;y=sy.Symbol('y');records=[];maxima={}
for cutoff,panels in [(10,4),(64,1)]:
 shifted=sy.Poly(sy.expand(p.as_expr().subs(k,y+cutoff)),y);rows=[];caps=[]
 for j in range(shifted.degree()+1):
  intervals=[]
  for i in range(panels):
   a=sy.Rational(3*i,8*panels);b=sy.Rational(3*(i+1),8*panels)
   cs=bernstein_univariate(shifted.nth(j),z,a,b)
   intervals.append(dict(z_domain_exact=[str(a),str(b)],coefficients_exact=list(map(str,cs)),upper_exact=str(max(cs))))
  cap=max(F(r['upper_exact'])for r in intervals);assert cap<=0;caps.append(cap)
  rows.append(dict(shifted_power=j,panels=intervals,global_upper_exact=str(cap)))
 assert caps[5]<0;maxima[cutoff]=caps;records.append(dict(cutoff=cutoff,records=rows))
assert shifted_certificate(p,10,4)is not None
assert shifted_certificate(p,10,1)is None
large_path=OUT/'large_tau_320.json';large=json.loads(large_path.read_text())
assert large['witness_sha256']==SHA and large['status']=='PASS_ANGLEupper95_UNBOUNDED_TAU'
assert list(map(F,large['crossing_negative_shifted_powers_exact']))==[-v for v in maxima[10]]
paths=[Path(__file__),WITNESS,CODE/'tail_upper95/angleupper95_campaign.py',CODE/'tail_upper95/angleupper95_models.py',CODE/'tail_upper95/certify_angleupper95_large_tau_moments.py',ROOT/'research/lower_g0_upper_g2_20261006/matrix_crossing_series.py',OUT/'prepared_components.json',large_path,OUT.parent/'audits/crossing_exact.json']
result=dict(status='PASS_EXACT_FOUR_PANEL_CROSSING_AND_UNCHANGED_K64_TAIL',witness_sha256=SHA,objective_exact=list(map(str,OBJECTIVE)),z_domain_exact=['0','3/8'],records=records,single_panel_cutoff10_inconclusive=True,k64_single_panel_unchanged=True,large_tau_320_coefficients_match=True,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in paths},scope='Exact polynomial coefficient bounds. Four adjacent panels cover the full z interval; their maximum is a uniform upper coefficient. This is not alone a support or complete UV bound.')
out=OUT/'fourpanel_crossing_review.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(result['status'])
