"""Exact source/crossing prerequisites for arbitrary-direction220-row tails."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,pickle,re,sys,time
ROOT=Path(__file__).resolve().parents[3];sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research')]
from tail_certify_joint_limit import derive
from matrix_source_norms import source_norms
from matrix_crossing_series import coefficient_polynomial,shifted_certificate,bernstein_univariate
import numpy as np
import sympy as sy
OUT=ROOT/'output/all_direction_certification_20261007/tail'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('witness',type=Path);ap.add_argument('--label',required=True);a=ap.parse_args();tic=time.time();w=a.witness.resolve()
    out=OUT/a.label;out.mkdir(parents=True,exist_ok=True);data=np.load(w);direction=list(map(lambda x:F(float(x)),data['direction']))
    cache=out/'joint_moments.pkl'
    if cache.exists():
        d=pickle.loads(cache.read_bytes());assert d['witness_sha256']==sha(w)
    else:
        d=dict(witness_sha256=sha(w),objective_exact=list(map(str,direction)),data=derive(w));cache.write_bytes(pickle.dumps(d))
    normfile=out/'source_norms.json'
    if normfile.exists():norm=json.loads(normfile.read_text());assert norm['sha256']==sha(w)
    else:norm=source_norms(w,16,16);normfile.write_text(json.dumps(norm,indent=2))
    p=coefficient_polynomial(w);k,z=p.gens;cert=shifted_certificate(p,10,1);assert cert is not None
    comparisons=[]
    for cut in [10,12,16,20,24,32,48,64]:
        # Coefficients in y=k-cut. If every positive-y coefficient is
        # nonpositive, the constant Bernstein maximum bounds all k>=cut.
        y=sy.Symbol('y');q=sy.Poly(sy.expand(p.as_expr().subs(k,cut+y)),y)
        for upper in [sy.Rational(3,143),sy.Rational(3,8)]:
            bounds=[max(bernstein_univariate(q.nth(j),z,0,upper)) for j in range(q.degree()+1)]
            comparisons.append(dict(k_min=cut,z_max=str(upper),higher_shifted_coefficients_nonpositive=all(v<=0 for v in bounds[1:]),constant_upper_exact=str(bounds[0])))
    files=[w,cache,normfile,Path(__file__),ROOT/'research/lower_g0_upper_g2_20261006/tail_certify_joint_limit.py',ROOT/'research/lower_g0_upper_g2_20261006/matrix_source_norms.py',ROOT/'research/lower_g0_upper_g2_20261006/matrix_crossing_series.py']
    record=dict(status='PASS',witness_sha256=sha(w),objective_exact=list(map(str,direction)),crossing_k10_certificate=cert,crossing_uniform_comparisons=comparisons,
        norm1_approx=norm['norm1_approx'],norm2_approx=norm['norm2_approx'],seconds=time.time()-tic,
        objective_warning='Existing joint/finiteUV kernels omit objective b; add the exact second objective coefficient to every resulting V/p.',
        sha256_manifest={str(p.relative_to(ROOT)):sha(p) for p in files})
    (out/'prepared_components.json').write_text(json.dumps(record,indent=2))
    print(json.dumps({k:v for k,v in record.items() if k not in ['crossing_k10_certificate','sha256_manifest']},indent=2))
if __name__=='__main__':main()
