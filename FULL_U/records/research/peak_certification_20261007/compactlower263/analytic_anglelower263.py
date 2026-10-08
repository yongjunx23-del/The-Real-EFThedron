from pathlib import Path
import sys
_ADAPT_ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(_ADAPT_ROOT/'research/peak_certification_20261007/tail_lower263'),str(_ADAPT_ROOT/'research/boundary_length_20261007/compact')]
"""Same reviewed compact envelope on smaller exact energy intervals.

Only qmin/qmax and the positive-Legendre denominator/mean endpoint change.
This keeps the physical-source, PV, crossing and derivative inequalities
unchanged while reducing the finite-spin bridge in each interval.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
from analytic_angle85 import *
import anglelower263_campaign as campaign
WITNESS=campaign.WITNESS;EXPECTED=campaign.SHA;OBJECTIVE=campaign.OBJECTIVE;CACHE=campaign.OUT
source_manifest=campaign.source_manifest
DEST=ROOT/'output/peak_certification_20261007/lower263/compact'
DEST.mkdir(parents=True,exist_ok=True)
PAIRS=[('4','19/4'),('19/4','5'),('5','11/2'),('11/2','6'),('6','13/2'),('13/2','7'),('7','15/2'),('15/2','8'),
       ('8','10'),('10','12'),('12','16'),('16','20'),('20','24'),('24','30'),('30','36'),('36','42'),('42','49')]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bits',type=int,default=400);a=ap.parse_args();ctx.prec=a.bits;tic=time.time();manifest=source_manifest()
    d=pickle.loads((CACHE/'joint_moments.pkl').read_bytes())['data'];cd=pickle.loads((CACHE/'compact_norms.pkl').read_bytes());assert cd['witness_sha256']==EXPECTED
    signpath=campaign.AUDITS/'roy_exact.json';sg=json.loads(signpath.read_text());assert sg['sha256']==EXPECTED and sg['passed']and F(sg['metadata']['ymin'])==2*C and F(sg['W_lower_rational'])>0
    pp=coefficient_polynomial(WITNESS);k,z=pp.gens;y=sy.Symbol('y');p=sy.Poly(sy.expand(pp.as_expr().subs(k,y+10)),y)
    coeff=[max(bernstein_univariate(p.nth(j),z,0,sy.Rational(3,8)))for j in range(p.degree()+1)];assert all(v<=0 for v in coeff)and coeff[5]<0;powers=[-F(v) for v in coeff]
    M=max([F(0)]+[F(max(bernstein_univariate(pp.as_expr().subs(k,j),z,0,sy.Rational(3,8))))for j in range(10)])
    co,base=forward_bounds(d);base+=OBJECTIVE[1];R,S=cd['roy_abs'],cd['strip_abs'];norm2=F(json.loads((CACHE/'source_norms.json').read_text())['norm2_rational']);records=[]
    for aa,bb in PAIRS:
        aa,bb=F(aa),F(bb);qmax=bb-C;selected=None
        for L in range(42,452,2):
            m=mean(L,bb-4)
            if m.lower()<=10:continue
            if bb<=8:
                ml=lo(m);r=ml/L;slope=qmax**3*(R*F(3,2)/C+S/C)
                constant=base+ratio_bound(co,L,bb-4)+qmax**3*(R*F(13,8)+S*F(15,8))+M
                value=constant+slope*L-sum((v*(ml-10)**j for j,v in enumerate(powers) if j>=1),F(0));derivative=slope-sum((j*v*r*(ml-10)**(j-1) for j,v in enumerate(powers) if j>=1),F(0))
                good=value<-1 and derivative<0;method='low';upval=value
            else:
                assert aa>=8;qmin=aa-C;dmin=F(83,60)
                outside=norm2/(1-(F(317,60)/qmin)**2)
                extra=qmax**3*(R*(F(1,8)+1/qmin+1/dmin)+sum((2*v/(B-A)for(A,B),v in cd['strip_D'].items()),F(0))+S*(F(1,8)+1/qmin+F(1,4)))
                aug=co.copy();aug[0]+=extra;aug[1]+=qmax**3*R/2
                value=ar(base+outside+M+ratio_bound(aug,L,bb-4))-sum((ar(v)*(m-10)**j for j,v in enumerate(powers) if j>=1),arb(0))
                good=value.is_finite()and value.upper()<-1;method='physical';upval=up(value);derivative=None
            if good:
                selected=dict(mu_domain_exact=[str(aa),str(bb)],minimum_even_spin=L,V_over_p_upper_exact=str(upval),
                  method=method,spin_derivative_upper_exact=None if derivative is None else str(derivative),mean_lower_exact=str(lo(m)),finite_bridge_spins=list(range(42,L,2)))
                break
        assert selected is not None,(aa,bb);records.append(selected);print('PASS',aa,bb,'L',L,'upper',float(upval),flush=True)
    for path in[Path(__file__),ROOT/'research/boundary_length_20261007/compact/analytic_angle85.py',ROOT/'research/peak_certification_20261007/tail_lower263/anglelower263_campaign.py',Path(forward_bounds.__code__.co_filename),CACHE/'prepared_components.json',CACHE/'compact_norms.pkl',signpath,WITNESS]:manifest[str(path.relative_to(ROOT))]=sha(path)
    out=dict(status='PASS_PIECEWISE_ALLSPIN_COMPACT_ENVELOPES',witness_sha256=EXPECTED,objective_exact=list(map(str,OBJECTIVE)),kappa_exact='1',precision_bits=a.bits,
      crossing_negative_shifted_powers_exact=list(map(str,powers)),crossing_low_k_cap_exact=str(M),records=records,sha256_manifest=manifest,seconds=time.time()-tic,
      proof='All nonnegative shifted powers j>=1 retained: Jensen bounds E[(k-10)_+^j] below by (E[k]-10)^j whenever E[k]>10. Uniform exact Bernstein bounds in z yield these coefficients; negative constant discarded. Reviewed sharp envelope restricted to exact subintervals. Low formula bounds all4<mu<=upper; physical formula applies onmu>=8 with qmin=lower-4/3,qmax=upper-4/3 and all other conservative constants unchanged. Each positive denominator monomial ratio decreases withspin, and the Legendre mean increases. This certificate does not include the listed finite bridge spins.')
    (DEST/f'analytic_piecewise_{a.bits}.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
