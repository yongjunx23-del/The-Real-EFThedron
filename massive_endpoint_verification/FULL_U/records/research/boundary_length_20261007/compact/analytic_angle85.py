"""Reviewed sharp compact envelope for the exact frozen85-degree witness."""
from pathlib import Path
import sys,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'research/boundary_length_20261007/tail'),str(ROOT/'research/all_direction_certification_20261007/tail')]
import angle85_campaign as campaign
from compact_highspin220_sharp import *
WITNESS=campaign.WITNESS;EXPECTED=campaign.SHA;OBJECTIVE=campaign.OBJECTIVE;CACHE=campaign.OUT
source_manifest=campaign.source_manifest
DEST=ROOT/'output/boundary_length_20261007/compact/angle85'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bits',type=int,default=400);ap.add_argument('--L',type=int,default=110)
    ap.add_argument('--kappa',default='1/4');ap.add_argument('--low-endpoint',default='9/2');ap.add_argument('--screen',action='store_true');a=ap.parse_args()
    ctx.prec=a.bits;tic=time.time();assert sha(WITNESS)==EXPECTED;manifest=source_manifest()
    d=pickle.loads((CACHE/'joint_moments.pkl').read_bytes())['data'];cd=pickle.loads((CACHE/'compact_norms.pkl').read_bytes());assert cd['witness_sha256']==EXPECTED
    signpath=ROOT/'output/full_direction_20261007/right/fresh_candidate_audits/angle85_220/roy_exact.json';sg=json.loads(signpath.read_text())
    assert sg['sha256']==EXPECTED and sg['passed']and F(sg['metadata']['ymin'])==2*C and F(sg['W_lower_rational'])>0
    pp=coefficient_polynomial(WITNESS);k,z=pp.gens;y=sy.Symbol('y');p=sy.Poly(sy.expand(pp.as_expr().subs(k,y+10)),y)
    coeff=[max(bernstein_univariate(p.nth(j),z,0,sy.Rational(3,8)))for j in range(p.degree()+1)]
    assert all(v<=0 for v in coeff)and coeff[5]<0;power5=-F(coeff[5])
    M=max([F(0)]+[F(max(bernstein_univariate(pp.as_expr().subs(k,j),z,0,sy.Rational(3,8))))for j in range(10)])
    co,base=forward_bounds(d);base+=OBJECTIVE[1];R,S=cd['roy_abs'],cd['strip_abs'];kappa=F(a.kappa)
    def low(L,endpoint):
        qmax=endpoint-C;m=lo(mean(L,endpoint-4));r=m/L;assert m>10
        slope=qmax**3*(R*F(3,2)/C+S/C)
        constant=base+ratio_bound(co,L,endpoint-4)+qmax**3*(R*F(13,8)+S*F(15,8))+M
        val=constant+slope*L-power5*(r*L-10)**5;derivative=slope-5*power5*r*(r*L-10)**4
        return dict(V_over_p_upper_exact=str(val),spin_derivative_upper_exact=str(derivative),positive_legendre_mean_lower_exact=str(m),constant_exact=str(constant),source_slope_exact=str(slope))
    low_endpoint=F(a.low_endpoint);assert 4<low_endpoint<=5
    low42=low(42,low_endpoint);lowL=low(a.L,F(8));qmin=F(20,3);qmax=F(143,3);dmin=F(83,60)
    norm2=F(json.loads((CACHE/'source_norms.json').read_text())['norm2_rational']);outside=norm2/(1-(F(317,60)/qmin)**2)
    extra=qmax**3*(R*(F(1,8)+1/qmin+1/dmin)+sum((2*v/(B-A)for(A,B),v in cd['strip_D'].items()),F(0))+S*(F(1,8)+1/qmin+F(1,4)))
    augmented=co.copy();augmented[0]+=extra;augmented[1]+=qmax**3*R/2;m=mean(a.L,F(45));assert m>10
    middle=ar(base+outside+M+ratio_bound(augmented,a.L,F(45)))-ar(power5)*(m-10)**5
    print('SCREEN',a.L,'low42',float(F(low42['V_over_p_upper_exact'])),'lowL',float(F(lowL['V_over_p_upper_exact'])),'middle',float(middle.upper()),flush=True)
    if a.screen:return
    assert middle.is_finite()and up(middle)<-kappa
    for row in [low42,lowL]:assert F(row['V_over_p_upper_exact'])<-kappa and F(row['spin_derivative_upper_exact'])<0
    for path in[Path(__file__),Path(forward_bounds.__code__.co_filename),CACHE/'prepared_components.json',CACHE/'compact_norms.pkl',signpath,WITNESS,ROOT/'research/boundary_length_20261007/tail/angle85_campaign.py']:
        manifest[str(path.relative_to(ROOT))]=sha(path)
    result=dict(status='PASS_ANGLE85_ANALYTIC_COMPACT_COERCIVITY',witness_sha256=EXPECTED,objective_exact=list(map(str,OBJECTIVE)),kappa_exact=str(kappa),precision_bits=a.bits,
      low_allspin42=dict(minimum_even_spin=42,mu_domain_exact=['4',str(low_endpoint)],**low42),
      high_allspin=dict(minimum_even_spin=a.L,mu_domain_exact=['4','49'],low_mu4_8=lowL,middle_mu8_49_upper_exact=str(up(middle))),
      crossing_power5_exact=str(power5),crossing_low_k_cap_exact=str(M),sha256_manifest=manifest,seconds=time.time()-tic,
      proof='Same reviewed sharp compact envelope, exact85 source contractions and exact objective b added. The g0 objective is included in alpha. Mean/monomial-ratio allspin monotonicity is unchanged.',
      scope=f'All even>=42 on4<mu<={low_endpoint} and all even>={a.L} on4<mu<=49; the finite compact bridge and UV remain separate.')
    DEST.mkdir(parents=True,exist_ok=True);(DEST/f'analytic_L{a.L}_{a.bits}.json').write_text(json.dumps(result,indent=2));print('PASS',a.bits,'seconds',result['seconds'],flush=True)
if __name__=='__main__':main()
