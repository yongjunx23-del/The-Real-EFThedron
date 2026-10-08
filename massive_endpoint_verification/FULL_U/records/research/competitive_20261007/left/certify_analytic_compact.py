"""Reviewed sharp compact envelope for the exact left180 functional."""
from pathlib import Path
import sys,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/all_direction_certification_20261007/tail'))
from compact_highspin220_sharp import *
CACHE=ROOT/'output/all_direction_certification_20261007/tail/left/refined'
DEST=ROOT/'output/competitive_20261007/left/compact'
WITNESS=ROOT/'output/full_direction_20261007/left/angle180_final.npz'
EXPECTED='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bits',type=int,default=400);a=ap.parse_args()
    ctx.prec=a.bits;tic=time.time();assert sha(WITNESS)==EXPECTED
    prep=json.loads((CACHE/'prepared_components.json').read_text());assert prep['objective_exact']==['-1','0']
    for p,h in prep['sha256_manifest'].items():assert sha(ROOT/p)==h,p
    d=pickle.loads((CACHE/'joint_moments.pkl').read_bytes())['data']
    cd=pickle.loads((CACHE/'compact_norms.pkl').read_bytes());assert cd['witness_sha256']==EXPECTED
    signpath=ROOT/'output/full_direction_20261007/left/audit180/roy_exact.json';sg=json.loads(signpath.read_text())
    assert sg['sha256']==EXPECTED and sg['passed'] and F(sg['metadata']['ymin'])==2*C and F(sg['W_lower_rational'])>0
    pp=coefficient_polynomial(WITNESS);k,z=pp.gens;y=sy.Symbol('y');p=sy.Poly(sy.expand(pp.as_expr().subs(k,y+10)),y)
    coeff=[max(bernstein_univariate(p.nth(j),z,0,sy.Rational(3,8)))for j in range(p.degree()+1)]
    assert all(v<=0 for v in coeff) and coeff[5]<0;power5=-F(coeff[5])
    M=max([F(0)]+[F(max(bernstein_univariate(pp.as_expr().subs(k,j),z,0,sy.Rational(3,8))))for j in range(10)])
    co,base=forward_bounds(d);R,S=cd['roy_abs'],cd['strip_abs'];kappa=F(10)
    def low(L,endpoint):
        qmax=endpoint-C;m=lo(mean(L,endpoint-4));r=m/L;assert m>10
        slope=qmax**3*(R*F(3,2)/C+S/C)
        constant=base+ratio_bound(co,L,endpoint-4)+qmax**3*(R*F(13,8)+S*F(15,8))+M
        val=constant+slope*L-power5*(r*L-10)**5
        derivative=slope-5*power5*r*(r*L-10)**4
        assert val<-kappa and derivative<0
        return dict(V_over_p_upper_exact=str(val),spin_derivative_upper_exact=str(derivative),
                    positive_legendre_mean_lower_exact=str(m),constant_exact=str(constant),source_slope_exact=str(slope))
    low42=low(42,F(5));low100=low(100,F(8))
    L=100;qmin=F(20,3);qmax=F(143,3);dmin=F(83,60)
    norm2=F(json.loads((CACHE/'source_norms.json').read_text())['norm2_rational'])
    outside=norm2/(1-(F(317,60)/qmin)**2)
    extra=qmax**3*(R*(F(1,8)+1/qmin+1/dmin)+sum((2*v/(B-A)for(A,B),v in cd['strip_D'].items()),F(0))+S*(F(1,8)+1/qmin+F(1,4)))
    augmented=co.copy();augmented[0]+=extra;augmented[1]+=qmax**3*R/2
    m=mean(L,F(45));assert m>10
    middle=ar(base+outside+M+ratio_bound(augmented,L,F(45)))-ar(power5)*(m-10)**5
    assert middle.is_finite() and up(middle)<-kappa
    manifest=dict(prep['sha256_manifest'])
    for path in [Path(__file__),Path(forward_bounds.__code__.co_filename),CACHE/'prepared_components.json',CACHE/'compact_norms.pkl',signpath,WITNESS]:
        manifest[str(path.relative_to(ROOT))]=sha(path)
    result=dict(status='PASS_LEFT_ANALYTIC_COMPACT_COERCIVITY',witness_sha256=EXPECTED,objective_exact=['-1','0'],
        kappa_exact='10',precision_bits=a.bits,
        low_allspin42=dict(minimum_even_spin=42,mu_domain_exact=['4','5'],**low42),
        high_allspin100=dict(minimum_even_spin=100,mu_domain_exact=['4','49'],low_mu4_8=low100,middle_mu8_49_upper_exact=str(up(middle))),
        crossing_power5_exact=str(power5),crossing_low_k_cap_exact=str(M),
        sha256_manifest=manifest,seconds=time.time()-tic,
        scope='V/p<-10 on all even spins>=42,4<mu<=5 and all even spins>=100,4<mu<=49. Finite42..98,5..49 and all UV remain separate obligations.',
        proof='Same reviewed sharp compact envelope, with fresh exact left-source contractions and objective b=0. The g0 objective is included in alpha. All-spin monotonicity of Legendre mean and denominator monomial ratios is unchanged.')
    DEST.mkdir(parents=True,exist_ok=True);(DEST/f'analytic_compact_{a.bits}.json').write_text(json.dumps(result,indent=2))
    print('PASS',a.bits,'low42',float(F(low42['V_over_p_upper_exact'])),'high100',float(middle.upper()),'seconds',result['seconds'],flush=True)
if __name__=='__main__':main()
