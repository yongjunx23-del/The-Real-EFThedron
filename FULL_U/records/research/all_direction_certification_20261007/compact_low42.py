"""All-spin low-energy closure using the reviewed sharp compact envelope."""
from pathlib import Path
import sys,argparse
sys.path.insert(0,str(Path(__file__).parent/'tail'))
from compact_highspin220_sharp import *


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bits',type=int,default=400);args=ap.parse_args()
    ctx.prec=args.bits;L=42;energy=F(5);tic=time.time()
    prep=json.loads((OUT/'prepared_components.json').read_text())
    for p,h in prep['sha256_manifest'].items():assert sha(ROOT/p)==h
    assert prep['objective_exact']==['0','1']
    d=pickle.loads((OUT/'joint_moments.pkl').read_bytes())['data']
    cd=pickle.loads((OUT/'compact_norms.pkl').read_bytes());assert cd['witness_sha256']==sha(WITNESS)
    signpath=ROOT/'output/competitive_20261007/robust/stage1_audit/roy_exact.json'
    sg=json.loads(signpath.read_text())
    assert sg['sha256']==sha(WITNESS) and sg['passed'] and F(sg['metadata']['ymin'])==2*C and F(sg['W_lower_rational'])>0
    co,base=forward_bounds(d)
    pp=coefficient_polynomial(WITNESS);k,z=pp.gens;y=sy.Symbol('y')
    pol=sy.Poly(sy.expand(pp.as_expr().subs(k,y+10)),y)
    cc=[max(bernstein_univariate(pol.nth(j),z,0,sy.Rational(3,8)))for j in range(pol.degree()+1)]
    assert all(v<=0 for v in cc) and cc[5]<0
    power5=-F(cc[5]);M=max([F(0)]+[F(max(bernstein_univariate(pp.as_expr().subs(k,j),z,0,sy.Rational(3,8))))for j in range(10)])
    R,S=cd['roy_abs'],cd['strip_abs'];qmax=energy-C
    mean_lower=lo(mean(L,energy-4));r=mean_lower/L;assert mean_lower>10
    slope=qmax**3*(R*F(3,2)/C+S/C)
    forward_ratio=ratio_bound(co,L,energy-4)
    constant=1+base+forward_ratio+qmax**3*(R*F(13,8)+S*F(15,8))+M
    # All quantities in these last inequalities are exact rational numbers.
    value=constant+slope*L-power5*(r*L-10)**5
    derivative=slope-5*power5*r*(r*L-10)**4
    assert value<-1 and derivative<0
    hashes=dict(prep['sha256_manifest'])
    files=[Path(__file__),Path(forward_bounds.__code__.co_filename),OUT/'compact_norms.pkl',signpath,ROOT/'research/crossing_kernels.py']
    hashes.update({str(p.relative_to(ROOT)):sha(p)for p in files})
    out=dict(status='PASS_ALL_SPIN_LOW_ENERGY_V',witness_sha256=sha(WITNESS),minimum_even_spin=L,
        energy_domain_exact=['4','5'],kappa_exact='1',precision_bits=args.bits,
        V_over_p_upper_exact=str(value),spin_derivative_upper_exact=str(derivative),
        positive_legendre_mean_lower_exact=str(mean_lower),crossing_power5_exact=str(power5),
        crossing_low_k_cap_exact=str(M),source_slope_exact=str(slope),constant_exact=str(constant),
        forward_ratio_upper_exact=str(forward_ratio),sha256_manifest=hashes,seconds=time.time()-tic,
        proof='Same reviewed compact low-energy envelope, replacing endpoint mu8 bymu5, L450/110 byL42. Uniform source resolvent constants13/8 and15/8 remain valid on4<mu<=5. Log-convexity gives m_l>=r*l; the exact upper envelope is decreasing for everyl>=42 by its negative derivative at42 and negative fifth-power term.',
        scope='Every even spin>=42 and4<mu<=5 only. Not a complete coupling certificate.')
    folder=ROOT/'output/competitive_20261007/compact';folder.mkdir(parents=True,exist_ok=True)
    (folder/f'low_allspin42_{args.bits}.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items()if k not in ['sha256_manifest','V_over_p_upper_exact','spin_derivative_upper_exact','constant_exact','forward_ratio_upper_exact']},indent=2))
    print('V upper',float(value),'derivative upper',float(derivative),flush=True)


if __name__=='__main__':main()
