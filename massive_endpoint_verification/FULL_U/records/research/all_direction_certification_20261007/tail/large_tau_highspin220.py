"""Unbounded joint-tau closure for the robust220 witness."""
from compact_highspin220 import *
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lmin',type=int,default=450);ap.add_argument('--precision',type=int,default=320);a=ap.parse_args();ctx.prec=a.precision;tic=time.time();L=a.lmin
 prep=json.loads((OUT/'prepared_components.json').read_text())
 for p,h in prep['sha256_manifest'].items():assert sha(ROOT/p)==h,p
 d=pickle.loads((OUT/'joint_moments.pkl').read_bytes())['data'];norm=json.loads((OUT/'source_norms.json').read_text())
 pp=coefficient_polynomial(WITNESS);k,z=pp.gens;y=sy.Symbol('y');p=sy.Poly(sy.expand(pp.as_expr().subs(k,y+10)),y)
 cs=[max(bernstein_univariate(p.nth(j),z,0,sy.Rational(3,8)))for j in range(p.degree()+1)]
 assert all(v<=0 for v in cs)and cs[5]<0;power5=-F(cs[5]);M=max([F(0)]+[F(max(bernstein_univariate(pp.as_expr().subs(k,j),z,0,sy.Rational(3,8))))for j in range(10)])
 # For |t|<=1 and0<=w<=1/45, |omega|<5/4 and |Q^3 F|<64.
 assert ((1+F(5,45))/(1+F(3,90)))**2<F(5,4)
 assert 50*((1+F(8,135))/(1+F(3,90)))**3<64
 co,base=forward_bounds(d,omega=F(5,4),complexFF=F(64));tau=F(256);K=20
 coeff=leg_coeffs(L);ak=coeff[K]*C**K
 forward=sum((v*tau**(j-K)/ak for j,v in enumerate(co)),F(0))
 aa=ar(C*tau);m=poly([j*v for j,v in enumerate(coeff)],aa)/poly(coeff,aa);assert m>10
 R1=(1+F(2,45))*(1+F(8,135))**2/((1-F(44,45))*(1-F(157,2700)))
 R2=1/(1-(F(317,60)/F(143,3))**2)
 source1=F(norm['norm1_rational'])*R1/(ak*tau**K);source2=F(norm['norm2_rational'])*R2
 value=ar(1+base+forward+source1+source2+M)-ar(power5)*(m-10)**5
 assert value.is_finite()and value.upper()<-1
 hashes=dict(prep['sha256_manifest'])
 for path in [Path(__file__),Path(forward_bounds.__code__.co_filename),ROOT/'research/crossing_kernels.py']:hashes[str(path.relative_to(ROOT))]=sha(path)
 out=dict(status='PASS_UNBOUNDED_TAU',witness_sha256=sha(WITNESS),lmin=L,mu_domain=['49','infinity'],tau_domain=['256','infinity'],kappa_exact='1',precision_bits=a.precision,
  V_over_p_upper_exact=str(up(value)),mean_lower_exact=str(lo(m)),positive_Legendre_term_degree=K,
  forward_ratio_upper_exact=str(forward),source1_ratio_upper_exact=str(source1),source2_upper_exact=str(source2),crossing_low_k_cap_exact=str(M),crossing_power5_exact=str(power5),
  proof='Positive Legendre coefficient stochastic ordering gives m_l(tau)>=m_L(256). Fifth-power Jensen crossing coercivity dominates source/forward terms. Positive term k20 bounds every tau^j/P, j<=3, uniformly for tau>=256. All finite mu>=49 are retained.',
  review_status='New analytic envelope; independent mathematical review requested.',seconds=time.time()-tic,sha256_manifest=hashes)
 (OUT/f'large_tau_L{L}_{a.precision}.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['sha256_manifest','forward_ratio_upper_exact','source1_ratio_upper_exact','source2_upper_exact']},indent=2))
if __name__=='__main__':main()
