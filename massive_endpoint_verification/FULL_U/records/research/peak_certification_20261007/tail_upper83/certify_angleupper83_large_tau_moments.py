"""All negative shifted crossing powers and Jensen, with the exact objective.

If F(k,z) <= M-sum_j a_j (k-10)_+^j, a_j>=0 and j>=1,
convexity gives E F <= M-sum_j a_j (E k-10)_+^j.
The exact Bernstein coefficient enclosure supplies the uniform a_j.
"""
from angleupper83_campaign import *
from compact_highspin220 import ar,arb,ctx,leg_coeffs,poly,forward_bounds,coefficient_polynomial,bernstein_univariate,up,lo,C
import sympy as sy
import pickle,time,argparse

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=320);a=ap.parse_args();ctx.prec=a.precision;tic=time.time();L=108;manifest=source_manifest()
 d=pickle.loads((OUT/'joint_moments.pkl').read_bytes())['data'];norm=json.loads((OUT/'source_norms.json').read_text())
 pp=coefficient_polynomial(WITNESS);k,z=pp.gens;y=sy.Symbol('y');p=sy.Poly(sy.expand(pp.as_expr().subs(k,y+10)),y)
 cs=[max(bernstein_univariate(p.nth(j),z,0,sy.Rational(3,8)))for j in range(p.degree()+1)];assert all(v<=0 for v in cs)and cs[5]<0
 power5=-F(cs[5]);M=max([F(0)]+[F(max(bernstein_univariate(pp.as_expr().subs(k,j),z,0,sy.Rational(3,8))))for j in range(10)])
 assert ((1+F(5,45))/(1+F(3,90)))**2<F(5,4)and 50*((1+F(8,135))/(1+F(3,90)))**3<64
 co,base=forward_bounds(d,omega=F(5,4),complexFF=F(64));tau=F(256);K=20;coeff=leg_coeffs(L);ak=coeff[K]*C**K
 forward=sum((v*tau**(j-K)/ak for j,v in enumerate(co)),F(0));aa=ar(C*tau);m=poly([j*v for j,v in enumerate(coeff)],aa)/poly(coeff,aa);assert m>10
 R1=(1+F(2,45))*(1+F(8,135))**2/((1-F(44,45))*(1-F(157,2700)));R2=1/(1-(F(317,60)/F(143,3))**2)
 source1=F(norm['norm1_rational'])*R1/(ak*tau**K);source2_absolute=F(norm['norm2_rational'])*R2
 cd=pickle.loads((OUT/'compact_norms.pkl').read_bytes());assert cd['witness_sha256']==SHA and cd['moments'][0]==d['Hs']
 zm=F(3,143);rr=F(317,60)*zm;N2=12;assert rr<1 and len(cd['moments'])>=N2
 source2=-cd['moments'][0]+sum((abs(cd['moments'][k])*zm**(2*k)for k in range(1,N2)),F(0))+F(norm['norm2_rational'])*rr**(2*N2)/(1-rr*rr)
 crossing_powers=[-F(v) for v in cs]
 value=ar(OBJECTIVE[1]+base+forward+source1+source2+M)-sum((ar(v)*(m-10)**j for j,v in enumerate(crossing_powers)if j>=1),arb(0))
 print('DIAGNOSTIC',str(value), 'mean',str(m),'source2',float(source2),'base',float(base),'forward',float(forward),'cap',float(M),flush=True)
 assert value.is_finite()and value.upper()<-1
 files=[Path(__file__),Path(forward_bounds.__code__.co_filename),Path(coefficient_polynomial.__code__.co_filename),ROOT/'research/crossing_kernels.py']
 manifest.update({str(p.relative_to(ROOT)):sha(p)for p in files})
 out=dict(status='PASS_ANGLEupper83_UNBOUNDED_TAU',witness_sha256=SHA,objective_exact=list(map(str,OBJECTIVE)),lmin=L,mu_domain=['49','infinity'],tau_domain=['256','infinity'],kappa_exact='1',precision_bits=a.precision,V_over_p_upper_exact=str(up(value)),crossing_negative_shifted_powers_exact=list(map(str,crossing_powers)),jensen_powers=list(range(1,len(cs))),crossing_low_k_cap_exact=str(M),mean_lower_exact=str(lo(m)),source1_upper_exact=str(source1),source2_upper_exact=str(source2),source2_absolute_upper_exact=str(source2_absolute),source2_moment_count=N2,source2_zmax_exact=str(zm),source2_remainder_ratio_exact=str(rr),source2_moments_exact=list(map(str,cd['moments'][:N2])),source2_enclosure='-m0 plus absolute higher moments plus geometric signed-source remainder; valid for 0<=z<=3/143',forward_upper_exact=str(forward),sha256_manifest=manifest,seconds=time.time()-tic)
 (OUT/f'large_tau_{a.precision}.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in['sha256_manifest','source1_upper_exact','source2_upper_exact','forward_upper_exact']},indent=2))

if __name__=='__main__':main()
