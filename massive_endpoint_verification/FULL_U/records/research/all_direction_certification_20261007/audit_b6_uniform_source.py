"""Independent exact finite-spin specialization of the b6 bivariate polynomials."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,math,hashlib,time
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent/'tail'))
import b6_uniform_uv as b
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def specialize(rows,eta):
    out=b.fmpq_poly([])
    for row in reversed(rows):out=out*b.fq(eta)+row
    return out
def main():
    tic=time.time();b.ctx.prec=512;v=b.B6Uniform();records=[]
    x=b.fmpq_poly([0,1]);pp=b.coefficient_polynomial(b.W);k,z=pp.gens
    ratio=(v.f.lam['X4']+F(107,5))/v.f.lam['X4'];pp=b.sy.Poly(pp.as_expr()*b.sy.Rational(ratio),k,z)
    for ell in [42,44,80,450]:
        T=ell*(ell+1);eta=F(1,T)
        aa=[F(math.comb(ell,j)*math.comb(ell+j,j),T**j)for j in range(41)]
        source=b.fmpq_poly([])
        for j,row in enumerate(v.src['moments']):
            source+=b.fmpq_poly([b.fq(F(0))]*j+[b.fq((-1)**m*aa[m]*co/T**j)for m,co in enumerate(row)])
        assert source==specialize(v.source,eta)
        q=1+b.fq(F(8,3*T))*x;cross=b.fmpq_poly([])
        for m in range(41):
            c=b.sy.Poly(pp.as_expr().subs(k,m),z)
            for j in range(max(0,c.degree())+1):
                cross+=b.fq(aa[m]*F(4,3)**m*F(c.nth(j))/T**j)*x**(m+j)*q**(v.D-j)
        assert cross==specialize(v.cross,eta)
        for tau in map(F,['0','1','2','4','8']):
            exactP=sum((F(math.comb(ell,j)*math.comb(ell+j,j))*F(4,3*T)**j*tau**j for j in range(ell+1)),F(0))
            bound=v.p_bounds(tau,tau,eta,eta)
            lower,upper=b.bounds(bound);assert lower<=exactP<=upper
            field=v.field(tau,tau,eta,eta);assert field.is_finite()
            records.append(dict(spin=ell,tau_exact=str(tau),P_exact_contained=True,field_upper_exact=str(b.upper(field))))
    out=dict(status='PASS_EXACT_FINITE_SPIN_SPECIALIZATION_B6',base_witness_sha256=sha(b.W),objective_exact=['1','6'],delta_X4_exact='107/5',
        precision_bits=512,source_and_crossing_polynomial_equalities_spins=[42,44,80,450],positive_Legendre_denominator_checks=records,
        domain_review='Every retained k<=40 coefficient is nonnegative for ell>=42. Omitted crossing k>=41 has exactly proved nonpositive coefficients; source and denominator tails use absolute positive Bessel-series majorants.',
        scope='Independent algebra and point checks support the construction; continuum negativity is established by the separate rational eta/tau covers.',
        seconds=time.time()-tic,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),Path(b.__file__),b.W]})
    p=ROOT/'output/competitive_20261007/upper_arc/b6/uniform_source_independent_audit.json';p.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items()if k not in ['sha256_manifest','positive_Legendre_denominator_checks']},indent=2))
if __name__=='__main__':main()
