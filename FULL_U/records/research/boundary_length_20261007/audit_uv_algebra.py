"""Fresh exact cache rebuild and independent finite/uniform specialization.

This is an algebra audit, not the all-energy negativity certificate.
"""
from pathlib import Path
from fractions import Fraction as F
import sys, argparse, importlib, pickle, json, time, hashlib
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research/boundary_length_20261007/tail'))


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--angle',choices=['85','255','265'],required=True);a=ap.parse_args()
    c=importlib.import_module(f'angle{a.angle}_campaign')
    models=importlib.import_module(f'angle{a.angle}_models')
    import finite_source_corrections as src
    import uniform_source_uv as u
    import finite_uv_bridge220 as f
    tic=time.time();manifest=c.source_manifest()
    source=c.OUT/'finite_source_N12_K80.pkl'
    old=pickle.loads(source.read_bytes())
    src.WITNESS=c.WITNESS;src.OUT=c.OUT
    fresh=src.build(old['N'],old['K'])
    keys=['witness_sha256','N','K','moments','remainder_coefficients','remainder_by_source','definition','denominator_lower','proof','builder_sha256']
    for key in keys:assert old[key]==fresh[key],key
    assert old['witness_sha256']==c.SHA
    u.ctx.prec=512;b=models.Uniform(450);finite=models.Finite();records=[]
    def specialize(rows,eta):
        out=u.fmpq_poly([])
        for row in reversed(rows):out=out*u.fq(eta)+row
        return out
    for ell in [450,452,600]:
        T,pp,sp,cn=finite.prepare(ell);eta=F(1,T)
        assert specialize(b.source,eta)==sp
        assert specialize(b.cross,eta)==cn
        for tau in map(F,['0','1','9/4','8','64','256']):
            direct=finite.field(ell,tau,tau)
            extended=b.field(tau,tau,eta,eta)
            assert direct.is_finite() and extended.is_finite() and direct.overlaps(extended)
            P=b.p_bounds(tau,tau,eta,eta);exactP=u.arb(pp(u.fq(tau)))
            assert P.overlaps(exactP)
            records.append(dict(spin=ell,tau_exact=str(tau),finite_enclosure=str(direct),uniform_enclosure=str(extended),P_encloses_exact=True))
    files=[Path(__file__),Path(src.__file__),Path(u.__file__),Path(f.__file__),Path(models.__file__),c.WITNESS,source]
    manifest.update({str(p.relative_to(ROOT)):sha(p)for p in files})
    result=dict(status='PASS_UV_EXACT_REBUILD_AND_ALGEBRA',witness_sha256=c.SHA,objective_exact=list(map(str,c.OBJECTIVE)),
        all_entries_equal=True,verified_keys=keys,moment_count=len(old['moments'])*len(old['moments'][0]),
        remainder_coefficient_count=sum(len(g['coefficient_abs_upper_exact'])for g in old['remainder_by_source']),
        source_polynomial_equality_spins=[450,452,600],crossing_polynomial_equality_spins=[450,452,600],
        precision_bits=512,checks=records,sha256_manifest=manifest,seconds=time.time()-tic,
        scope='Fresh exact-rational source cache rebuild using the reviewed builder, plus independent algebraic specialization of uniform polynomials to finite-spin polynomials and 18 point enclosure comparisons. Does not replace the separate complete interval cover.')
    path=c.OUT/'exact_rebuild_and_algebra_audit.json';assert not path.exists()
    path.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items()if k not in ['checks','sha256_manifest']},indent=2))


if __name__=='__main__':main()
