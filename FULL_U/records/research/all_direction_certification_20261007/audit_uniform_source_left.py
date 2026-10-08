"""Independent exact specialization of the bivariate left UV construction."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent/'tail'))
import uniform_source_uv as u
import finite_uv_bridge220 as f

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def specialize(rows,eta):
    out=u.fmpq_poly([])
    for p in reversed(rows):out=out*u.fq(eta)+p
    return out

def main():
    u.ctx.prec=512;tic=time.time();b=u.UniformSource(450)
    f.WITNESS=u.LEFT;f.OUT=u.FOLDER;c=f.FiniteUV();records=[]
    for ell in [450,452,600]:
        T,pp,sp,cn=c.prepare(ell);eta=F(1,T)
        assert specialize(b.source,eta)==sp
        assert specialize(b.cross,eta)==cn
        for tau in map(F,['0','1','9/4','8','64','256']):
            # FiniteUV has a hardcoded g2 objective +1; the left objective is0.
            direct=c.field(ell,tau,tau)-1
            extended=b.field(tau,tau,eta,eta)
            assert direct.is_finite() and extended.is_finite()
            assert direct.overlaps(extended),(ell,tau,str(direct),str(extended))
            P=b.p_bounds(tau,tau,eta,eta)
            exactP=u.arb(pp(u.fq(tau)))
            assert P.overlaps(exactP)
            records.append(dict(spin=ell,tau_exact=str(tau),finite_enclosure=str(direct),uniform_enclosure=str(extended),P_encloses_exact=True))
    out=dict(status='PASS_EXACT_ETA_SPECIALIZATION_AND_INTERVAL_COMPARISON',witness_sha256=sha(u.LEFT),
        source_polynomial_equality_spins=[450,452,600],crossing_polynomial_equality_spins=[450,452,600],
        precision_bits=512,checks=records,seconds=time.time()-tic,
        scope='Algebra and point-enclosure audit only. Complete continuum negativity requires separate tau/eta cell certificates.',
        sha256_manifest={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),Path(u.__file__),Path(f.__file__),u.LEFT]})
    path=ROOT/'output/competitive_20261007/left/uniform_source_independent_audit.json';path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['checks','sha256_manifest']},indent=2))
if __name__=='__main__':main()
