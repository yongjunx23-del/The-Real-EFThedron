"""Save exact endpoints hidden by Arb's deliberately coarse interval string.

The profile hook observes the local vn ball at the unchanged cap method's
return. It never changes that ball or the evaluated formula.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys,time
from certify_head_constant import FiniteConstantCentered,driver

ROOT=driver.ROOT
WITNESS=ROOT/'output/boundary_length_20261007/search/angle259p2941_from255.npz'
BASE=ROOT/'output/boundary_length_20261007/head/angle259'

def main():
    out=BASE/'constant_cap_exact_sign.json';assert not out.exists()
    tic=time.time();records=[]
    for bits in [256,384]:
        driver.ctx.prec=bits;model=FiniteConstantCentered(WITNESS,2,BASE/'moments128',N=128)
        observed=[]
        def observe(frame,event,arg):
            if event=='return' and frame.f_code is driver.FiniteSpin.threshold_cap.__code__ and arg is not None:
                value=frame.f_locals['vn'];assert value.upper()<0
                observed.append(dict(lower_dyadic=list(map(str,value.lower().man_exp())),
                    upper_dyadic=list(map(str,value.upper().man_exp())),lower_float=float(value.lower()),
                    upper_float=float(value.upper()),display_string=str(value)))
        sys.setprofile(observe)
        try:cost=model.threshold_cap(F(1,256))
        finally:sys.setprofile(None)
        assert len(observed)==1
        records.append(dict(bits=bits,beta_domain_exact=['0','1/256'],V_over_p=observed[0],
                            cap_upper_dyadic=list(map(str,cost.upper().man_exp()))))
    deps=driver.code_hashes();deps[str(Path(__file__).relative_to(ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(status='PASS_EXACT_NEGATIVE_CONSTANT_CAP_ENDPOINTS',witness_sha256=hashlib.sha256(WITNESS.read_bytes()).hexdigest(),
                source_manifest_sha256=hashlib.sha256((BASE/'moments128/manifest.json').read_bytes()).hexdigest(),
                spin=2,records=records,dependency_sha256=deps,seconds=time.time()-tic,
                scope='Exact endpoints of the existing strictly negative finite-width cap; supplemental representation, not a new proof formula.')
    out.write_text(json.dumps(result,indent=2));print(json.dumps(dict(status=result['status'],records=records),indent=2))

if __name__=='__main__':main()
