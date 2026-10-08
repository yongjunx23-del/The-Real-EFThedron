"""Freeze exact repaired83 head inputs; existing mathematical kernels unchanged."""
import os,sys,json,time,hashlib,subprocess
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/peak_certification_20261007/upper83/head'
WITNESS=ROOT/'output/weak_regions_20261007/upper/angle83_repaired.npz'
EXPECTED='68e39591f3f9a9d3b5cdfe3ade3924fecb227a329de4a9543729f0607e8f663a'
OLD=ROOT/'research/boundary_length_20261007/head'
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert sha(WITNESS)==EXPECTED
    assert not BASE.exists(),'Preserve existing proof campaign; use a fresh path.'
    BASE.mkdir(parents=True);commands=[];start=time.time()
    for label in ('moments128','moments128_replay'):
        command=[sys.executable,str(ROOT/'research/all_direction_certification_20261007/head/build_moments.py'),'--witness',str(WITNESS),'--out',str(BASE/label),'--lmax','40','--N','128','--Nuv','32','--bits','256']
        tic=time.time();log=BASE/(label+'.log')
        with log.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
        commands.append(dict(command=command,log=str(log.relative_to(ROOT)),exit_code=r.returncode,elapsed_seconds=time.time()-tic))
        assert r.returncode==0,commands[-1]
        print(label,'complete',time.time()-tic,flush=True)
    manifest=json.loads((BASE/'moments128/manifest.json').read_text());checks=[]
    for t in manifest['tables']:
        for suffix,key in [('_moments.json','moment_sha256'),('_uv_moments.json','uv_sha256')]:
            name=t['label']+suffix
            assert sha(BASE/'moments128'/name)==sha(BASE/'moments128_replay'/name)==t[key]
            checks.append(dict(file=name,sha256=t[key],bit_for_bit_rebuild=True))
    sys.path.insert(0,str(OLD))
    from certify_head_constant import FiniteConstantCentered,driver
    numeric_path=ROOT/'output/weak_regions_20261007/upper/repaired83_audit/head384/angle83_repaired.json'
    numeric=json.loads(numeric_path.read_text());assert numeric['sha256']==EXPECTED
    rows=[];preflight=[]
    for ell in range(0,41,2):
        kind='centered'if ell<4 else'fast'
        if ell in [0,2,4,6,8,14,40]:
            for bits in [256,384]:
                driver.ctx.prec=bits;model=driver.MODELS[kind](WITNESS,ell,BASE/'moments128',N=128)
                try:gaps=model.gaps(30);bounds=[model.gap_bound(lo,hi)for lo,hi in gaps]
                except ValueError:
                    if ell!=2:raise
                    kind='constant_centered';model=driver.MODELS[kind](WITNESS,ell,BASE/'moments128',N=128)
                    gaps=model.gaps(30);bounds=[model.gap_bound(lo,hi)for lo,hi in gaps]
                assert all(x.is_finite()and x>=0 for x in bounds)
                preflight.append(dict(spin=ell,model=kind,precision_bits=bits,gaps=[[str(lo),str(hi)]for lo,hi in gaps],gap_upper_dyadic=[list(map(str,x.upper().man_exp()))for x in bounds],threshold=model.cap_status))
                print('CAP_PASS',ell,bits,kind,flush=True)
        slack={0:'.00025',2:'.00008',4:'.00003',6:'.00002'}.get(ell,'.00001')
        rows.append(dict(spin=ell,model=kind,numerical=str(numeric['per_spin'][str(ell)]),slack=slack,target=str(Decimal(str(numeric['per_spin'][str(ell)]))+Decimal(slack))))
    budget=dict(witness_sha256=EXPECTED,rows=rows,slack_sum_exact=str(sum((Fraction(r['slack'])for r in rows),Fraction(0))),scope='Numerical values choose quality targets only; actual interval enclosures determine validity.')
    (BASE/'budgets.json').write_text(json.dumps(budget,indent=2))
    definition=json.loads((ROOT/'output/weak_regions_20261007/upper/angle83_repaired_exact_coefficients.json').read_text())
    (BASE/'exact_witness_definition.json').write_text(json.dumps(definition,indent=2))
    (BASE/'preparation.json').write_text(json.dumps(dict(status='PASS_PREPARED_FINITE_HEAD_INPUTS',witness_sha256=EXPECTED,source_rebuild_checks=checks,preflight=preflight,commands=commands,script_sha256=sha(__file__),dependencies_sha256=driver.code_hashes(),numerical_budget_source=str(numeric_path.relative_to(ROOT)),numerical_budget_source_sha256=sha(numeric_path),elapsed_seconds=time.time()-start),indent=2))
    print('PREPARED',time.time()-start,flush=True)
if __name__=='__main__':main()
