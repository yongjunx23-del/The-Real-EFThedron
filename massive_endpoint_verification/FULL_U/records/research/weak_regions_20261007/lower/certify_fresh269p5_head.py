"""Finite head only for the new269.5 witness, using unchanged proof kernels."""
import os,sys,json,time,hashlib,subprocess,concurrent.futures
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/weak_regions_20261007/lower'
BASE=OUT/'head_certification/angle269p5'
WITNESS=OUT/'search/angle269p5.npz'
EXPECTED='f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae'
OLD=ROOT/'research/boundary_length_20261007/head'
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(WITNESS)==EXPECTED
assert not BASE.exists(),'Preserve existing proof campaign; use a fresh path.'
BASE.mkdir(parents=True)
commands=[];start=time.time()
def run(command,log):
    tic=time.time()
    with log.open('w') as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
    entry=dict(command=command,log=str(log.relative_to(ROOT)),exit_code=r.returncode,elapsed_seconds=time.time()-tic)
    assert r.returncode==0,entry
    return entry
for label in ('moments128','moments128_replay'):
    command=[sys.executable,str(ROOT/'research/all_direction_certification_20261007/head/build_moments.py'),'--witness',str(WITNESS),'--out',str(BASE/label),'--lmax','40','--N','128','--Nuv','32','--bits','256']
    commands.append(run(command,BASE/(label+'.log')))
    print(label,'complete',flush=True)
manifest=json.loads((BASE/'moments128/manifest.json').read_text());checks=[]
for t in manifest['tables']:
    for suffix,key in [('_moments.json','moment_sha256'),('_uv_moments.json','uv_sha256')]:
        name=t['label']+suffix
        assert sha(BASE/'moments128'/name)==sha(BASE/'moments128_replay'/name)==t[key]
        checks.append(dict(file=name,sha256=t[key],bit_for_bit_rebuild=True))
sys.path.insert(0,str(OLD))
from certify_head_constant import FiniteConstantCentered,driver
numeric=json.loads((OUT/'head384/angle269p5.json').read_text());assert numeric['sha256']==EXPECTED
rows=[];preflight=[]
for ell in range(0,41,2):
    kind='centered' if ell<4 else 'fast'
    if ell==2:kind='constant_centered'
    slack={0:'.000020',2:'.000006',4:'.000005'}.get(ell,'.000001')
    rows.append(dict(spin=ell,model=kind,numerical=str(numeric['per_spin'][str(ell)]),slack=slack,target=str(Decimal(str(numeric['per_spin'][str(ell)]))+Decimal(slack))))
    if ell not in [0,2,4,6,8,14,40]:continue
    for bits in [256,384]:
        driver.ctx.prec=bits;model=driver.MODELS[kind](WITNESS,ell,BASE/'moments128',N=128)
        gaps=model.gaps(30);bounds=[model.gap_bound(lo,hi) for lo,hi in gaps]
        assert all(x.is_finite() and x>=0 for x in bounds)
        preflight.append(dict(spin=ell,model=kind,precision_bits=bits,gaps=[[str(lo),str(hi)]for lo,hi in gaps],gap_upper_dyadic=[list(map(str,x.upper().man_exp())) for x in bounds],threshold=model.cap_status))
budget=dict(witness_sha256=EXPECTED,rows=rows,slack_sum_exact=str(sum((Fraction(r['slack']) for r in rows),Fraction(0))),scope='Numerical values only choose quality targets; validity uses actual directed enclosures.')
(BASE/'budgets.json').write_text(json.dumps(budget,indent=2))
(BASE/'preparation.json').write_text(json.dumps(dict(status='PASS_PREPARED_FINITE_HEAD_INPUTS',witness_sha256=EXPECTED,source_rebuild_checks=checks,preflight=preflight,commands=commands,script_sha256=sha(Path(__file__)),dependencies_sha256=driver.code_hashes(),elapsed_seconds=time.time()-start),indent=2))
print('PREPARED',time.time()-start,flush=True)
def work(row):
    ell=row['spin'];folder=BASE/f'spin{ell}';folder.mkdir()
    command=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'head_l{ell}.json'),'--spin',str(ell),'--model',row['model'],'--target',row['target'],'--max-cells','7000','--bits','256']
    runs=[run(command,folder/f'head_l{ell}.log')]
    command=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'replay_l{ell}.json'),'--replay',str(folder/f'head_l{ell}.json'),'--bits','384']
    runs.append(run(command,folder/f'replay_l{ell}.log'))
    (folder/'commands.json').write_text(json.dumps(runs,indent=2));print('REPLAYED',ell,flush=True)
    return dict(spin=ell,commands=runs)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:jobs=list(pool.map(work,rows))
(BASE/'jobs.json').write_text(json.dumps(jobs,indent=2))
command=[sys.executable,str(OLD/'assemble.py'),'--witness',str(WITNESS),'--base',str(BASE)]
commands.append(run(command,BASE/'assembly.log'))
(BASE/'campaign.json').write_text(json.dumps(dict(status='COMPLETE_REPLAYED_FINITE_HEAD_ONLY',witness_sha256=EXPECTED,jobs=jobs,commands=commands,elapsed_seconds=time.time()-start,scope='Every energy, even spins0..40 only. Same-witness spin42..infinity tail is still required.'),indent=2))
print((BASE/'head_summary.json').read_text(),flush=True)
