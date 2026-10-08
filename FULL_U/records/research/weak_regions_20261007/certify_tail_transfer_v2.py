"""Fail-closed transfer of completed all-spin tail proof by negative g2 shift.

The source witness, every old proof dependency and all finite-tail coverage
are checked. No quadrature, source moment or nonpositive-tail assumption is
introduced by this algebraic transfer. No finite head is claimed here.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,os,subprocess,sys,time
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research/boundary_length_20261007'))
from assemble_complete_v3 import Inputs,sha,dyadic,ceiling_decimal
from banded_finite_cost import dyadic_ceiling
BASE=ROOT/'output/boundary_length_20261007'
OUT=ROOT/'output/weak_regions_20261007/tail_transfer'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--tag',default='from268_to269p04');a=ap.parse_args();started=time.time()
    p=OUT/'definitions_v2.json';probes=json.loads(p.read_text());case=next(c for c in probes['cases']if c['tag']==a.tag)
    alias=case['base'];assert alias in ['85','261p3814','268'];wp=ROOT/case['witness'];oldwp=ROOT/case['base_witness']
    assert sha(wp)==case['witness_sha256'] and sha(oldwp)==case['base_sha256']
    old=np.load(oldwp);new=np.load(wp)
    assert len(new['names'])==len(new['multipliers'])==220
    assert np.array_equal(new['names'],old['names'])
    assert np.array_equal(new['multipliers'],old['multipliers'])
    oldd=list(map(lambda x:F(float(x)),old['direction']));newd=list(map(lambda x:F(float(x)),new['direction']))
    assert newd[0]==oldd[0] and newd[1]<oldd[1];delta=newd[1]-oldd[1]
    assert [str(x)for x in newd]==case['direction_exact'] and str(delta)==case['g2_shift_exact']
    target=OUT/(a.tag+'_tail_certificate.json');assert not target.exists()
    recert=OUT/('base'+alias+'_reassembled_for_transfer.json')
    cmd=[sys.executable,str(ROOT/'research/boundary_length_20261007/assemble_complete_v3.py'),'--angle',alias,
         '--finite-banded',str(BASE/f'real/angle{alias}/finite_banded_cost.json'),'--out',str(recert)]
    if alias=='268':cmd+=['--head',str(BASE/'head/angle268/head_summary_final.json')]
    if not recert.exists():
        with (OUT/('base'+alias+'_reassembled_for_transfer.log')).open('w')as fh:
            subprocess.run(cmd,cwd=ROOT,stdout=fh,stderr=subprocess.STDOUT,check=True)
    inputs=Inputs();base=inputs.read(recert);assert base['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
    assert base['witness_sha256']==case['base_sha256'] and base['objective_exact']==[str(x)for x in oldd]
    assert F(base['support_upper_exact'])==sum(map(F,base['components_exact'].values()),F(0))
    inputs.bind(wp);inputs.bind(oldwp);inputs.bind(Path(__file__))
    for key,digest in base['sha256_manifest'].items():inputs.bind(key,digest)
    fin=inputs.read(BASE/f'real/angle{alias}/finite_banded_cost.json');inf=inputs.read(BASE/f'real/angle{alias}/tail_split_replay.json')
    assert fin['witness_sha256']==inf['witness_sha256']==case['base_sha256']
    assert fin['status']=='PASS_BANDED_COERCIVITY_REAL_COST' and inf['status']=='PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY'
    assert (fin['first_even_spin'],fin['last_even_spin'])==(42,448)
    assert (inf['first_even_spin'],inf['last_even_spin'])==(450,'infinity')
    assert fin['energy_domain']==['4','infinity']
    ledger=[];fcost=F(0)
    for cell in fin['cells']:
        oldk=F(cell['kappa_exact']);k=oldk-delta;assert k>oldk>0
        A=F(cell['A_upper_exact']);w=dyadic_ceiling(1/(2*k));assert w>=1/(2*k)
        cost=A*w;fcost+=cost
        ledger.append(dict(lo=cell['lo'],hi=cell['hi'],region=cell['region'],A_upper_exact=str(A),
            base_kappa_exact=str(oldk),transferred_kappa_exact=str(k),reciprocal_weight_upper_exact=str(w),cost_upper_exact=str(cost)))
    assert ledger[0]['lo']=='0' and ledger[-1]['hi']=='1'
    assert all(F(x['hi'])==F(y['lo'])for x,y in zip(ledger,ledger[1:]))
    constants=base['coercivity_constants_exact'];ck=F(constants['compact']);ik=F(constants['infinite_UV'])
    components={'finite_42_448_all_energy':fcost};tail=[]
    for region,k in [('compact',ck),('remainder',min(ck,ik))]:
        A=F(inf['A_upper_exact'][region]);newk=k-delta;w=dyadic_ceiling(1/(2*newk));cost=A*w
        components['infinite_450_up_'+region]=cost
        tail.append(dict(region=region,A_upper_exact=str(A),base_kappa_exact=str(k),transferred_kappa_exact=str(newk),
            reciprocal_weight_upper_exact=str(w),cost_upper_exact=str(cost)))
    total=sum(components.values(),F(0));assert total<=F(base['tail_upper_exact'])
    result=dict(status='PASS_COMPLETE_TRANSFERRED_CONTINUUM_TAIL',witness=str(wp.relative_to(ROOT)),witness_sha256=sha(wp),
        objective_exact=[str(x)for x in newd],base_witness=str(oldwp.relative_to(ROOT)),base_witness_sha256=sha(oldwp),
        base_completed_certificate=str(recert.relative_to(ROOT)),base_completed_certificate_sha256=sha(recert),
        g2_shift_exact=str(delta),multipliers_bitwise_identical=True,g0_objective_exactly_identical=True,
        field_identities=['U_new=U_base','V_new=V_base+delta*p','p=P_ell(1+8/(3*(mu-4)))/(mu-4/3)^3>0'],
        coercivity_identity='V_base<=-kappa*p implies V_new<=-(kappa-delta)*p, delta<0',
        support_inequality='sqrt(U^2+V^2)+V<=U^2/(2*(kappa-delta)*p)',
        first_even_spin=42,last_even_spin='infinity',energy_domain=['4','infinity'],lower_energy_endpoint_open=True,
        finite_tail_cells=ledger,infinite_tail_regions=tail,components_exact={k:str(v)for k,v in components.items()},
        components_float={k:float(v)for k,v in components.items()},tail_upper_exact=str(total),tail_upper_float=float(total),
        source_reassembly_command=cmd,sha256_manifest=inputs.manifest,elapsed_seconds=time.time()-started,
        scope='Every even spin>=42, every mu>4. The new witness finite head0..40 must be separately certified; this tail packet alone is not a coupling bound.')
    target.write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k]for k in ['status','tail_upper_float','components_float','elapsed_seconds']},indent=2),flush=True)
if __name__=='__main__':main()
