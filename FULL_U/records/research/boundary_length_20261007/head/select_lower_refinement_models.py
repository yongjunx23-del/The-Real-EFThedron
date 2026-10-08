"""Select existing cap methods; persist paired checks and exact sign endpoints."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, sys, time
from certify_head_constant import FiniteConstantCentered, driver
from prepare_lower_refinement import CONFIGS

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--campaign',choices=CONFIGS,required=True); a=ap.parse_args()
    root=driver.ROOT; base=root/f'output/boundary_length_20261007/head/angle{a.campaign}'
    witness=root/f'output/boundary_length_20261007/search/angle{a.campaign}.npz'
    assert sha(witness)==CONFIGS[a.campaign]['sha']
    assert not (base/'model_selection.json').exists()
    preflight=json.loads((base/'preflight.json').read_text())
    failures=[r['spin'] for r in preflight['rows'] if not r['status'].startswith('PASS')]
    assert failures in [[],[2]], failures
    budget=json.loads((base/'budgets.json').read_text()); tic=time.time()
    dependencies=driver.code_hashes(); dependencies[str(Path(__file__).relative_to(root))]=sha(Path(__file__))
    components={'preflight.json':sha(base/'preflight.json')}
    if failures:
        cap_rows=[]; endpoint_rows=[]
        for bits in [256,384]:
            driver.ctx.prec=bits; model=FiniteConstantCentered(witness,2,base/'moments128',N=128)
            gaps=model.gaps(30); bounds=[model.gap_bound(lo,hi) for lo,hi in gaps]
            assert all(v.is_finite() and v>=0 for v in bounds)
            cap_rows.append(dict(bits=bits,gaps=[[str(lo),str(hi)] for lo,hi in gaps],threshold=model.cap_status,
                gap_upper=[str(v.upper()) for v in bounds]))
            observed=[]
            def observe(frame,event,arg):
                if event=='return' and frame.f_code is driver.FiniteSpin.threshold_cap.__code__ and arg is not None:
                    vn=frame.f_locals['vn']; assert vn.upper()<0
                    observed.append(dict(lower_dyadic=list(map(str,vn.lower().man_exp())),upper_dyadic=list(map(str,vn.upper().man_exp())),
                        lower_float=float(vn.lower()),upper_float=float(vn.upper()),display_string=str(vn)))
            sys.setprofile(observe)
            try: cost=model.threshold_cap(gaps[0][1])
            finally: sys.setprofile(None)
            assert len(observed)==1
            endpoint_rows.append(dict(bits=bits,beta_domain_exact=[str(gaps[0][0]),str(gaps[0][1])],
                V_over_p=observed[0],cap_upper_dyadic=list(map(str,cost.upper().man_exp()))))
        assert cap_rows[0]['gaps']==cap_rows[1]['gaps']
        shared=dict(spin=2,witness_sha256=sha(witness),source_manifest_sha256=sha(base/'moments128/manifest.json'),dependency_sha256=dependencies)
        cap=dict(status='PASS_CONSTANT_CAP_PREFLIGHT',**shared,rows=cap_rows,seconds=time.time()-tic,
            scope='Only finite-width gaps with the unchanged reviewed constant cap; not the full head.')
        exact=dict(status='PASS_EXACT_NEGATIVE_CONSTANT_CAP_ENDPOINTS',**shared,records=endpoint_rows,
            scope='Exact directed endpoints observed from the unchanged threshold-cap return.')
        for name,data in [('constant_cap_preflight.json',cap),('constant_cap_exact_sign.json',exact)]:
            assert not (base/name).exists(); (base/name).write_text(json.dumps(data,indent=2)); components[name]=sha(base/name)
        for row in budget['rows']:
            if row['spin']==2: row['model']='constant_centered'
    (base/'budgets_selected.json').write_text(json.dumps(budget,indent=2))
    selection=dict(status='PASS_SELECTED_CAP_PREFLIGHT',witness_sha256=sha(witness),original_failed_spins=failures,
        constant_cap_spins=failures,components_sha256=components,budgets_selected_sha256=sha(base/'budgets_selected.json'),
        dependencies_sha256=dependencies,seconds=time.time()-tic,
        scope='All preflight-listed spins pass their selected existing method. No kernel changes or complete bound claimed.')
    (base/'model_selection.json').write_text(json.dumps(selection,indent=2)); print(json.dumps(selection,indent=2))

if __name__=='__main__': main()
