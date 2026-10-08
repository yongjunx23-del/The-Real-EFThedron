"""Independent saved-ledger accounting adapted from the lower peer review; no interval kernels rerun."""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left,bisect_right
import json,hashlib
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'output/six_targets_20261007/upper95_direct';hashes={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    d=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
    for rel,h in d.get('sha256_manifest',{}).items():assert sha(ROOT/rel)==h;hashes[rel]=h
    return d
old=ROOT/'research/weak_regions_20261007/banded_fresh269p5.py';new=ROOT/'research/six_targets_20261007/banded_upper95.py'
assert old.read_text().replace('output/weak_regions_20261007/fresh269p5','output/six_targets_20261007/upper95_direct').replace("choices=['269p5']","choices=['upper95']").replace('complete_compact_coercivity_stronger.json','complete_compact_coercivity.json')==new.read_text()
hashes[str(old.relative_to(ROOT))]=sha(old);hashes[str(new.relative_to(ROOT))]=sha(new)
u=read(BASE/'tail/uv_complete.json');assert u['status']=='PASS_COMPLETE_DIRECTIONAL_UV_COERCIVITY'
source=read(BASE/'tail/exact_rebuild_and_algebra_audit.json');assert source['status']=='PASS_UV_EXACT_REBUILD_AND_ALGEBRA'and source['all_entries_equal']and source['moment_count']==972 and source['remainder_coefficient_count']==9
assert source['source_polynomial_equality_spins']==source['crossing_polynomial_equality_spins']==[450,452,600]
assert all(r['P_encloses_exact']for r in source['checks'])
compact=read(BASE/'compact/complete_compact_coercivity.json');assert F(compact['kappa_exact'])==1
finite=read(BASE/'tail/finite_uv_complete.json');assert finite['spins']==list(range(42,450,2))
uv={};paired_uv_cells=0
for pair in finite['paired_files']:
    p=read(ROOT/pair['producer']);r=read(ROOT/pair['replay'])
    assert p['witness_sha256']==r['witness_sha256']==u['witness_sha256']and p['objective_exact']==r['objective_exact']==u['objective_exact']
    assert p['precision_bits']==400 and r['precision_bits']==480 and r['reference_certificate_sha256']==sha(ROOT/pair['producer'])
    assert p['spins']==r['spins']and len(p['records'])==len(r['records'])==len(p['spins'])
    for ell,a,b in zip(p['spins'],p['records'],r['records']):
        assert ell==a['spin']==b['spin']and ell not in uv
        cells=[];cursor=F(0);assert len(a['cells'])==len(b['cells'])
        for c,d in zip(a['cells'],b['cells']):
            lo,hi=F(c['a']),F(c['b']);assert lo==F(d['a'])==cursor<hi==F(d['b'])
            cap=max(F(c['upper_exact']),F(d['upper_exact']));assert cap<-F(1,200);cursor=hi;cells.append((lo,hi,-cap));paired_uv_cells+=1
        assert cursor==min(F(256),F(ell*(ell+1),45));uv[ell]=cells
assert sorted(uv)==list(range(42,450,2))
starts={ell:[v[0]for v in cells]for ell,cells in uv.items()};ends={ell:[v[1]for v in cells]for ell,cells in uv.items()}
real=read(BASE/'real/finite_split_replay.json');cost=read(BASE/'real/finite_banded_cost.json')
assert cost['witness_sha256']==real['witness_sha256']==u['witness_sha256']and cost['objective_exact']==u['objective_exact']
assert len(real['cells'])==len(cost['cells']);total=F(0);cursor=F(0)
for a,c in zip(real['cells'],cost['cells']):
    lo,hi=F(a['lo']),F(a['hi']);assert lo==F(c['lo'])==cursor<hi==F(c['hi']);cursor=hi
    mass=max(F(a['first_upper_exact']),F(a['replay_upper_exact']));assert mass==F(a['retained_upper_exact'])==F(c['A_upper_exact'])>=0
    margins=[];expected_ranges=[];expected_large=[]
    if lo*lo<F(45,49):margins.append(F(1))
    for ell in range(42,450,2):
        if hi*hi<=F(45,49):continue
        T=ell*(ell+1);ta=F(T)*(1-hi*hi)/(4*hi*hi)
        tb=min(F(T,45),F(T)*(1-lo*lo)/(4*lo*lo))if lo else F(T,45)
        if tb>256:assert ell>=108;expected_large.append(ell);margins.append(F(1))
        tb=min(tb,F(256))
        if ta<=tb:
            touched=list(range(bisect_left(ends[ell],ta),bisect_right(starts[ell],tb)))
            assert touched and touched==list(range(touched[0],touched[-1]+1))
            assert uv[ell][touched[0]][0]<=ta<=uv[ell][touched[0]][1]and uv[ell][touched[-1]][1]>=tb
            expected_ranges.append([ell,touched[0],touched[-1]+1]);margins.extend(uv[ell][i][2]for i in touched)
    assert expected_ranges==c['finite_tau_cell_index_ranges']and expected_large==c['large_tau_spins']
    k=min(margins);assert k==F(c['kappa_exact'])>0
    w=F(c['reciprocal_weight_upper_exact']);assert w>=1/(2*k)
    contribution=mass*w;assert contribution==F(c['cost_upper_exact']);total+=contribution
assert cursor==1 and total==F(cost['cost_upper_exact'])
for name in ['assemble_uv_peak.py','review_upper_tail.py']:
    p=Path(__file__).with_name(name);hashes[str(p.relative_to(ROOT))]=sha(p)
report=dict(status='PASS_UPPER_TAIL_INDEPENDENT_COVERAGE_AND_ACCOUNTING',witness_sha256=u['witness_sha256'],paired_finite_uv_cells=paired_uv_cells,real_cells_checked=len(cost['cells']),spins_per_real_cell=204,cost_upper_exact=str(total),cost_upper_float=float(total),sha256_manifest=hashes,review='Generic UV assembly checked for exact finite42..448 and eta0..1/(450*451) infinite coverage, paired400/480 partitions, tau256 large-region closure and witness identities. Banded adaptation changes paths/case and uses the original compact kappa1 packet. This separate ledger check reconstructs every touching UV interval for every real cell and spin without importing the banded coverage helper, then sums exact outward weighted costs.',scope='Internal exact accounting check using an independently constructed closed-interval index; physical kernels and interval evaluators share the established implementation. This is not a complete directional certificate without head and infinite-spin real costs.')
(BASE/'tail/independent_coverage_accounting.json').write_text(json.dumps(report,indent=2));print(report['status'],paired_uv_cells,'pairedUVcells',len(cost['cells']),'realcells',float(total))
