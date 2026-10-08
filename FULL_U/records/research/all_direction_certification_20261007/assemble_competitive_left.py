"""Strict complete assembly for the frozen left180 witness and its exact tails."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,argparse
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research'))
from angle_benchmarks import support_report
SHA='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyadic(v):m,e=map(int,v);return F(m)*F(2)**e
def upward(v):n=-((-v.numerator*10**6)//v.denominator);return f'{n//10**6}.{n%10**6:06d}'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--finite-cost',type=Path,default=ROOT/'output/competitive_20261007/left/real/banded_cost.json')
    ap.add_argument('--head',type=Path,default=ROOT/'output/all_direction_certification_20261007/head/left/head_summary.json')
    args=ap.parse_args();manifest={};seen=set()
    def bind(p,h=None):
        p=Path(p);p=p if p.is_absolute()else ROOT/p;p=p.resolve();key=str(p.relative_to(ROOT));actual=sha(p)
        if h is not None:assert actual==h,key
        if key in manifest:assert manifest[key]==actual,key
        manifest[key]=actual
        return p,key
    def read(p,expected=None):
        p,key=bind(p,expected);d=json.loads(p.read_text())
        if key not in seen:
            seen.add(key)
            for sub,h in d.get('sha256_manifest',{}).items():
                q,_=bind(sub,h)
                if q.suffix=='.json':read(q,h)
        return d
    base=ROOT/'output/all_direction_certification_20261007/head/left';head=read(args.head.resolve())
    assert head['status']=='COMPLETE_REPLAYED_FINITE_HEAD'and head['witness_sha256']==SHA
    assert head['direction']==[-1,0]and head['included_spins']==list(range(0,41,2))
    assert [r['spin']for r in head['records']]==list(range(0,41,2))
    assert not head['missing_spins']and not head['missing_replays']
    for p,h in head['assembly_and_source_code_sha256'].items():bind(p,h)
    for r in head['source_checks']:
        assert r['bit_for_bit_rebuild_match'];p,_=bind(r['file'],r['sha256']);bind(base/'moments128_replay'/p.name,r['sha256'])
    source=read(base/'moments128/manifest.json');assert source['witness_sha256']==SHA
    for t in source['tables']:
        m=t['moment_metadata'];assert 2*m['order']-1>=max(m['energy_degree_bound'],m['angular_degree_bound'])
    first=F(0);second=F(0);retained=F(0)
    for row in head['records']:
        r=read(row['path'],row['sha256']);q=read(row['replay'],row['replay_sha256'])
        assert r['witness_sha256']==q['witness_sha256']==SHA and r['spin']==q['spin']==row['spin']
        assert q['passed']and q['partition_covers_beta_zero_to_one_exactly']
        assert q['source_record_sha256']==row['sha256']and r['source_manifest_sha256']==sha(base/'moments128/manifest.json')
        for group in [r['dependencies'],q['dependency_sha256']]:
            for p,h in group.items():bind(p,h)
        a,b=dyadic(r['upper_dyadic']),dyadic(q['upper_dyadic'])
        assert a==F(row['upper_exact'])and b==F(row['replayed_upper_exact'])
        first+=a;second+=b;retained+=max(a,b)
    assert first==F(head['upper_exact'])and second==F(head['replayed_partial_upper_exact'])
    compact=read(ROOT/'output/competitive_20261007/left/compact/complete_compact_coercivity.json')
    assert compact['status']=='PASS_COMPLETE_COMPACT_COERCIVITY'and compact['witness_sha256']==SHA
    assert compact['objective_exact']==['-1','0']and compact['minimum_even_spin']==42
    assert compact['mu_domain_exact']==['4','49']and compact['lower_energy_endpoint_open']
    ck=F(compact['kappa_exact']);assert ck>0
    uv=read(ROOT/'output/all_direction_certification_20261007/tail/left/refined/infinite_uv_complete.json')
    assert uv['status']=='PASS_LEFT_INFINITE_SPIN_UV'and uv['witness_sha256']==SHA
    assert uv['objective_exact']==['-1','0']and uv['lmin']==450 and uv['mu_domain']==['49','infinity']
    uk=F(uv['kappa_exact']);assert 0<uk<=ck
    rebuilt=read(ROOT/'output/all_direction_certification_20261007/tail/left/refined/finite_source_exact_rebuild.json')
    assert rebuilt['status']=='PASS_EXACT_RATIONAL_REBUILD'and rebuilt['witness_sha256']==SHA
    assert rebuilt['all_entries_equal']and rebuilt['moment_count']==972 and rebuilt['remainder_coefficient_count']==9
    finite=read(args.finite_cost.resolve());assert finite['witness_sha256']==SHA and finite['objective_exact']==['-1','0']
    assert finite['first_even_spin']==42 and finite['last_even_spin']==448 and finite['energy_domain']==['4','infinity']
    assert finite['status']in ['PASS_BANDED_COERCIVITY_REAL_COST','PASS_LOCAL_COERCIVITY_WEIGHTED_REAL_COST']
    fsum=F(0);cursor=F(0)
    for c in finite['cells']:
        lo,hi=F(c['lo']),F(c['hi']);assert cursor==lo<hi<=1;cursor=hi
        if finite['status']=='PASS_BANDED_COERCIVITY_REAL_COST':
            weight=F(c['reciprocal_weight_upper_exact']);assert weight>=1/(2*F(c['kappa_exact']))>0
            value=F(c['A_upper_exact'])*weight;assert value==F(c['cost_upper_exact'])
        else:
            value=max(F(c['first_upper_exact']),F(c['replay_upper_exact']));assert value==F(c['retained_upper_exact'])
        assert value>=0;fsum+=value
    assert cursor==1 and fsum==F(finite['cost_upper_exact'])
    real=read(ROOT/'output/competitive_20261007/left/real/tail450_split_replay.json')
    assert real['status']=='PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY'and real['witness_sha256']==SHA
    assert real['mode']=='tail'and real['first_even_spin']==450 and real['last_even_spin']=='infinity'
    cut=F(real['compact_beta_cutoff_exact']);assert cut**2<=F(45,49)
    sums={'compact':F(0),'remainder':F(0)};cursor=F(0)
    for c in real['cells']:
        lo,hi=F(c['lo']),F(c['hi']);assert cursor==lo<hi<=1;cursor=hi;k=c['region']
        assert (hi<=cut)if k=='compact'else(lo>=cut)
        a,b,value=map(F,[c['first_upper_exact'],c['replay_upper_exact'],c['retained_upper_exact']]);assert a>=0 and b>=0 and value==max(a,b)
        sums[k]+=value
    assert cursor==1 and all(sums[k]==F(real['A_upper_exact'][k])for k in sums)
    components={'head_0_40':retained,'finite_42_448_allenergy':fsum,'tail450_compact':sums['compact']/(2*ck),'tail450_UV_and_sliver':sums['remainder']/(2*uk)}
    total=sum(components.values(),F(0));decimal=upward(total);assert F(decimal)>=total
    check=read(ROOT/'output/competitive_20261007/left/real/independent_field_check.json')
    assert check['status']=='PASS_INDEPENDENT_REAL_KERNEL_COMPARISON'and check['witness_sha256']==SHA
    for p in [Path(__file__),ROOT/'research/all_direction_certification_20261007/check_real220_fields.py',ROOT/'research/all_direction_certification_20261007/REVIEW_UNIFORM_SOURCE_LEFT.md',ROOT/'research/all_direction_certification_20261007/head/REVIEW_LOCAL_COERCIVITY.md']:bind(p)
    dual=support_report(180);primal=support_report(180,'primal')
    out=dict(status='PASS_COMPLETE_CONTINUUM_LEFT_BOUND',witness_sha256=SHA,objective_exact=['-1','0'],
        inequality='g0 >= -'+decimal,support_upper_exact=str(total),support_upper_outward_6dp=decimal,
        components_exact={k:str(v)for k,v in components.items()},components_float={k:float(v)for k,v in components.items()},
        tail_upper_exact=str(total-retained),tail_upper_float=float(total-retained),
        coverage='Every even spin ell>=0 and every energy mu>4, including finite-width threshold and source endpoint caps, infinite energy and infinite-spin remainders.',
        graphical_comparison=dict(dual=dual,primal=primal,relative_to_dual_estimate_percent=100*(float(F(decimal))/dual['support']-1)),
        comparison_warning='Published graph extraction, not exact author optimizer data; scientific assumptions must match before treating this as universal dominance.',
        proof_scope='Internal computer-assisted directed interval certificate under the manuscript analytic/dispersion/crossing/unitarity assumptions. Exact witness, source identities and remainders, rational partitions, precision replays. Not formal verification, external review, or an optimality claim.',
        sha256_manifest=manifest)
    p=ROOT/'output/competitive_20261007/left/complete_left_certificate.json';p.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items()if k not in ['sha256_manifest','components_exact','support_upper_exact','tail_upper_exact','graphical_comparison']},indent=2))
if __name__=='__main__':main()
