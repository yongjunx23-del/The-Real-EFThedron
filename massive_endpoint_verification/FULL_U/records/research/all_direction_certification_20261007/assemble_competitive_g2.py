"""Assemble only complete, matching-witness all-spin/all-energy components."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research'))
from angle_benchmarks import support_report
SHA='933ee2f7889a52b8ccabbf58836c256161bc5c59d3763cdbca2b40bb01ac6f92'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def upward(x,n=6):
    v=-((-x.numerator*10**n)//x.denominator)
    return f'{v//10**n}.{v%10**n:0{n}d}'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--strong-uv',action='store_true');ap.add_argument('--refined-head',action='store_true')
    ap.add_argument('--head',type=Path,default=ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1/final_head_audit.json')
    ap.add_argument('--out',type=Path,default=ROOT/'output/competitive_20261007/complete_g2_certificate.json')
    args=ap.parse_args();manifest={}
    def bind(path):
        path=Path(path);key=str(path.relative_to(ROOT));h=sha(path)
        if key in manifest:assert manifest[key]==h,key
        manifest[key]=h
    def read(path):
        bind(path);d=json.loads(Path(path).read_text());assert d['witness_sha256']==SHA
        for p,h in d.get('sha256_manifest',{}).items():
            path=Path(p);path=path if path.is_absolute() else ROOT/path
            assert sha(path)==h,p;bind(path)
        return d
    head=read(args.head.resolve());assert head['status']=='COMPLETE_REPLAYED_FINITE_HEAD'
    assert head['spins_included']==list(range(0,41,2)) and not head['missing']
    for p,h in head['all_head_code_sha256'].items():assert sha(ROOT/p)==h;bind(ROOT/p)
    hsum=F(0)
    for row in head['records']:
        p=ROOT/row['record'];q=ROOT/row['replay'];assert sha(p)==row['record_sha256'] and sha(q)==row['replay_sha256'];bind(p);bind(q)
        record=json.loads(p.read_text());replay=json.loads(q.read_text())
        assert record['witness_sha256']==replay['witness_sha256']==SHA and replay['passed']
        assert replay['partition_covers_beta_zero_to_one_exactly'] and replay['source_record_sha256']==sha(p)
        for sub in [record.get('dependencies',{}),replay.get('dependency_sha256',{})]:
            for name,digest in sub.items():assert sha(ROOT/name)==digest;bind(ROOT/name)
        m,e=map(int,replay['upper_dyadic']);value=F(m)*F(2)**e
        assert value==F(row['replayed_upper_exact']);hsum+=value
    assert hsum==F(head['upper_exact'])
    if args.refined_head:
        base=args.head.resolve().parent
        selected=read(base/'head_refined_summary.json')
        assert selected['status']=='COMPLETE_REPLAYED_REFINED_FINITE_HEAD'
        assert selected['original_replay_audit_sha256']==sha(args.head.resolve())
        record_path=base/'low0_refined/head_l0.json';replay_path=base/'low0_refined/replay_l0.json'
        record=read(record_path);replay=read(replay_path)
        assert replay['passed'] and replay['partition_covers_beta_zero_to_one_exactly']
        assert replay['source_record_sha256']==sha(record_path)==selected['records'][0]['sha256']
        assert replay['source_record_sha256']==sha(ROOT/selected['records'][0]['path'])
        assert selected['refined_spin0_replay_sha256']==sha(replay_path)
        for sub in [record['dependencies'],replay['dependency_sha256']]:
            for name,digest in sub.items():assert sha(ROOT/name)==digest;bind(ROOT/name)
        m,e=map(int,replay['upper_dyadic']);value=F(m)*F(2)**e
        old0=next(row for row in head['records'] if row['spin']==0)
        hsum+=value-F(old0['replayed_upper_exact'])
        assert hsum==F(selected['replay_upper_exact'])
    compact=read(ROOT/'output/competitive_20261007/compact/complete_compact_coercivity.json')
    assert compact['status']=='PASS_COMPLETE_COMPACT_COERCIVITY' and F(compact['kappa_exact'])==1
    assert compact['objective_exact']==['0','1'] and compact['mu_domain_exact']==['4','49']
    assert compact['minimum_even_spin']==42 and compact['lower_energy_endpoint_open']
    uvbase=ROOT/'output/all_direction_certification_20261007/tail/g2_robust_stage1'
    uv=read(uvbase/('strong_uv_coercivity.json' if args.strong_uv else 'complete_uv_coercivity.json'))
    assert uv['objective_exact']==['0','1'] and uv['energy_domain']==['49','infinity']
    if args.strong_uv:
        assert uv['status']=='PASS_STRONG_FINITE_UV_COERCIVITY' and uv['finite_spin_domain']==[42,448]
        eta_finite=F(uv['finite_kappa_exact']);eta_tail=F(uv['infinite_kappa_exact'])
    else:
        assert uv['status']=='PASS_COMPLETE_UV_COERCIVITY'
        eta_finite=eta_tail=F(uv['kappa_exact'])
    assert 0<eta_finite<=1 and 0<eta_tail<=1
    components={'head_0_40':hsum};norm_records=[]
    for mode,name,first,last,eta in [('finite','finite42_448_split_replay.json',42,448,eta_finite),('tail','tail450_split_replay.json',450,'infinity',eta_tail)]:
        q=read(ROOT/'output/competitive_20261007/real'/name)
        assert q['status']=='PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY' and q['mode']==mode
        assert q['first_even_spin']==first and q['last_even_spin']==last
        cutoff=F(q['compact_beta_cutoff_exact']);assert cutoff**2<=F(45,49)
        totals={'compact':F(0),'remainder':F(0)};cursor=F(0)
        for c in q['cells']:
            lo,hi=F(c['lo']),F(c['hi']);assert cursor==lo<hi<=1
            region=c['region'];assert region in totals
            assert (hi<=cutoff) if region=='compact' else (lo>=cutoff)
            a,b,retained=map(F,[c['first_upper_exact'],c['replay_upper_exact'],c['retained_upper_exact']])
            assert a>=0 and b>=0 and retained==max(a,b)
            totals[region]+=retained;cursor=hi
        assert cursor==1 and all(totals[k]==F(q['A_upper_exact'][k]) for k in totals)
        components[mode+'_compact']=totals['compact']/2
        components[mode+'_UV_and_cutoff_sliver']=totals['remainder']/(2*eta)
        norm_records.append(dict(mode=mode,eta_compact_exact='1',eta_UV_exact=str(eta),A_exact={k:str(v) for k,v in totals.items()}))
    review=read(ROOT/'output/competitive_20261007/real/highspin_v_audit.json')
    assert review['status']=='PASS_HIGH_SPIN_PACKET_AND_INDEPENDENT_JET_COMPARISON'
    for p in [Path(__file__),ROOT/'research/all_direction_certification_20261007/REVIEW_HIGHSPIN220.md',ROOT/'research/all_direction_certification_20261007/REVIEW_FINITE_UV220.md',ROOT/'research/all_direction_certification_20261007/head/REVIEW_REAL220.md']:
        bind(p)
    total=sum(components.values(),F(0));decimal=upward(total);rounded=F(decimal);assert rounded>=total
    dual=support_report(90);primal=support_report(90,'primal')
    result=dict(status='PASS_COMPLETE_CONTINUUM_G2_BOUND',witness_sha256=SHA,objective_exact=['0','1'],
        inequality=f'g2 <= {decimal}',support_upper_exact=str(total),support_upper_outward_6dp=decimal,
        components_exact={k:str(v) for k,v in components.items()},components_float={k:float(v) for k,v in components.items()},
        real_norm_components=norm_records,finite_UV_kappa_exact=str(eta_finite),infinite_UV_kappa_exact=str(eta_tail),
        refined_head_used=args.refined_head,
        coverage='Every even spin ell>=0 and every mu>4, with explicit threshold, source endpoint, infinite-energy and infinite-spin remainders.',
        proof='Disk support H(U,V)=sqrt(U²+V²)+V. Finite0..40 enclosed directly. For higher spins, V<=-eta*p implies H<=U²/(2*eta*p). The real norm partitions and coercivity proofs cover all energies and spins. The tiny rational-cutoff sliver belowmu49 uses the weaker UV eta, which is implied by the compact eta1.',
        graphical_comparison=dict(dual=dual,primal=primal,certified_minus_dual_estimate=float(rounded)-dual['support'],relative_to_dual_estimate_percent=100*(float(rounded)/dual['support']-1)),
        comparison_warning='The published curve is extracted graphically from the supplied PDF; it is not a supplied exact dual witness or independently certified numerical interval.',
        verification_scope='Exact-rational witness interpretation, directed Arb interval arithmetic, analytic continuum tails, higher-precision replays and internal independent mathematical/kernel reviews. Conditional on the dispersion, crossing, analyticity and unitarity assumptions used to derive the null equations. Not formal theorem-prover verification or an external audit.',
        sha256_manifest=manifest)
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['sha256_manifest','components_exact','support_upper_exact','real_norm_components','graphical_comparison']},indent=2))


if __name__=='__main__':main()
