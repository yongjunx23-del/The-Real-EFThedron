"""Assemble the fixed-witness all-spin upper bound with exact rational totals.

Default: independently reevaluate every finite-head cell at 280-bit precision.
--assemble-only checks stored enclosures and fresh endpoint caps instead.
The separate all-spin tail verifier and finite-V verifier supply the analytic
tail components; their hashes and exact domain partitions are checked here.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,argparse,time
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from flint import arb,ctx
from head_certification_spin0 import SpinZero
from head_certification_lowspin import LowSpin

EXPECTED='bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c'
WIT=ROOT/'output/angles_resummed/angle_000p000.npz'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(pair):return F(int(pair[0]))*F(2)**int(pair[1])
def upper(x):
    assert x.is_finite()
    return dyad(x.upper().man_exp())
def read(p):return json.loads(Path(p).read_text())
def relative_path(p):return ROOT/str(p).replace('\\','/')
def partition(parts,left=F(0),right=F(1)):
    parts=sorted(parts);assert parts and parts[0][0]==left and parts[-1][1]==right
    assert all(a<b for a,b in parts)
    assert all(b==c for (a,b),(c,d) in zip(parts,parts[1:]))
def outward_decimal(value,digits=6):
    factor=10**digits;v=value*factor;n=-(-v.numerator//v.denominator)
    return f'{n//factor}.{n%factor:0{digits}d}'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--assemble-only',action='store_true')
    ap.add_argument('--v-certificate',type=Path)
    ap.add_argument('--v-replay-record',type=Path,help='Independent exact all-spin leaf replay record, bound to the selected certificate')
    ap.add_argument('--out',type=Path,default=ROOT/'output/certification_head/complete_upper_certificate.json')
    args=ap.parse_args();ctx.prec=280;tic=time.time();assert sha(WIT)==EXPECTED
    inputs={str(WIT.relative_to(ROOT)):sha(WIT)}
    def use(path):
        path=Path(path).resolve();inputs[str(path.relative_to(ROOT))]=sha(path);return read(path)
    mc=use(ROOT/'output/certification_head/all_moment_tables_verified.json')
    assert mc['status']=='PASS' and mc['witness_sha256']==EXPECTED
    verified_tables={}
    for r in mc['tables']:
        p=relative_path(r['path']);assert sha(p)==r['sha256'];inputs[r['path']]=r['sha256']
        verified_tables[p.resolve()]=r['sha256']
    tail=use(ROOT/'output/certification_tail/fixed_tail_certificate.json')
    assert tail['witness_sha256']==EXPECTED and tail['replayed_components_this_run']
    for p,h in tail['sha256_manifest'].items():
        assert sha(relative_path(p))==h,('tail hash',p);inputs[p]=h
    tail_A=F(tail['A40_upper_exact']);tail_kappa=-F(tail['V_over_p_upper_rational'])
    assert tail_A>0 and tail_kappa>0
    tail_cost=tail_A/(2*tail_kappa)
    assert tail_cost==F(tail['disk_support_tail_upper_exact']) and tail_cost<F(1,62000)
    vpath=args.v_certificate
    if vpath is None:
        vpath=ROOT/'output/certification_head/V8_40_kappa14_exact.json'
        if not vpath.exists():vpath=ROOT/'output/certification_head/V8_40_exact.json'
    vpath=Path(vpath).resolve()
    vv=use(vpath);assert vv['witness_sha256']==EXPECTED and vv['spins']==list(range(8,42,2))
    kappa=F(vv['kappa']);assert kappa>0
    replaypath=args.v_replay_record or vpath.with_name(vpath.stem.replace('_exact','')+'_independent_replay.json')
    vr=use(replaypath)
    assert vr['status']=='PASS' and vr['certificate_sha256']==sha(vpath)
    assert vr['witness_sha256']==EXPECTED and F(vr['kappa'])==kappa
    assert vr.get('cells',vr['leaves'])==len(vv['leaves'])
    assert vr['spin_cell_inequalities']==17*len(vv['leaves']) and vr['exact_gapless_domain_partitions']
    if 'every_spin_finite_and_below_threshold' in vr:assert vr['every_spin_finite_and_below_threshold'] is True
    replay_hashes={p.replace('\\','/'):h for p,h in vr['sha256_manifest'].items()}
    assert replay_hashes['research/tail_certification/certified_finite_v.py']==sha(ROOT/'research/tail_certification/certified_finite_v.py')
    assert replay_hashes['research/tail_certification/certify_finite_v_domain.py']==sha(ROOT/'research/tail_certification/certify_finite_v_domain.py')
    for p,h in vr['sha256_manifest'].items():
        pp=relative_path(p).resolve();assert sha(pp)==h;inputs[str(pp.relative_to(ROOT))]=h
    for name in ('roy_moments.json','strip_12_24_moments.json','strip_24_48_moments.json'):
        pp=(ROOT/'output/certification_head'/name).resolve();assert pp in verified_tables
    for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
        cells=[c for c in vv['leaves'] if c['variable']==kind]
        partition([(F(c['a']),F(c['b'])) for c in cells],lo,hi)
        assert all(dyad(c['upper_dyadic'])<-kappa for c in cells)
    assert sum(c['variable'] in ('mu','w') for c in vv['leaves'])==len(vv['leaves'])
    aa=use(ROOT/'output/certification_tail/finite_A8_40_tighter.json')
    assert aa['witness_sha256']==EXPECTED and aa['spins']==vv['spins']
    partition([(F(c['a']),F(c['b'])) for c in aa['partition']])
    A=dyad(aa['A_upper_dyadic']);assert sum((dyad(c['mass_upper']) for c in aa['partition']),F(0))<=A
    assert A<F(235,1000)
    block=A/(2*kappa)
    block_method='uniform coercivity'
    local_slack=F(0)
    localpath=ROOT/'output/certification_head/block8_40_local_coercivity.json'
    if localpath.exists() and kappa==14:
        local=use(localpath)
        assert local['status']=='PASS' and local['witness_sha256']==EXPECTED and local['spins']==vv['spins']
        for p,h in local['sha256'].items():assert sha(relative_path(p))==h;inputs[p]=h
        assert local['sha256']['output/certification_head/V8_40_kappa14_exact.json']==sha(vpath)
        lr=use(ROOT/'output/certification_head/V8_40_kappa14_local_replay.json')
        assert lr['status']=='PASS' and lr['certificate_sha256']==sha(vpath) and lr['witness_sha256']==EXPECTED
        assert lr['spins']==vv['spins'] and lr['cells']==len(vv['leaves']) and lr['spin_cell_inequalities']==17*len(vv['leaves'])
        assert lr['every_spin_finite'] and lr['every_spin_upper_le_saved_leaf_upper_plus_slack'] and lr['exact_gapless_domain_partitions']
        for p,h in lr['sha256_manifest'].items():assert sha(relative_path(p))==h;inputs[p]=h
        local_slack=F(lr['uniform_slack_exact']);assert local_slack==F(1,1000000)
        segments=[]
        for c in vv['leaves']:
            a,b=F(c['a']),F(c['b']);kap=-dyad(c['upper_dyadic'])
            lo,hi=(a,b) if c['variable']=='mu' else (4+1/b,4+1/a if a else None)
            segments.append((lo,hi,kap))
        segments.sort();assert segments[0][0]==4 and segments[-1][1] is None
        assert all(a[1]==b[0] for a,b in zip(segments,segments[1:]))
        acells=sorted(aa['partition'],key=lambda c:F(c['a']));lcells=local['partition']
        assert len(acells)==len(lcells);j=0;subtotal=F(0);relaxed_total=F(0)
        for ac,lc in zip(acells,lcells):
            a,b=F(ac['a']),F(ac['b']);assert (a,b)==(F(lc['a']),F(lc['b']))
            lo,hi=4/(1-a*a),4/(1-b*b) if b<1 else None
            while segments[j][1] is not None and segments[j][1]<=lo:j+=1
            stop=j
            while stop+1<len(segments) and (hi is None or segments[stop+1][0]<hi):stop+=1
            kap=min(s[2] for s in segments[j:stop+1])
            assert F(lc['coercivity_exact'])==kap
            mass=dyad(ac['mass_upper']);cost=dyad(lc['support_upper_dyadic'])
            assert mass>=0 and cost>=mass/(2*kap)
            subtotal+=cost
            relaxed_kappa=kap-local_slack;assert relaxed_kappa>0
            rational_cost=mass/(2*relaxed_kappa)
            relaxed_total+=upper(arb(rational_cost.numerator)/rational_cost.denominator)
        localupper=dyad(local['support_upper_dyadic']);assert subtotal<=localupper<block
        assert localupper<=relaxed_total<block
        block=relaxed_total;block_method='local coercivity with independently replayed per-leaf bounds, explicit slack, exact overlap and directed division'
    head=[]
    for ell in (0,2,4,6):
        name='spin0_component_tight.json' if ell==0 else f'spin{ell}_tight_component.json'
        p=ROOT/'output/certification_head'/name;r=use(p)
        assert r['witness_sha256']==EXPECTED and r['target_met']
        code=ROOT/'research'/('head_certification_spin0.py' if ell==0 else 'head_certification_lowspin.py')
        assert r['source_sha256']==sha(code),('head source hash',ell)
        if ell==0:assert r['real_tail_source_sha256']==sha(ROOT/'research/certify_real_tail.py')
        inputs[str(code.relative_to(ROOT))]=sha(code)
        for dep,h in r.get('dependency_sha256',{}).items():
            assert sha(ROOT/'research'/dep)==h;inputs['research/'+dep]=h
        for t in r.get('source_tables',[]):
            for kind in ('moment','uv'):
                # Metadata records the original workstation path; bind its
                # output-relative suffix for a relocated standalone bundle.
                suffix=t[kind+'_path'].replace('\\','/').split('/output/',1)[1]
                pp=(ROOT/'output'/suffix).resolve()
                assert sha(pp)==t[kind+'_sha256'] and pp in verified_tables
        gaps=[(F(a),F(b)) for a,b in r['gap_beta_intervals']]
        cells=[(F(c['a']),F(c['b'])) for c in r['cells_data']]
        partition(cells+gaps)
        model=SpinZero(WIT,r['N']) if ell==0 else LowSpin(WIT,ell)
        gap=sum((upper(model.gap_bound(a,b)) for a,b in gaps),F(0))
        if args.assemble_only:
            bulk=sum((dyad(c['upper_dyadic']) for c in r['cells_data']),F(0))
        else:
            bulk=F(0)
            for i,(a,b) in enumerate(cells):
                value,_=model.trapezoid(a,b);mass=upper(value)
                assert mass>=0
                bulk+=mass
                if (i+1)%1000==0:print(json.dumps({'spin':ell,'replayed_cells':i+1}),flush=True)
        cost=bulk+gap;assert cost<F(r['target']),(ell,float(cost),r['target'])
        if ell:
            assert dyad(r['threshold']['V_over_p_upper_dyadic'])<0
        head.append({'spin':ell,'upper_exact':str(cost),'upper_outward_12dp':outward_decimal(cost,12),
                     'cells':len(cells),'endpoint_upper_exact':str(gap)})
        print(json.dumps({'spin':ell,'certified_upper':outward_decimal(cost,12)}),flush=True)
    total=sum((F(h['upper_exact']) for h in head),F(0))+block+tail_cost
    assert total<F('2.72322')
    for name in ('verify_complete_upper_bound.py','verify_moment_tables.py','certified_roy_moments.py','certified_strip_moments.py','certified_uv_moments.py','validated_quadrature.py'):
        p=ROOT/'research'/name;inputs[str(p.relative_to(ROOT))]=sha(p)
    for name in ('certified_finite_v.py','certify_finite_v_domain.py','certify_finite_real.py','FINITE_HEAD_COMPONENTS.txt'):
        p=ROOT/'research/tail_certification'/name;inputs[str(p.relative_to(ROOT))]=sha(p)
    out={'status':'PASS','claim':'g0 is at most the outward-rounded bound under the manuscript amplitude assumptions.',
         'witness_sha256':EXPECTED,'head_cells_replayed_this_run':not args.assemble_only,'precision_bits':ctx.prec,
         'finite_V_full_independent_replay_record':str(Path(replaypath).resolve().relative_to(ROOT)),
         'head_components':head,'block8_40_A_upper_exact':str(A),'block8_40_kappa':str(kappa),
         'block8_40_support_upper_exact':str(block),'tail42_infinity_upper_exact':str(tail_cost),
         'block8_40_method':block_method,
         'block8_40_local_V_upper_slack_exact':str(local_slack),
         'complete_support_upper_exact':str(total),'certified_g0_upper_outward_6dp':outward_decimal(total),
         'certified_g0_upper_outward_12dp':outward_decimal(total,12),'manuscript_reference':'2.72322','published_reference':'2.7272',
         'covers':'all even physical spins and every mu>4, including threshold and infinite energy',
         'new_266_row_witness_certified':False,'seconds':time.time()-tic,'sha256_manifest':inputs}
    args.out.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='sha256_manifest'},indent=2))

if __name__=='__main__':main()
