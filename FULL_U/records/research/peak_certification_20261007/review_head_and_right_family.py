"""Independent source/adaptation and exact-accounting checks for peer work."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/peak_certification_20261007';hashes={}
def sha(p):
    p=Path(p);h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def bind(rel,expected=None):
    p=Path(rel);p=p if p.is_absolute()else ROOT/p
    name=str(p.relative_to(ROOT));actual=hashes.get(name)
    if actual is None:actual=sha(p);hashes[name]=actual
    assert expected is None or actual==expected,(name,actual,expected)
    return p
def read(rel):
    p=bind(rel);d=json.loads(p.read_text())
    for key in ['sha256_manifest','dependency_sha256','dependencies_sha256']:
        for name,h in d.get(key,{}).items():bind(name,h)
    return d
def dy(x):return F(int(x[0]))*F(2)**int(x[1])

# The upper orchestration uses original guarded models; verify current prepared
# witness, duplicated sources and exact coefficient export independently.
prep=read(OUT/'upper83/head/preparation.json');wp=bind('output/weak_regions_20261007/upper/angle83_repaired.npz',prep['witness_sha256']);w=np.load(wp)
definition=read(OUT/'upper83/head/exact_witness_definition.json');assert definition['witness_sha256']==sha(wp)
assert list(map(F,definition['direction_exact']))==list(map(lambda x:F(float(x)),w['direction']))
assert len(definition['coefficients'])==len(w['names'])==len(w['multipliers'])==220
for r,n,v in zip(definition['coefficients'],w['names'],w['multipliers']):
    assert r['row']==str(n) and F(r['exact_binary64_rational'])==F(float(v)) and float(v).hex()==r['hex']
for c in prep['source_rebuild_checks']:
    bind(OUT/'upper83/head/moments128'/c['file'],c['sha256']);bind(OUT/'upper83/head/moments128_replay'/c['file'],c['sha256']);assert c['bit_for_bit_rebuild']
budget=read(OUT/'upper83/head/budgets.json');assert [r['spin']for r in budget['rows']]==list(range(0,41,2))
assert [r['model']for r in budget['rows']]==['centered','centered']+['fast']*19
assert sum((F(r['slack'])for r in budget['rows']),F(0))==F(budget['slack_sum_exact'])
for r in prep['preflight']:
    assert r['precision_bits']in[256,384] and all(dy(v)>=0 for v in r['gap_upper_dyadic'])
for name in ['prepare_head.py','run_head.py']:bind(ROOT/'research/peak_certification_20261007/upper83'/name)
assert prep['script_sha256']==sha(ROOT/'research/peak_certification_20261007/upper83/prepare_head.py')

# Exactly the documented driver constants change in finite-source recombination.
source=ROOT/'research/full_direction_20261007/combine_free_point_finite.py'
skip=('CASES=','    ctx.prec=','    here=','    path=ROOT/')
canonical=[x for x in source.read_text().splitlines()if not x.startswith(skip)]
for name in ['combine_right_family.py','replay_right_finite.py']:
    p=ROOT/'research/peak_certification_20261007'/name
    assert [x for x in p.read_text().splitlines()if not x.startswith(skip)]==canonical
    bind(p)
bind(source)
family=read(OUT/'right_family/complete_right_family.json');assert family['status']=='PASS_EIGHT_FULL_CONTINUUM_DIRECTIONAL_BOUNDS'
base=read('output/full_direction_20261007/certification/complete_fifteen_direction_certificate.json')
tail=read('output/certification_tail/fixed_tail_certificate.json');assert tail['replayed_components_this_run']
wp=bind('output/angles_resummed/angle_000p000.npz',base['witness_sha256']);w=np.load(wp)
assert tail['witness_sha256']==sha(wp) and F(tail['V_over_p_upper_rational'])==-F(31,50)
finite=[read(OUT/'right_family'/f)for f in ['finite_family.json','finite_family_replay.json']]
assert [r['precision_bits']for r in finite]==[280,340]
slopes=list(map(F,['-20','-24','-28','-40','-48','-96','-192','-256']))
head={ell:[read(OUT/'right_family'/f'head_p{bits}_spin{ell}.json')for bits in[280,340]]for ell in[0,2,4,6]}
results=[]
for rel,embedded,b in zip(family['certificates'],family['results'],slopes):
    d=read(rel);assert d==embedded and d['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
    assert list(map(F,d['objective_exact']))==[F(1),b]
    definition=read(d['witness']);bind(d['witness'],d['witness_sha256'])
    assert definition['format']=='exact_objective_shift_of_frozen_witness' and definition['objective_exact']==d['objective_exact']
    assert definition['names']==w['names'].tolist() and len(definition['names'])==150
    assert list(map(F,definition['multipliers_exact']))==[F(float(v))for v in w['multipliers']]
    assert F(definition['crossing_increment_exact'])==0 and definition['base_witness_sha256']==sha(wp)
    components={}
    for ell,pair in head.items():
        vals=[];parts=[]
        for bits,h in zip([280,340],pair):
            assert h['status']=='PASS'and h['spin']==ell and h['precision_bits']==bits and h['witness_sha256']==sha(wp)
            bind('research/full_direction_20261007/head_family_extended.py',h['source_sha256'])
            row=next(r for r in h['directions']if F(r['slope_exact'])==b);assert F(row['crossing_increment_exact'])==0
            cells=row['cells']+row['gaps'];partition=sorted((F(c['a']),F(c['b']))for c in cells)
            assert partition[0][0]==0 and partition[-1][1]==1 and all(a<z for a,z in partition)
            assert all(x[1]==y[0]for x,y in zip(partition,partition[1:]));parts.append(partition)
            val=sum((F(c['upper_exact'])for c in cells),F(0));assert val==F(row['head_upper_exact']) and all(F(c['upper_exact'])>=0 for c in cells)
            if ell:assert dy(row['threshold']['new_V_over_p_upper_dyadic'])<0
            vals.append(val)
        assert parts[0]==parts[1];components[f'head_spin{ell}']=max(vals)
    vals=[]
    for packet in finite:
        assert packet['status']=='PASS' and packet['witness_sha256']==sha(wp)and packet['spins']==list(range(8,42,2))
        row=next(r for r in packet['results']if F(r['b'])==b);assert row['status']=='PASS'and F(row['delta_X4'])==0
        assert sorted(map(int,row['support_upper_by_spin_dyadic']))==list(range(8,42,2))
        val=F(row['finite_support_upper_exact']);assert val>=sum(map(dy,row['support_upper_by_spin_dyadic'].values()),F(0));vals.append(val)
    components['finite_spins8_40']=max(vals);k=F(31,50)-b;assert k==F(d['tail_kappa_exact'])and k>0
    components['tail_spins42_infinity']=F(tail['A40_upper_exact'])/(2*k)
    assert components=={k:F(v)for k,v in d['components_exact'].items()}
    total=sum(components.values(),F(0));assert total==F(d['support_upper_exact'])
    outward=F(d['support_upper_outward_12dp']);normalized=F(d['support_normalized_upper_outward_12dp'])
    norm2=1+b*b;assert total<=outward<total+F(1,10**12)and normalized**2*norm2>=outward**2
    scale=2**256;normlo=F(math.isqrt(norm2.numerator*scale*scale//norm2.denominator),scale)
    assert outward/normlo<=normalized<outward/normlo+F(1,10**12)
    results.append(dict(certificate=rel,sha256=sha(ROOT/rel),objective_exact=d['objective_exact'],normalized_upper=str(normalized)))
assert len(results)==8
bind(Path(__file__));bind(ROOT/'research/peak_certification_20261007/assemble_right_family.py')
report=dict(status='PASS_PEER_SOURCE_AND_EXACT_ACCOUNTING_REVIEW',upper_head_preparation='PASS_WITNESS_SOURCES_MODELS_AND_EXACT_COEFFICIENT_EXPORT; running head integrals not predeclared complete',right_family=results,checked_file_count=len(hashes),sha256_manifest=hashes,scope='Independent source review and exact rational accounting, paired partitions, witness identities and file bindings. No independent rederivation of physical kernels or rerun of every interval integration. Frozen LowSpin normalized shortcut is restricted below beta1/2; no new Centered source-domain shortcut is introduced.')
(OUT/'peer_head_and_right_family_review.json').write_text(json.dumps(report,indent=2));print(report['status'],len(results),'right-family certificates',len(hashes),'file hashes')
