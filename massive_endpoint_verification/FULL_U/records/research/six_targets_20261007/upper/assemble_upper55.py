"""Exact positive-mixture full certificate with only spins0 and2 reintegrated."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,math,sys,re
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'research'),str(ROOT/'tmp/certification_deps')]
from flint import arb,ctx
from angle_benchmarks import support_report
from upper55_common import align_gaps
OUT=ROOT/'output/six_targets_20261007/upper55'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dy(p):return F(int(p[0]))*F(2)**int(p[1])
def decimal_ball_upper(text):
    m=re.fullmatch(r'\[(.+) \+/- (.+)\]',text);assert m is not None
    return F(m[1])+F(m[2])
def upward(v,n=9):
    q=-((-v.numerator*10**n)//v.denominator);return f'{q//10**n}.{q%10**n:0{n}d}'
def partition(cells,gaps):
    z=sorted(cells+gaps);assert z[0][0]==0 and z[-1][1]==1 and all(a<b for a,b in z)and all(a[1]==b[0]for a,b in zip(z,z[1:]))
class Inputs:
    def __init__(self):self.manifest={};self.seen=set()
    def bind(self,p,h=None):
        p=Path(p);p=(p if p.is_absolute()else ROOT/p).resolve();key=str(p.relative_to(ROOT));actual=sha(p)
        assert h is None or actual==h,key
        assert key not in self.manifest or self.manifest[key]==actual,key
        self.manifest[key]=actual;return p
    def read(self,p,h=None):
        p=self.bind(p,h);d=json.loads(p.read_text())
        if p not in self.seen:
            self.seen.add(p)
            if isinstance(d,dict):
                for q,digest in d.get('sha256_manifest',{}).items():
                    qq=self.bind(q,digest)
                    if qq.suffix=='.json':self.read(qq,digest)
        return d
    def dependencies(self,d):
        for label in ['dependencies','dependency_sha256','assembly_and_source_code_sha256','all_head_code_sha256']:
            for p,h in d.get(label,{}).items():self.bind(p,h)

def headgeneric(i,c,path):
    d=i.read(path);assert d['status']=='COMPLETE_REPLAYED_FINITE_HEAD'and d['witness_sha256']==c['witness_sha256'];i.dependencies(d)
    assert d['included_spins']==list(range(0,41,2));assert not d['missing_spins']and not d['unfinished_spins']and not d['missing_replays']
    costs={};records={}
    for r in d['records']:
        p=i.read(r['path'],r['sha256']);q=i.read(r['replay'],r['replay_sha256']);assert p['witness_sha256']==q['witness_sha256']==c['witness_sha256'];assert p['spin']==q['spin']==r['spin']
        assert p['status']=='PASS_FINITE_SPIN_ENCLOSURE'and p['final'];assert q['status']=='PASS_FULL_PARTITION_REPLAY'and q['passed']and q['partition_covers_beta_zero_to_one_exactly'];assert q['source_record_sha256']==sha(ROOT/r['path']);assert q['precision_bits']>p['precision_bits']
        i.dependencies(p);i.dependencies(q);costs[r['spin']]=max(dy(p['upper_dyadic']),dy(q['upper_dyadic']));records[r['spin']]=p
    assert sorted(costs)==list(range(0,41,2));assert sum(costs.values(),F(0))==F(c['components_exact']['head_0_40'])==F(d['conservative_max_of_original_replay_upper_exact']);return costs,records


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert not a.out.exists();ctx.prec=512;i=Inputs()
    paths=['output/peak_certification_20261007/complete_upper83_certificate.json','output/boundary_length_20261007/complete_angle85_banded_certificate.json'];c=[i.read(p)for p in paths]
    assert all(cc['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'for cc in c)
    for cc in c:assert sum(map(F,cc['components_exact'].values()),F(0))==F(cc['support_upper_exact'])
    wp=['output/weak_regions_20261007/upper/angle83_repaired.npz','output/full_direction_20261007/right/angle85_220.npz'];w=[np.load(i.bind(p,cc['witness_sha256']))for p,cc in zip(wp,c)];assert np.array_equal(w[0]['names'],w[1]['names'])and len(w[0]['names'])==220
    direction_file='output/six_targets_20261007/upper55/objective.npz';target=np.load(i.bind(direction_file));direction=list(map(lambda q:F(float(q)),target['direction']));bd=[list(map(F,cc['objective_exact']))for cc in c]
    for ww,dd in zip(w,bd):assert list(map(lambda q:F(float(q)),ww['direction']))==dd
    det=bd[0][0]*bd[1][1]-bd[0][1]*bd[1][0]
    weights=[(direction[0]*bd[1][1]-direction[1]*bd[1][0])/det,(bd[0][0]*direction[1]-bd[0][1]*direction[0])/det]
    assert min(weights)>0 and[sum(weights[j]*bd[j][k]for j in[0,1])for k in[0,1]]==direction
    multipliers=[sum(weights[j]*F(float(w[j]['multipliers'][k]))for j in[0,1])for k in range(220)]
    hc=[];source=[]
    for cc,path in zip(c,['output/peak_certification_20261007/upper83/head/head_summary_selected.json','output/boundary_length_20261007/head/angle85/head_summary.json']):
        costs,recs=headgeneric(i,cc,path);hc.append(costs);source.append(recs)
    moments=['output/peak_certification_20261007/upper83/head/moments128','output/boundary_length_20261007/head/angle85/moments128']
    for n,m in enumerate(moments):
        mp=i.bind(Path(m)/'manifest.json');md=i.read(mp);assert md['witness_sha256']==c[n]['witness_sha256']
        for ell in[0,2]:assert source[n][ell]['source_manifest_sha256']==sha(mp)
        for table in md['tables']:
            label=table['label'];p=i.bind(Path(m)/(label+'_moments.json'),table['moment_sha256']);q=i.bind(Path(m)/(label+'_uv_moments.json'),table['uv_sha256']);mm=json.loads(p.read_text())['metadata'];um=json.loads(q.read_text())['metadata'];assert mm['witness_sha256']==um['witness_sha256']==c[n]['witness_sha256'];assert 2*mm['order']-1>=max(mm['energy_degree_bound'],mm['angular_degree_bound'])
    updates=[];retained=[]
    for ell in[0,2]:
        aligned,coarsening=align_gaps([source[j][ell]for j in[0,1]],ell)
        pair=[i.read(OUT/f'head_l{ell}_{bits}.json')for bits in[256,384]]
        geometry=None
        for r,bits in zip(pair,[256,384]):
            assert r['status']=='PASS_EXACT_POSITIVE_MIXTURE_HEAD_COMPONENT'and r['spin']==ell and r['precision_bits']==bits;i.dependencies(r)
            assert all(cell['method']!='mixture_negative_V'or F(cell['b'])<F(1,2)or(ell and F(cell['b'])**2<=F(79,119))for cell in r['cells'])
            assert list(map(F,r['objective_exact']))==direction and list(map(F,r['weights_exact']))==weights and list(map(F,r['multipliers_exact']))==multipliers
            assert r['names']==w[0]['names'].tolist()and r['base_witness_sha256']==[x['witness_sha256']for x in c];assert [str(Path(x))for x in r['base_certificate_paths']]==[str(Path(x))for x in paths]
            expected_sources=[ROOT/f'output/peak_certification_20261007/upper83/head/{"spin2_refined"if ell==2 else "spin0"}/head_l{ell}.json',ROOT/f'output/boundary_length_20261007/head/angle85/spin{ell}/head_l{ell}.json']
            assert [(ROOT/Path(p)).resolve()for p in r['source_record_paths']]==[p.resolve()for p in expected_sources]
            assert r['coarsening']==coarsening
            gaps=[tuple(map(F,g))for g in r['gap_intervals']];assert gaps==[tuple(map(F,g))for g in aligned[0]['gap_beta_intervals']]==[tuple(map(F,g))for g in aligned[1]['gap_beta_intervals']]
            gap_components=list(map(F,r['gap_components_exact']))
            for j in[0,1]:
                s=aligned[j];oldgap=dy(s['gap_upper_dyadic'])if 'gap_upper_dyadic'in s else decimal_ball_upper(s['gap_upper']);assert gap_components[j]>=oldgap
            assert F(r['gap_upper_exact'])==sum(weights[j]*gap_components[j]for j in[0,1]);cells=[(F(q['a']),F(q['b']))for q in r['cells']];partition(cells,gaps)
            if geometry is None:geometry=cells
            else:assert geometry==cells
            assert all(F(q['upper_exact'])>=0 and q['method']in['convex_base_cells','mixture_riemann','mixture_convex_endpoints','mixture_negative_V']for q in r['cells'])
            assert F(r['upper_exact'])==F(r['gap_upper_exact'])+sum((F(q['upper_exact'])for q in r['cells']),F(0))
        raw=max(F(r['upper_exact'])for r in pair);convex=sum(weights[j]*hc[j][ell]for j in[0,1]);keep=min(raw,convex);retained.append(keep);updates.append(dict(spin=ell,paired_direct_upper_exact=str(raw),base_retained_components_exact=[str(h[ell])for h in hc],old_convex_upper_exact=str(convex),retained_upper_exact=str(keep),cells=len(geometry),gaps=len(pair[0]['gap_intervals'])))
    old_total=sum(weights[j]*F(c[j]['support_upper_exact'])for j in[0,1]);subtract=sum(weights[j]*sum((hc[j][ell]for ell in[0,2]),F(0))for j in[0,1]);remainder=old_total-subtract;assert remainder>=0;total=remainder+sum(retained,F(0));assert 0<=total<=old_total
    norm2=sum(x*x for x in direction);scale=2**256;normlo=F(math.isqrt(norm2.numerator*scale*scale//norm2.denominator),scale);assert 0<normlo and normlo**2<=norm2<(normlo+F(1,scale))**2
    outward=upward(total);normal=upward(F(outward)/normlo);angle=math.degrees(math.atan2(float(direction[1]),float(direction[0])))%360
    definition={'witness_format':'exact_rational_positive_combination','base_certificate_paths':paths,'base_witness_paths':wp,'base_witness_sha256':[cc['witness_sha256']for cc in c],'weights_exact':list(map(str,weights)),'objective_exact':list(map(str,direction)),'names':w[0]['names'].tolist(),'multipliers_exact':list(map(str,multipliers))};dp=OUT/'mixture_definition.json';assert not dp.exists();dp.write_text(json.dumps(definition,indent=2));i.bind(dp)
    i.read(OUT/'preparation.json');i.read(OUT/'assembly_adaptation.json');runs=i.read(OUT/'run_manifest.json');assert runs['status']=='COMPLETE_PAIRED_LOW_SPIN_MIXTURE'and len(runs['jobs'])==4 and all(j['exit_code']==0 for j in runs['jobs'])
    for job in runs['jobs']:i.read(OUT/f"command_l{job['spin']}_{job['bits']}.json")
    i.bind(ROOT/'research/six_targets_20261007/upper/run_upper55.py',runs['script_sha256'])
    i.bind(ROOT/'research/six_targets_20261007/upper/upper55_common.py')
    i.bind(Path(__file__));i.bind(ROOT/'research/weak_regions_20261007/REVIEW_MIXTURE_HEAD.md');i.bind(ROOT/'research/six_targets_20261007/REVIEW_UPPER55_MIXTURE.md');result={'status':'PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND','witness':str(dp.relative_to(ROOT)),'witness_format':definition['witness_format'],'witness_sha256':sha(dp),'objective_exact':list(map(str,direction)),'nominal_angle_degrees':angle,'weights_exact':list(map(str,weights)),'support_upper_exact':str(total),'support_upper_outward_9dp':outward,'support_normalized_upper_outward_9dp':normal,'support_normalized_float':float(normal),'direction_norm_squared_exact':str(norm2),'direction_norm_lower_exact':str(normlo),'inequality':f'({direction[0]}) g0+({direction[1]}) g2 <= {outward}','old_convex_bound_exact':str(old_total),'removed_components_exact':str(subtract),'retained_other_spins_exact':str(remainder),'refined_components':updates,'components_exact':{'all_other_spins_and_tails':str(remainder),'spin0':str(retained[0]),'spin2':str(retained[1])},'coverage':'Every even spin and every energy above threshold. Spins0,2 integrated as an exact mixture; all other head and tail pieces bounded by positive homogeneity and subadditivity of disk support.','proof_scope':'Exact rational mixture of already completed witnesses, two directed-interval precisions for reintegrated spins, exact verification of base additive decomposition before subtraction. Internal shared-code review; not external or formal verification.','graphical_comparison':{'dual':support_report(angle),'primal':support_report(angle,'primal')},'sha256_manifest':i.manifest};a.out.write_text(json.dumps(result,indent=2));print(result['status'],normal,'oldconvex',float(old_total),'exactcomponents', {k:float(F(v))for k,v in result['components_exact'].items()})
if __name__=='__main__':main()

