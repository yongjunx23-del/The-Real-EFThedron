"""Tighten completed lower15/45 packets by their exactly bisected spin0 proof."""
from pathlib import Path
from fractions import Fraction as F
import json,sys,bisect,math
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/weak_regions_20261007'))
from assemble_mixture86 import Inputs,sha,dy,upward,partition
OUT=ROOT/'output/six_targets_20261007/lower'
def main(tag):
 output=OUT/f'complete_{tag}_refined_certificate.json';assert not output.exists();i=Inputs()
 oldpath=OUT/f'complete_{tag}_certificate.json';old=i.read(oldpath);assert old['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
 dp=ROOT/old['witness'];definition=i.read(dp,old['witness_sha256']);weights=list(map(F,definition['weights_exact']))
 oldhead=i.read(OUT/tag/'head_l0_256.json');assert oldhead['definition_sha256']==sha(dp)
 expected=[]
 for c in oldhead['cells']:
  lo,hi=F(c['a']),F(c['b']);expected.extend([(lo,(lo+hi)/2),((lo+hi)/2,hi)])
 gaps=[tuple(map(F,g))for g in oldhead['gap_intervals']];partition(expected,gaps)
 base_records=[i.read(p)for p in oldhead['source_record_paths']]
 base_cells=[[(F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in r['cells_data']]for r in base_records]
 starts=[[c[0]for c in row]for row in base_cells]
 pair=[i.read(OUT/'spin0_refined'/tag/f'head_l0_{bits}.json')for bits in[256,384]]
 for r,bits in zip(pair,[256,384]):
  i.dependencies(r);assert r['status']=='PASS_EXACT_POSITIVE_MIXTURE_HEAD_COMPONENT'and r['case']==tag and r['spin']==0 and r['precision_bits']==bits
  assert r['definition_sha256']==sha(dp)
  for key in['objective_exact','weights_exact','multipliers_exact','names','base_certificate_paths','base_witness_sha256','source_record_paths','gap_intervals','gap_components_exact','gap_upper_exact']:
   assert r[key]==oldhead[key],key
  assert [(F(c['a']),F(c['b']))for c in r['cells']]==expected
  for c in r['cells']:
   lo,hi=F(c['a']),F(c['b']);assert F(c['upper_exact'])>=0
   assert c['method']in['convex_base_cells','mixture_riemann','mixture_convex_endpoints','mixture_negative_V']
   if c['method']=='mixture_negative_V':assert hi<F(1,2)
   if c['method']=='convex_base_cells':
    fallback=F(0)
    for j in[0,1]:
     index=bisect.bisect_right(starts[j],lo)-1;a,b,v=base_cells[j][index];assert a<=lo<hi<=b;fallback+=weights[j]*v
    assert fallback==F(c['upper_exact'])
  assert F(r['upper_exact'])==F(r['gap_upper_exact'])+sum((F(c['upper_exact'])for c in r['cells']),F(0))
 campaign=i.read(OUT/'spin0_refined/campaign.json');assert campaign['status']=='COMPLETE_PAIRED_SPIN0_SUBDIVISION'
 for bits in[256,384]:
  matches=[j for j in campaign['jobs']if j['case']==tag and j['bits']==bits];assert len(matches)==1 and matches[0]['exit_code']==0
  assert matches[0]==i.read(OUT/'spin0_refined'/tag/f'command_l0_{bits}.json')
 i.bind(Path(__file__).with_name('mixture_head_spin0_refined.py'),campaign['source_sha256']);i.bind(Path(__file__).with_name('run_spin0_refinement.py'),campaign['driver_sha256'])
 parts={k:F(v)for k,v in old['components_exact'].items()};oldtotal=sum(parts.values(),F(0));assert oldtotal==F(old['support_upper_exact'])
 oldspin=parts['spin0'];newspin=max(F(r['upper_exact'])for r in pair);assert 0<=newspin<oldspin
 parts['spin0']=newspin;total=sum(parts.values(),F(0));assert total==oldtotal-oldspin+newspin
 norm2=sum(F(v)**2 for v in old['objective_exact']);normlo=F(old['direction_norm_lower_exact'])
 assert norm2==F(old['direction_norm_squared_exact'])and normlo>0 and normlo**2<=norm2<(normlo+F(1,2**256))**2
 upper=upward(total);normal=upward(F(upper)/normlo)
 updates=[]
 for row in old['refined_components']:
  row=dict(row)
  if row['spin']==0:
   row.update(previous_retained_upper_exact=str(oldspin),paired_direct_upper_exact=str(newspin),retained_upper_exact=str(newspin),cells=len(expected),exact_subdivision_factor=2)
  updates.append(row)
 for p in[Path(__file__),ROOT/'research/weak_regions_20261007/assemble_mixture86.py',ROOT/'research/six_targets_20261007/REVIEW_LOWER_SPIN0_REFINEMENT.md']:i.bind(p)
 result=dict(old);result.update(support_upper_exact=str(total),support_upper_outward_9dp=upper,support_normalized_upper_outward_9dp=normal,support_normalized_float=float(normal),inequality=f"({old['objective_exact'][0]}) g0+({old['objective_exact'][1]}) g2 <= {upper}",components_exact={k:str(v)for k,v in parts.items()},components_float={k:float(v)for k,v in parts.items()},refined_components=updates,subdivision_refinement=dict(previous_complete_certificate=str(oldpath.relative_to(ROOT)),previous_complete_sha256=sha(oldpath),old_spin0_upper_exact=str(oldspin),new_spin0_upper_exact=str(newspin),improvement_exact=str(oldspin-newspin),improvement_float=float(oldspin-newspin),method='Exact bisection of every original spin0 mixture ordinary cell; identical gap bounds and field formulas; paired256/384 enclosures.'),sha256_manifest=i.manifest)
 output.write_text(json.dumps(result,indent=2));print(tag,normal,'extra support improvement',float(oldspin-newspin),flush=True)
if __name__=='__main__':
 for tag in['lower15','lower45']:main(tag)
