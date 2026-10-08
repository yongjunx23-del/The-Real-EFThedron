"""Replace only authenticated convex spin4 cost in completed lower60 packet."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,bisect
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(Path(__file__).parent))
from assemble import base_head,Inputs,sha,dy,upward,partition
OUT=ROOT/'output/six_targets_20261007/lower';EXTRA=OUT/'lower60_spin4'
def main():
 output=OUT/'complete_lower60_refined_certificate.json';assert not output.exists();i=Inputs()
 oldpath=OUT/'complete_lower60_certificate.json';old=i.read(oldpath)
 assert old['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
 dp=ROOT/old['witness'];d=i.read(dp,old['witness_sha256']);weights=list(map(F,d['weights_exact']))
 assert min(weights)>0 and d['tag']=='lower60'and d['objective_exact']==old['objective_exact']
 certs=[i.read(p,h)for p,h in zip(d['base_certificate_paths'],d['base_certificate_sha256'])]
 info=[base_head(i,c,p)for c,p in zip(certs,d['base_head_directories'])]
 costs=[x[0]for x in info];records=[x[1][4]for x in info];paths=[x[2][4]for x in info]
 assert [{str(k):str(v)for k,v in c.items()}for c in costs]==d['base_retained_costs']
 oldgaps=[[tuple(map(F,g))for g in r['gap_beta_intervals']]for r in records]
 assert oldgaps[0][1:]==oldgaps[1][1:]and all(g[0][0]==0 for g in oldgaps)
 gaps=[(F(0),max(g[0][1]for g in oldgaps))]+oldgaps[0][1:]
 oldcells=[[(F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in r['cells_data']]for r in records]
 starts=[[x[0]for x in row]for row in oldcells]
 charged=[[(a,b,v)for a,b,v in row if any(max(a,ga)<min(b,gb)for ga,gb in gaps)]for row in oldcells]
 original_gap_cost=[dy(r['gap_upper_dyadic'])for r in records]
 gapcost=[original_gap_cost[j]+sum((v for a,b,v in charged[j]),F(0))for j in[0,1]]
 knots=sorted(set([x for row in oldcells for a,b,v in row for x in[a,b]]+[x for g in gaps for x in g]))
 expected=[(a,b)for a,b in zip(knots,knots[1:])if not any(a>=ga and b<=gb for ga,gb in gaps)]
 partition(expected,gaps)
 pair=[i.read(EXTRA/f'head_l4_{bits}.json')for bits in[256,384]]
 for r,bits in zip(pair,[256,384]):
  i.dependencies(r)
  assert r['status']=='PASS_EXACT_POSITIVE_MIXTURE_HEAD_COMPONENT'and r['case']=='lower60'and r['spin']==4 and r['precision_bits']==bits
  assert r['definition_sha256']==sha(dp)and(ROOT/r['definition_path']).resolve()==dp.resolve()
  for key in['objective_exact','weights_exact','multipliers_exact','names','base_witness_sha256']:assert r[key]==d[key],key
  assert [str(Path(p))for p in r['base_certificate_paths']]==[str(Path(p))for p in d['base_certificate_paths']]
  assert [str(Path(p))for p in r['source_record_paths']]==[str(Path(p))for p in paths]
  assert [tuple(map(F,g))for g in r['gap_intervals']]==gaps
  coarsen=r['gap_coarsening'];assert coarsen['original_gap_intervals']==[rec['gap_beta_intervals']for rec in records]
  assert list(map(F,coarsen['original_gap_components_exact']))==original_gap_cost
  assert coarsen['charged_whole_cells']==[[dict(a=str(a),b=str(b),upper_exact=str(v))for a,b,v in row]for row in charged]
  assert list(map(F,r['gap_components_exact']))==gapcost
  assert F(r['gap_upper_exact'])==sum(weights[j]*gapcost[j]for j in[0,1])
  assert [(F(c['a']),F(c['b']))for c in r['cells']]==expected
  for c in r['cells']:
   a,b=F(c['a']),F(c['b']);v=F(c['upper_exact']);assert v>=0
   assert c['method']in['convex_base_cells','mixture_riemann','mixture_convex_endpoints','mixture_negative_V']
   if c['method']=='mixture_negative_V':assert b<F(1,2)
   if c['method']=='convex_base_cells':
    value=F(0)
    for j in[0,1]:
     ind=bisect.bisect_right(starts[j],a)-1;lo,hi,cost=oldcells[j][ind];assert lo<=a<b<=hi;value+=weights[j]*cost
    assert value==v
  assert F(r['upper_exact'])==F(r['gap_upper_exact'])+sum((F(c['upper_exact'])for c in r['cells']),F(0))
 campaign=i.read(EXTRA/'campaign.json');assert campaign['status']=='COMPLETE_PAIRED_LOWER60_SPIN4_UNION_GAPS'
 for bits in[256,384]:
  matches=[j for j in campaign['jobs']if j['bits']==bits];assert len(matches)==1 and matches[0]['case']=='lower60'and matches[0]['spin']==4 and matches[0]['exit_code']==0
  assert matches[0]==i.read(EXTRA/f'command_l4_{bits}.json')
 source=Path(__file__).with_name('mixture_head_lower60_spin4.py');i.bind(source,campaign['source_sha256']);i.bind(Path(__file__).with_name('run_lower60_spin4.py'),campaign['driver_sha256'])
 adaptation=i.read(EXTRA/'adaptation.json');i.bind(adaptation['original_source'],adaptation['original_sha256']);i.bind(adaptation['source'],adaptation['source_sha256'])
 text=(ROOT/adaptation['original_source']).read_text()
 for a,b in adaptation['replacements']:assert text.count(a)==1;text=text.replace(a,b)
 assert text==source.read_text()
 parts={k:F(v)for k,v in old['components_exact'].items()};assert set(parts)=={'all_other_spins_and_tails','spin0','spin2'}
 assert [r['spin']for r in old['refined_components']]==[0,2]
 original_convex=sum(weights[j]*F(certs[j]['support_upper_exact'])for j in[0,1])
 original_removed=sum(weights[j]*(costs[j][0]+costs[j][2])for j in[0,1])
 assert F(old['old_convex_bound_exact'])==original_convex and F(old['removed_components_exact'])==original_removed
 assert parts['all_other_spins_and_tails']==F(old['retained_other_spins_exact'])==original_convex-original_removed
 previous=sum(parts.values(),F(0));assert previous==F(old['support_upper_exact'])
 convex=sum(weights[j]*costs[j][4]for j in[0,1]);direct=max(F(r['upper_exact'])for r in pair);assert 0<=direct<convex
 parts['all_other_spins_and_tails']-=convex;assert parts['all_other_spins_and_tails']>=0;parts['spin4']=direct
 total=sum(parts.values(),F(0));assert total==previous-convex+direct
 norm2=sum(F(x)**2 for x in old['objective_exact']);normlo=F(old['direction_norm_lower_exact'])
 assert norm2==F(old['direction_norm_squared_exact'])and 0<normlo and normlo**2<=norm2<(normlo+F(1,2**256))**2
 upper=upward(total);normal=upward(F(upper)/normlo)
 for p in[Path(__file__),Path(__file__).with_name('assemble.py'),Path(__file__).with_name('adapt_lower60_spin4.py'),ROOT/'research/six_targets_20261007/REVIEW_LOWER60_SPIN4.md']:i.bind(p)
 result=dict(old);result.update(support_upper_exact=str(total),support_upper_outward_9dp=upper,support_normalized_upper_outward_9dp=normal,support_normalized_float=float(normal),inequality=f"({old['objective_exact'][0]}) g0+({old['objective_exact'][1]}) g2 <= {upper}",components_exact={k:str(v)for k,v in parts.items()},components_float={k:float(v)for k,v in parts.items()},removed_components_exact=str(F(old['removed_components_exact'])+convex),retained_other_spins_exact=str(parts['all_other_spins_and_tails']),refined_components=old['refined_components']+[dict(spin=4,paired_direct_upper_exact=str(direct),old_convex_upper_exact=str(convex),retained_upper_exact=str(direct),cells=len(expected),gaps=len(gaps))],spin4_refinement=dict(previous_complete_certificate=str(oldpath.relative_to(ROOT)),previous_complete_sha256=sha(oldpath),old_spin4_convex_upper_exact=str(convex),new_spin4_paired_upper_exact=str(direct),improvement_exact=str(convex-direct),improvement_float=float(convex-direct),coarsened_gap_whole_cell_counts=list(map(len,charged)),method='Exact threshold-gap union with complete original gaps plus every intersecting whole ordinary base cell; unchanged spin4 interval formulas at256/384 bits.'),coverage='Every even spin and every energy above threshold. Spins0,2,4 are integrated for the same exact rational mixture; other head/tail pieces retain the authenticated positive-homogeneous inherited costs.',sha256_manifest=i.manifest)
 output.write_text(json.dumps(result,indent=2));print(normal,'extra spin4 gain',float(convex-direct),flush=True)
if __name__=='__main__':main()
