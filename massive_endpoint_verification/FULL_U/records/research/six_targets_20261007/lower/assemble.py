"""Exact complete lower certificates by positive mixtures and low-spin substitution."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,math,bisect
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'research/weak_regions_20261007'),str(ROOT/'research')]
from assemble_mixture86 import Inputs,sha,dy,upward,partition
from angle_benchmarks import support_report
OUT=ROOT/'output/six_targets_20261007/lower'

def base_head(i,cert,folder):
 s=i.read(Path(folder)/'head_summary.json');i.dependencies(s)
 assert s['status']=='COMPLETE_REPLAYED_FINITE_HEAD'and s['witness_sha256']==cert['witness_sha256']
 assert s['included_spins']==list(range(0,41,2))and not any(s[k]for k in['missing_spins','unfinished_spins','missing_replays'])
 assert s['direction_exact']==cert['objective_exact'];costs={};records={};paths={}
 for row in s['records']:
  p=i.read(row['path'],row['sha256']);q=i.read(row['replay'],row['replay_sha256']);i.dependencies(p);i.dependencies(q)
  ell=row['spin'];assert p['spin']==q['spin']==ell and p['witness_sha256']==q['witness_sha256']==cert['witness_sha256']
  assert p['status']=='PASS_FINITE_SPIN_ENCLOSURE'and p['final']
  assert q['status']=='PASS_FULL_PARTITION_REPLAY'and q['passed']and q['partition_covers_beta_zero_to_one_exactly']
  assert q['source_record_sha256']==sha(ROOT/row['path'])and q['precision_bits']>p['precision_bits']
  assert p['direction_exact']==q['direction_exact']==cert['objective_exact']
  assert p['source_manifest_sha256']==q['source_manifest_sha256']==s['source_manifest_sha256']
  cells=[(F(c['a']),F(c['b']))for c in p['cells_data']];gaps=[tuple(map(F,g))for g in p['gap_beta_intervals']];partition(cells,gaps)
  assert q['cells_replayed']==len(cells)and q['gaps_replayed']==len(gaps)
  assert dy(p['gap_upper_dyadic'])+sum((dy(c['upper_dyadic'])for c in p['cells_data']),F(0))<=dy(p['upper_dyadic'])
  costs[ell]=max(dy(p['upper_dyadic']),dy(q['upper_dyadic']));records[ell]=p;paths[ell]=row['path']
 assert sum(costs.values(),F(0))==F(cert['components_exact']['head_0_40'])==F(s['conservative_max_of_original_replay_upper_exact'])
 manifest=i.read(Path(folder)/'moments128/manifest.json',s['source_manifest_sha256'])
 assert manifest['witness_sha256']==cert['witness_sha256']
 for table in manifest['tables']:
  for suffix,key in[('_moments.json','moment_sha256'),('_uv_moments.json','uv_sha256')]:
   p=Path(folder)/'moments128'/(table['label']+suffix);data=i.read(p,table[key]);assert data['metadata']['witness_sha256']==cert['witness_sha256']
 return costs,records,paths

def assemble(tag):
 folder=OUT/tag;output=OUT/f'complete_{tag}_certificate.json';assert not output.exists();i=Inputs()
 dp=folder/'definition.json';d=i.read(dp);assert d['witness_format']=='exact_rational_positive_combination'and d['tag']==tag
 cp=d['base_certificate_paths'];certs=[i.read(p,h)for p,h in zip(cp,d['base_certificate_sha256'])]
 for c in certs:
  assert c['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
  assert sum(map(F,c['components_exact'].values()),F(0))==F(c['support_upper_exact'])
 ws=[np.load(i.bind(p,h))for p,h in zip(d['base_witness_paths'],d['base_witness_sha256'])]
 assert all(c['witness_sha256']==h for c,h in zip(certs,d['base_witness_sha256']))
 assert np.array_equal(ws[0]['names'],ws[1]['names'])and ws[0]['names'].tolist()==d['names']and len(d['names'])==220
 weights=list(map(F,d['weights_exact']));direction=list(map(F,d['objective_exact']));assert len(weights)==2 and min(weights)>0
 base_d=[list(map(F,c['objective_exact']))for c in certs]
 assert all([F(float(x))for x in w['direction']]==bd for w,bd in zip(ws,base_d))
 assert [sum(weights[j]*base_d[j][k]for j in[0,1])for k in[0,1]]==direction
 multipliers=[sum(weights[j]*F(float(ws[j]['multipliers'][k]))for j in[0,1])for k in range(220)]
 assert list(map(str,multipliers))==d['multipliers_exact']
 info=[base_head(i,c,p)for c,p in zip(certs,d['base_head_directories'])]
 costs=[r[0]for r in info];records=[r[1]for r in info];paths=[r[2]for r in info]
 assert [{str(k):str(v)for k,v in r.items()}for r in costs]==d['base_retained_costs']
 refined_spins=[0,2,4]if tag in['lower15','lower45']else[0,2]
 updates=[];retained=[]
 for ell in refined_spins:
  pairs=[i.read(folder/f'head_l{ell}_{bits}.json')for bits in[256,384]];geometry=None
  for rec,bits in zip(pairs,[256,384]):
   i.dependencies(rec)
   assert rec['status']=='PASS_EXACT_POSITIVE_MIXTURE_HEAD_COMPONENT'and rec['case']==tag and rec['spin']==ell and rec['precision_bits']==bits
   assert rec['definition_sha256']==sha(dp)and (ROOT/rec['definition_path']).resolve()==dp.resolve()
   assert rec['objective_exact']==d['objective_exact']and rec['weights_exact']==d['weights_exact']and rec['multipliers_exact']==d['multipliers_exact']and rec['names']==d['names']
   assert [str(Path(p))for p in rec['base_certificate_paths']]==[str(Path(p))for p in cp]and rec['base_witness_sha256']==d['base_witness_sha256']
   assert [str(Path(p))for p in rec['source_record_paths']]==[str(Path(r[ell]))for r in paths]
   gaps=[tuple(map(F,g))for g in rec['gap_intervals']]
   assert all(gaps==[tuple(map(F,g))for g in r[ell]['gap_beta_intervals']]for r in records)
   gap_components=list(map(F,rec['gap_components_exact']))
   assert all(gap_components[j]>=dy(records[j][ell]['gap_upper_dyadic'])for j in[0,1])
   assert F(rec['gap_upper_exact'])==sum(weights[j]*gap_components[j]for j in[0,1])
   cells=[(F(c['a']),F(c['b']))for c in rec['cells']];partition(cells,gaps)
   if geometry is None:geometry=cells
   else:assert geometry==cells
   oldcells=[[(F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in r[ell]['cells_data']]for r in records]
   starts=[[c[0]for c in row]for row in oldcells]
   for cell,(lo,hi)in zip(rec['cells'],cells):
    assert F(cell['upper_exact'])>=0 and cell['method']in['convex_base_cells','mixture_riemann','mixture_convex_endpoints','mixture_negative_V']
    if cell['method']=='mixture_negative_V':assert hi<F(1,2)
    if cell['method']=='convex_base_cells':
     fallback=F(0)
     for j in[0,1]:
      index=bisect.bisect_right(starts[j],lo)-1;a,b,value=oldcells[j][index];assert a<=lo<hi<=b;fallback+=weights[j]*value
     assert fallback==F(cell['upper_exact'])
   assert F(rec['upper_exact'])==F(rec['gap_upper_exact'])+sum((F(c['upper_exact'])for c in rec['cells']),F(0))
  direct=max(F(r['upper_exact'])for r in pairs);convex=sum(weights[j]*costs[j][ell]for j in[0,1]);value=min(direct,convex)
  retained.append(value);updates.append(dict(spin=ell,paired_direct_upper_exact=str(direct),old_convex_upper_exact=str(convex),retained_upper_exact=str(value),cells=len(geometry),gaps=len(gaps)))
  commands=i.read(folder/(f'commands_l{ell}_v2.json'if ell==4 else f'commands_l{ell}.json'));assert len(commands)==2 and all(c['exit_code']==0 for c in commands)
 old=sum(weights[j]*F(certs[j]['support_upper_exact'])for j in[0,1]);assert old==F(d['old_convex_total_exact'])
 removed=sum(weights[j]*sum((costs[j][ell]for ell in refined_spins),F(0))for j in[0,1]);remaining=old-removed;assert remaining>=0
 total=remaining+sum(retained,F(0));assert 0<=total<=old
 norm2=sum(x*x for x in direction);scale=2**256;normlo=F(math.isqrt(norm2.numerator*scale*scale//norm2.denominator),scale)
 assert normlo>0 and normlo**2<=norm2<(normlo+F(1,scale))**2
 upper=upward(total);normal=upward(F(upper)/normlo);angle=math.degrees(math.atan2(float(direction[1]),float(direction[0])))%360
 campaign=i.read(OUT/'head_campaign.json');assert campaign['status']=='COMPLETE_PAIRED_PRECISION_MIXTURE_HEADS'
 i.bind(Path(__file__).with_name('mixture_head.py'),campaign['source_sha256']);i.bind(Path(__file__).with_name('run_mixtures.py'),campaign['driver_sha256'])
 if 4 in refined_spins:
  campaign4=i.read(OUT/'head4_campaign_v2.json');assert campaign4['status']=='COMPLETE_PAIRED_PRECISION_MIXTURE_SPIN4'
  i.bind(Path(__file__).with_name('mixture_head_spin4_v2.py'),campaign4['source_sha256']);i.bind(Path(__file__).with_name('run_spin4_v2.py'),campaign4['driver_sha256'])
 for p in[Path(__file__),Path(__file__).with_name('prepare.py'),ROOT/'research/weak_regions_20261007/assemble_mixture86.py',ROOT/'research/six_targets_20261007/REVIEW_LOWER_MIXTURES.md']:i.bind(p)
 components={'all_other_spins_and_tails':remaining,**{f'spin{ell}':v for ell,v in zip(refined_spins,retained)}}
 result=dict(status='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND',witness=str(dp.relative_to(ROOT)),witness_sha256=sha(dp),witness_format=d['witness_format'],objective_exact=d['objective_exact'],nominal_angle_degrees=angle,target_station=d['target_station'],weights_exact=d['weights_exact'],support_upper_exact=str(total),support_upper_outward_9dp=upper,support_normalized_upper_outward_9dp=normal,support_normalized_float=float(normal),direction_norm_squared_exact=str(norm2),direction_norm_lower_exact=str(normlo),inequality=f'({direction[0]}) g0+({direction[1]}) g2 <= {upper}',old_convex_bound_exact=str(old),removed_components_exact=str(removed),retained_other_spins_exact=str(remaining),refined_components=updates,components_exact={k:str(v)for k,v in components.items()},components_float={k:float(v)for k,v in components.items()},coverage=f'Every even spin and every energy above threshold. Spins{refined_spins} are integrated for the exact rational mixture; other head and tail pieces follow from positive homogeneity and subadditivity of disk support.',proof_scope='Same exact positive mixture, paired256/384-bit low-spin interval enclosures, complete exact partitions and exact additive decomposition. Internal computer-assisted certificate under established physical assumptions; not formal verification or external review.',graphical_comparison=dict(dual=support_report(angle),primal=support_report(angle,'primal')),sha256_manifest=i.manifest)
 output.write_text(json.dumps(result,indent=2));print(tag,normal,'gain against old polygon',d['target_station']['old_certified_support']-float(normal),flush=True)
if __name__=='__main__':
 for case in['lower15','lower45','lower60']:assemble(case)
