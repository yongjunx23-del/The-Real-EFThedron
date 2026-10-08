"""Build explicit per-cut original-assembly and transitive proof-source routes."""
from pathlib import Path
import collections,json
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output/FIGURE_REPRO_20261008'

def main():
 reg=json.loads((OUT/'data/selection_registry.json').read_text());dep=json.loads((OUT/'DEPENDENCIES.json').read_text());known={r['path']for r in dep['files']}
 graph=collections.defaultdict(set)
 for b in dep['record_bindings']:graph[b['parent']].add(b['path'])
 routes=[]
 for row in reg['records']:
  i=row['cut_index'];source=row['source_certificate'];out='fresh_'+row['cut_id']+'.json';stash=[];args=[];script=None;produced=out
  if i in [0,1,2,3,4,5,6,8]:script='research/full_direction_20261007/consolidate_certificates.py';produced=source
  elif i in [15,16,17]:script='research/full_direction_20261007/assemble_free_point.py';produced=source
  elif i==20:script='research/all_direction_certification_20261007/assemble_competitive_g2.py';args=['--strong-uv','--refined-head','--out',out]
  elif i==21:script='research/all_direction_certification_20261007/assemble_competitive_left.py';args=['--head','output/all_direction_certification_20261007/head/left/head_refined_summary.json'];produced=source
  elif i in [29,31,33,34,36,37]:
   angle={29:'255',31:'259',33:'261p3814',34:'265',36:'268',37:'85'}[i];script='research/boundary_length_20261007/assemble_complete_v3.py'
   head='head_summary_final.json'if i==36 else 'head_summary.json'
   args=['--angle',angle,'--head',f'output/boundary_length_20261007/head/angle{angle}/{head}','--finite-banded',f'output/boundary_length_20261007/real/angle{angle}/finite_banded_cost.json','--out',out]
  elif i==39:
   script='research/weak_regions_20261007/assemble_shifted_family_v2.py';args=['--witness','output/weak_regions_20261007/tail_transfer/from268_to269p04041748833.npz','--head-dir','output/weak_regions_20261007/head','--tail','output/weak_regions_20261007/tail_transfer/from268_to269p04041748833_tail_certificate.json','--out',out]
  elif i==41:script='research/weak_regions_20261007/assemble_fresh269p5.py';args=['--out',out]
  elif 43<=i<=50:
   script='research/peak_certification_20261007/assemble_right_family.py';produced=source
   p='output/peak_certification_20261007/right_family/'
   stash=[p+'complete_right_family.json']+[p+f'{tag}_b{b}.json'for b in [20,24,28,40,48,96,192,256]for tag in ['complete','witness']]
  elif i==51:script='research/peak_certification_20261007/assemble_extra_upper.py';args=['--out',out];stash=['output/peak_certification_20261007/extra_upper85p842/mixture_definition.json']
  elif i==52:script='research/peak_certification_20261007/assemble_lower263.py';args=['--out',out]
  elif i==53:script='research/peak_certification_20261007/assemble_upper83_banded.py';args=['--out',out]
  elif i==54:script='research/six_targets_20261007/assemble_upper86_banded_v2.py';args=['--out',out]
  elif i==55:script='research/six_targets_20261007/upper/assemble_upper55.py';args=['--out',out];stash=['output/six_targets_20261007/upper55/mixture_definition.json']
  elif i==56:script='research/six_targets_20261007/assemble_upper95_direct_refined.py';args=['--out',out]
  elif i in [57,58]:
   script='research/six_targets_20261007/lower/assemble_spin0_refined.py';produced=source
   stash=[f'output/six_targets_20261007/lower/complete_lower{k}_refined_certificate.json'for k in [15,45]]
  elif i==59:script='research/six_targets_20261007/lower/assemble_lower60_spin4.py';produced=source;stash=[source]
  elif i!=28:raise AssertionError(i)
  seen=set();todo=[source]if source else []
  while todo:
   p=todo.pop()
   if p in seen:continue
   seen.add(p);todo.extend(graph[p]-seen)
  if script:assert script in known,script
  paths=sorted(seen)
  route={'cut_id':row['cut_id'],'cut_index':i,'source_certificate':source,'source_result_index':row.get('source_result_index'),'proof_family':row['proof_family'],'original_assembly':None if not script else {'script':script,'arguments':args,'produced_record':produced,'stash_existing_derived_outputs':stash,'scope':'Re-executes original exact component assembly/checks from frozen interval records; does not rerun all interval kernels.'},'proof_files':paths,'head_and_replay_records':[p for p in paths if p.endswith('.json')and any(t in p for t in ['/head/','head_spin','refined_spin'])],'source_moment_tables':[p for p in paths if p.endswith('.json')and'moment'in p],'compact_tail_and_finite_records':[p for p in paths if p.endswith('.json')and any(t in p for t in ['/real/','/compact','/tail','finite_','norm_','coercivity'])],'original_proof_code':[p for p in paths if p.endswith('.py')],'original_method_documents':[p for p in paths if p.endswith(('.md','.txt'))],'fresh_all_domain_replay_scope':'To recompute the complete proof, regenerate all moment inputs, every selected head and cap partition, finite-spin sign/mass blocks and infinite-energy/infinite-spin tail records via these preserved producers, then reassemble. The finite-grid optimizer is not a proof step.'}
  routes.append(route)
 commands=[r['path']for r in dep['files']if 'recorded-producer-invocation'in r['reasons']]
 (OUT/'data/proof_replay_routes.json').write_text(json.dumps({'scope':'All39 selected cuts;38 numeric original assemblies plus analytic C28. Transitive source SHA bindings may include shared family ancestors. Numeric assembly routes are executable in disposable copies, with exact output comparison by bundle_reassemble.py. Full fresh interval campaign was not rerun for export.','routes':routes,'recorded_producer_invocations':commands},indent=2)+'\n')
 print('ROUTES',len(routes),'NUMERIC',sum(r['original_assembly']is not None for r in routes),'INVOCATIONS',len(commands))

if __name__=='__main__':main()
