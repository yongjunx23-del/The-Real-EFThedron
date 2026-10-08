"""Audit final selected lower certificates, preserving baseline reviews unchanged."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/six_targets_20261007/lower'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
baseline=json.loads((OUT/'final_lower_review.json').read_text());assert baseline['status'].startswith('PASS')
records=[]
for tag in['lower15','lower45','lower60']:
 p=OUT/f'complete_{tag}{"_refined"if tag!="lower60"else""}_certificate.json';c=json.loads(p.read_text());assert c['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
 for rel,h in c['sha256_manifest'].items():assert sha(ROOT/rel)==h,rel
 parts={k:F(v)for k,v in c['components_exact'].items()};assert min(parts.values())>=0 and sum(parts.values(),F(0))==F(c['support_upper_exact'])
 norm2=sum(F(v)**2 for v in c['objective_exact']);normlo=F(c['direction_norm_lower_exact']);assert normlo>0 and normlo**2<=norm2<(normlo+F(1,2**256))**2
 assert F(c['support_normalized_upper_outward_9dp'])*normlo>=F(c['support_upper_exact'])
 base=next(r for r in baseline['records']if r['tag']==tag);added=[]
 if tag!='lower60':
  old=json.loads((ROOT/base['packet']).read_text());ref=c['subdivision_refinement'];assert ref['previous_complete_sha256']==sha(ROOT/base['packet'])
  assert F(old['components_exact']['spin0'])-F(c['components_exact']['spin0'])==F(ref['improvement_exact'])>0
  assert F(old['support_upper_exact'])-F(c['support_upper_exact'])==F(ref['improvement_exact'])
  assert all(old['components_exact'][k]==v for k,v in c['components_exact'].items()if k!='spin0')
  for bits in[256,384]:
   job=json.loads((OUT/'spin0_refined'/tag/f'command_l0_{bits}.json').read_text());assert job['exit_code']==0;added.append(dict(bits=bits,elapsed_seconds=job['elapsed_seconds']))
 new=float(c['support_normalized_upper_outward_9dp']);before=base['old_polygon_support'];dual=c['graphical_comparison']['dual']['support']
 assert new<=base['new_certified_support']<before
 records.append(dict(tag=tag,packet=str(p.relative_to(ROOT)),packet_sha256=sha(p),witness=str(c['witness']),witness_sha256=c['witness_sha256'],angle=c['nominal_angle_degrees'],old_polygon_support=before,initial_new_support=base['new_certified_support'],selected_certified_support=new,absolute_support_gain=before-new,relative_improvement_percent=100*(before-new)/before,digitized_dual_support=dual,relative_gap_to_digitized_dual_percent=100*(new/dual-1),components_exact=c['components_exact'],components_float={k:float(v)for k,v in parts.items()},initial_jobs=base['timings'],refinement_jobs=added,sum_job_elapsed_seconds=base['sum_job_elapsed_seconds']+sum(x['elapsed_seconds']for x in added),dependencies_checked=len(c['sha256_manifest'])))
sources={str(Path(__file__).with_name(n).relative_to(ROOT)):sha(Path(__file__).with_name(n))for n in['mixture_head.py','mixture_head_spin4_v2.py','mixture_head_spin0_refined.py','assemble.py','assemble_spin0_refined.py','audit_selected.py']}
r=dict(status='PASS_SELECTED_LOWER_CERTIFICATES_EXACT_ACCOUNTING_AND_HASH_REVIEW',records=records,source_sha256=sources,baseline_review=str((OUT/'final_lower_review.json').relative_to(ROOT)),baseline_review_sha256=sha(OUT/'final_lower_review.json'),timing_note='Summed per-job elapsed seconds measure computational work, not wall time; independent jobs and campaigns overlap. No new220-multiplier numerical optimization was used.')
out=OUT/'final_selected_lower_review.json';assert not out.exists();out.write_text(json.dumps(r,indent=2))
lines=['# Final selected lower-target certificates','', '| Boundary station | Normal angle | Previous support | Selected certified support | Improvement | Digitized dual |','|---|---:|---:|---:|---:|---:|']
for d in records:lines.append(f"| {d['tag']} | {d['angle']:.10f} | {d['old_polygon_support']:.12f} | {d['selected_certified_support']:.9f} | {d['relative_improvement_percent']:.4f}% | {d['digitized_dual_support']:.12f} |")
lines+=['','Smaller support is stronger. These are exact positive mixtures with interval reintegration, not independently reoptimized220-component witnesses. Published comparisons use digitized curves. Full continuum coverage and its assumptions are inherited from the authenticated parents plus the selected interval head records.','']
for d in records:
 lines+=[f"## {d['tag']}",'',f"Packet: `{d['packet']}`",f"SHA-256: `{d['packet_sha256']}`",f"Witness SHA-256: `{d['witness_sha256']}`",f"Verified bound dependencies: {d['dependencies_checked']}.",'','| Component | Retained upper contribution |','|---|---:|']
 for k,v in d['components_float'].items():lines.append(f'| {k} | {v:.17g} |')
 lines+=['','| Job | Producer seconds | Replay seconds |','|---|---:|---:|']
 for t in d['initial_jobs']:lines.append(f"| Initial spin{t['spin']} | {t['producer_seconds']:.3f} | {t['replay_seconds']:.3f} |")
 if d['refinement_jobs']:lines.append(f"| Extra spin0 bisection | {d['refinement_jobs'][0]['elapsed_seconds']:.3f} | {d['refinement_jobs'][1]['elapsed_seconds']:.3f} |")
 lines+=['',f"Sum of overlapping job elapsed times: {d['sum_job_elapsed_seconds']:.3f} seconds. This is not elapsed wall time.",'']
lines+=['## Source hashes','']+[f'- `{p}`: `{h}`'for p,h in sources.items()]
lines+=['','The exact definitions contain all multipliers, row identities and weights. See `research/six_targets_20261007/lower/REPRODUCTION.md` for complete commands, constraints, approximations, source windows, partition handling and proof limitations. Baseline complete packets and their audit remain unchanged.','']
(OUT/'FINAL_SELECTED_LOWER_RESULTS.md').write_text('\n'.join(lines));print(json.dumps(records,indent=2))
