"""Final exact arithmetic/hash checks and compact quantitative handoff."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/six_targets_20261007/lower'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
reports=[]
for tag in['lower15','lower45','lower60']:
 p=OUT/f'complete_{tag}_certificate.json';d=json.loads(p.read_text());assert d['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
 for rel,h in d['sha256_manifest'].items():assert sha(ROOT/rel)==h,rel
 definition=json.loads((ROOT/d['witness']).read_text());assert sha(ROOT/d['witness'])==d['witness_sha256']
 assert definition['objective_exact']==d['objective_exact']and definition['weights_exact']==d['weights_exact']
 parts={k:F(v)for k,v in d['components_exact'].items()};total=sum(parts.values(),F(0));assert total==F(d['support_upper_exact'])and min(parts.values())>=0
 assert F(d['retained_other_spins_exact'])==F(d['old_convex_bound_exact'])-F(d['removed_components_exact'])==parts['all_other_spins_and_tails']
 assert [r['spin']for r in d['refined_components']]==([0,2,4]if tag!='lower60'else[0,2])
 for r in d['refined_components']:
  assert F(r['retained_upper_exact'])==min(F(r['paired_direct_upper_exact']),F(r['old_convex_upper_exact']))==parts[f"spin{r['spin']}"]
 norm2=sum(F(v)**2 for v in d['objective_exact']);nl=F(d['direction_norm_lower_exact']);assert nl>0 and nl**2<=norm2<(nl+F(1,2**256))**2
 assert F(d['support_normalized_upper_outward_9dp'])*nl>=total
 new=float(d['support_normalized_upper_outward_9dp']);old=d['target_station']['old_certified_support'];dual=d['graphical_comparison']['dual']['support'];assert new<old
 timings=[]
 for row in d['refined_components']:
  ell=row['spin'];commands=json.loads((OUT/tag/(f'commands_l{ell}_v2.json'if ell==4 else f'commands_l{ell}.json')).read_text());assert len(commands)==2 and all(c['exit_code']==0 for c in commands)
  timings.append(dict(spin=ell,producer_seconds=commands[0]['elapsed_seconds'],replay_seconds=commands[1]['elapsed_seconds'],pair_elapsed_seconds=sum(c['elapsed_seconds']for c in commands)))
 reports.append(dict(tag=tag,packet=str(p.relative_to(ROOT)),packet_sha256=sha(p),witness_sha256=d['witness_sha256'],nominal_angle_degrees=d['nominal_angle_degrees'],old_polygon_support=old,new_certified_support=new,absolute_support_gain=old-new,relative_improvement_percent=100*(old-new)/old,digitized_dual_support=dual,relative_gap_to_digitized_dual_percent=100*(new/dual-1),components_exact=d['components_exact'],components_float={k:float(v)for k,v in parts.items()},refined_components=d['refined_components'],timings=timings,sum_job_elapsed_seconds=sum(x['pair_elapsed_seconds']for x in timings),dependencies_checked=len(d['sha256_manifest'])))
source_names=['prepare.py','mixture_head.py','mixture_head_spin4_v2.py','run_mixtures.py','run_spin4_v2.py','assemble.py','audit_and_report.py']
source_hashes={str(Path(__file__).with_name(n).relative_to(ROOT)):sha(Path(__file__).with_name(n))for n in source_names}
campaigns=[json.loads((OUT/n).read_text())for n in['head_campaign.json','head4_campaign_v2.json']]
report=dict(status='PASS_THREE_LOWER_CERTIFICATES_EXACT_ACCOUNTING_AND_HASH_REVIEW',records=reports,source_sha256=source_hashes,campaign_elapsed_seconds=dict(spins0_2=campaigns[0]['elapsed_seconds'],spin4=campaigns[1]['elapsed_seconds']),timing_note='Per-job elapsed times are summed only as computational work, not serial wall time. The six spin0/2 jobs and two spin4 jobs ran in overlapping campaigns; each job recomputed its full partition at256 then384bits.',scope='All dependency hashes and exact witness/decomposition/component/normalization arithmetic rechecked. Shared interval kernels are not independently rederived here.')
out=OUT/'final_lower_review.json';assert not out.exists();out.write_text(json.dumps(report,indent=2))
lines=['# Completed lower-target results','', 'All three completed continuum certificates passed a final exact arithmetic and dependency-hash audit. A smaller support upper bound is stronger. These mixtures use the existing220-row constraint bank and inherit every omitted-spin and infinite-energy proof from complete parents.','', '| Target station | Normal angle | Previous polygon support | New certified support | Improvement | Digitized dual |','|---|---:|---:|---:|---:|---:|']
for r in reports:lines.append(f"| {r['tag']} | {r['nominal_angle_degrees']:.10f} | {r['old_polygon_support']:.12f} | {r['new_certified_support']:.9f} | {r['relative_improvement_percent']:.4f}% | {r['digitized_dual_support']:.12f} |")
lines+=['','The comparison uses support in exactly the recorded target direction. The parent campaign separately recomputes boundary-length geometry. Published values are digitized graphical comparisons, not reconstructed author certificates.','']
for r in reports:
 lines += [f"## {r['tag']}",'',f"Packet: `{r['packet']}`",f"SHA-256: `{r['packet_sha256']}`",f"Exact witness SHA-256: `{r['witness_sha256']}`",f"Bound dependencies checked: {r['dependencies_checked']}.",'','| Component | Upper contribution |','|---|---:|']
 for k,v in r['components_float'].items():lines.append(f'| {k} | {v:.17g} |')
 lines += ['','| Reintegrated spin | Producer seconds | Replay seconds |','|---|---:|---:|']
 for t in r['timings']:lines.append(f"| {t['spin']} | {t['producer_seconds']:.3f} | {t['replay_seconds']:.3f} |")
 lines += ['',f"Sum of job elapsed times: {r['sum_job_elapsed_seconds']:.3f} seconds. This is not wall time; jobs overlap.",'']
lines += ['## Reproduction and source identity','', 'Read `research/six_targets_20261007/lower/REPRODUCTION.md` for all construction, domain, precision, fallback, and inherited-constraint details. Every rational multiplier and both exact positive weights are stored in the referenced definition. The failed first spin4 attempt remains diagnostic and is excluded from the selected proof records.','',f"Spin0/2 campaign wall time: {campaigns[0]['elapsed_seconds']:.3f} seconds. Separate overlapping spin4 campaign wall time: {campaigns[1]['elapsed_seconds']:.3f} seconds.",'']
for name,h in source_hashes.items():lines.append(f'- `{name}`: `{h}`')
(OUT/'FINAL_LOWER_RESULTS.md').write_text('\n'.join(lines)+'\n')
print(json.dumps([{k:r[k]for k in['tag','packet_sha256','new_certified_support','relative_improvement_percent','relative_gap_to_digitized_dual_percent','sum_job_elapsed_seconds']}for r in reports],indent=2))
