"""Prepare and execute one fresh, relocated all-energy fixed-spin replay.

Uses copies only; the canonical records are not changed. Runtime is roughly a
few minutes. This is a limited executable portability test, not an all-cut replay.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[2]
PACKAGE=ROOT/'output/FIGURE_REPRO_20261008'
WORK=ROOT/'tmp/figure_repro_smoke_20261008'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--assembly-only',action='store_true');a=ap.parse_args()
 assert a.assembly_only or not WORK.exists(),'Use a fresh scratch path; this script never deletes earlier results.'
 records=PACKAGE/'records'
 if not a.assembly_only:
  for folder in ['research','tmp/certification_deps','output/boundary_length_20261007/head/angle268']:
   shutil.copytree(records/folder,WORK/folder)
  for name in ['output/boundary_length_20261007/search/angle268.npz']:
   p=WORK/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(records/name,p)
 evidence=PACKAGE/'verification';evidence.mkdir(exist_ok=True)
 command=[sys.executable,str(PACKAGE/'bundle_run.py'),'--package',str(PACKAGE),'--workdir',str(WORK),'research/boundary_length_20261007/head/certify_head.py','--witness','output/boundary_length_20261007/search/angle268.npz','--moments','output/boundary_length_20261007/head/angle268/moments128','--out','fresh_replay_spin12.json','--replay','output/boundary_length_20261007/head/angle268/spin12/head_l12.json','--bits','384']
 start=time.time()
 if not a.assembly_only:
  p=subprocess.run(command,capture_output=True,text=True)
  (evidence/'fresh_relocated_spin12.log').write_text(p.stdout+p.stderr)
  assert p.returncode==0,p.stderr
 fresh=json.loads((WORK/'fresh_replay_spin12.json').read_text());assert fresh['status']=='PASS_FULL_PARTITION_REPLAY'and fresh['passed']
 shutil.copy2(WORK/'fresh_replay_spin12.json',evidence/'fresh_relocated_spin12.json')
 original=WORK/'output/boundary_length_20261007/head/angle268/head_summary_final.json';prior=json.loads(original.read_text());canonical_hash_before=sha(records/original.relative_to(WORK));assert sha(original)==canonical_hash_before
 assembly=[sys.executable,str(PACKAGE/'bundle_run.py'),'--package',str(PACKAGE),'--workdir',str(WORK),'research/boundary_length_20261007/head/assemble_lower_final.py','--witness','output/boundary_length_20261007/search/angle268.npz','--base','output/boundary_length_20261007/head/angle268']
 p=subprocess.run(assembly,capture_output=True,text=True);(evidence/'fresh_relocated_head_assembly.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
 after=json.loads(original.read_text())
 fields=['status','witness_sha256','direction_exact','included_spins','missing_spins','unfinished_spins','missing_replays','upper_exact','conservative_max_of_original_replay_upper_exact']
 for key in fields:assert prior[key]==after[key],key
 assert sha(records/original.relative_to(WORK))==canonical_hash_before
 source='output/boundary_length_20261007/head/angle268/spin12/head_l12.json'
 report={'status':'PASS_FRESH_RELOCATED_SINGLE_SPIN_REPLAY_AND_HEAD_REASSEMBLY','selected_cut_index':36,'witness':'output/boundary_length_20261007/search/angle268.npz','saved_partition':source,'saved_partition_sha256':sha(records/source),'fresh_replay_source':'research/boundary_length_20261007/head/certify_head.py','fresh_replay_source_sha256':sha(records/'research/boundary_length_20261007/head/certify_head.py'),'precision_bits':fresh['precision_bits'],'cells_recomputed':fresh['cells_replayed'],'gaps_recomputed':fresh['gaps_replayed'],'spin':12,'fixed_spin_all_energy_replay':True,'all_selected_cuts_freshly_replayed':False,'head_assembly_fields_identical':fields,'head_spins_assembled':after['included_spins'],'current_invocation_seconds':time.time()-start,'fresh_interval_seconds':fresh['elapsed_seconds'],'assembly_only_invocation':a.assembly_only,'commands':[command,assembly],'scope':'Fresh directed interval recomputation of exactly one selected fixed-spin saved partition (all energies), followed by exact reassembly of the saved21-spin head. The other20 head integrals and all tails are not freshly interval-recomputed by this test. Source records remained unchanged; only a disposable copy was used.'}
 (evidence/'fresh_relocated_smoke_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
