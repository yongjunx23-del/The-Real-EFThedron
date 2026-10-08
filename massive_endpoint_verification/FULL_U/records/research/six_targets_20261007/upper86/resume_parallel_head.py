"""Adopt running proof children; change only scheduler width3->6."""
from prepare_head import *
import ctypes,re,concurrent.futures
from ctypes import wintypes
K=ctypes.WinDLL('kernel32',use_last_error=True)
K.OpenProcess.argtypes=[wintypes.DWORD,wintypes.BOOL,wintypes.DWORD];K.OpenProcess.restype=wintypes.HANDLE
K.WaitForSingleObject.argtypes=[wintypes.HANDLE,wintypes.DWORD];K.WaitForSingleObject.restype=wintypes.DWORD
K.GetExitCodeProcess.argtypes=[wintypes.HANDLE,ctypes.POINTER(wintypes.DWORD)]
K.TerminateProcess.argtypes=[wintypes.HANDLE,wintypes.UINT];K.CloseHandle.argtypes=[wintypes.HANDLE]
def inventory():
 cmd="Get-CimInstance Win32_Process -Filter \"Name = 'python.exe'\" | Select-Object ProcessId,ParentProcessId,CommandLine,CreationDate | ConvertTo-Json -Compress"
 rows=json.loads(subprocess.check_output(['powershell','-NoProfile','-Command',cmd],text=True))
 return rows if isinstance(rows,list)else[rows]
def matched(r):return 'six_targets_20261007' in(r.get('CommandLine')or'')and 'upper86'in r['CommandLine']
start=time.time();assert not(BASE/'campaign.json').exists()and not(BASE/'orchestration_switch.json').exists()
parents=[r for r in inventory()if matched(r)and ('upper86/run_head.py'in r['CommandLine']or'upper86\\run_head.py'in r['CommandLine'])]
assert len(parents)==1;parent=parents[0];handles={};observed={}
def adopt(rows):
 for r in rows:
  if not matched(r)or'certify_head_constant.py'not in r['CommandLine']:continue
  m=re.search(r'spin(\d+)[\\/](?:head|replay)_l',r['CommandLine']);assert m,r
  ell=int(m[1]);phase='replay'if' --replay 'in r['CommandLine']else'producer'
  key=(ell,phase)
  if key in handles:continue
  h=K.OpenProcess(0x100000|0x1000,False,r['ProcessId']);assert h,('Cannot adopt child',r,ctypes.get_last_error())
  handles[key]=h;observed[key]=r
adopt(inventory());assert handles
parent_handle=K.OpenProcess(0x100000|0x1000|1,False,parent['ProcessId']);assert parent_handle
assert K.TerminateProcess(parent_handle,49166);assert K.WaitForSingleObject(parent_handle,10000)==0
exit_code=wintypes.DWORD();assert K.GetExitCodeProcess(parent_handle,ctypes.byref(exit_code));K.CloseHandle(parent_handle);assert exit_code.value==49166
adopt(inventory())
switch=dict(status='OWNED_ORCHESTRATOR_REPLACED_ACTIVE_PROOF_CHILDREN_PRESERVED',old_orchestrator=parent,old_orchestrator_exit_code=exit_code.value,old_orchestrator_sha256=sha(Path(__file__).with_name('run_head.py')),new_orchestrator=str(Path(__file__).relative_to(ROOT)),new_orchestrator_sha256=sha(__file__),active_workers=[dict(r,spin=e,phase=p)for(e,p),r in observed.items()],old_workers=3,new_workers=6,quality_budgets_unchanged=True,budgets_sha256=sha(BASE/'budgets.json'),started_unix=start,scope='Only task-owned scheduler terminated. Proof children retain their original arguments; opened process handles retain genuine exit codes after termination.')
(BASE/'orchestration_switch.json').write_text(json.dumps(switch,indent=2)+'\n');print('ADOPTED',list(handles),flush=True)
rows=json.loads((BASE/'budgets.json').read_text())['rows']
def run(command,log):
 tic=time.time();assert not log.exists()
 with log.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
 result=dict(command=command,log=str(log.relative_to(ROOT)),exit_code=r.returncode,elapsed_seconds=time.time()-tic)
 assert r.returncode==0,result;return result
def waited(key,command,log):
 tic=time.time();h=handles[key]
 while K.WaitForSingleObject(h,1000)==258:pass
 code=wintypes.DWORD();assert K.GetExitCodeProcess(h,ctypes.byref(code));K.CloseHandle(h);handles[key]=None
 row=dict(command=command,log=str(log.relative_to(ROOT)),exit_code=code.value,observed_wait_seconds=time.time()-tic,adopted_original_worker=True,observed_original_process=observed[key])
 assert code.value==0,row;return row
def work(row):
 ell=row['spin'];folder=BASE/f'spin{ell}';saved=folder/'commands.json'
 if saved.exists():
  commands=json.loads(saved.read_text());assert len(commands)==2 and all(c['exit_code']==0 for c in commands);return dict(spin=ell,commands=commands)
 maxcells=6000 if ell<=4 else 2500
 producer=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'head_l{ell}.json'),'--spin',str(ell),'--model',row['model'],'--target',row['target'],'--max-cells',str(maxcells),'--bits','256']
 replay=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'replay_l{ell}.json'),'--replay',str(folder/f'head_l{ell}.json'),'--bits','384']
 if (ell,'replay')in observed:
  d=json.loads((folder/f'head_l{ell}.json').read_text());assert d['final']and d['witness_sha256']==EXPECTED and d['spin']==ell and d['target']==row['target']
  p=dict(command=producer,log=str((folder/f'head_l{ell}.log').relative_to(ROOT)),exit_code=0,exit_code_evidence='The unchanged original scheduler only launches the observed replay after its producer subprocess returns0; the final producer packet independently has PASS status.',elapsed_seconds=d['elapsed_seconds'],completed_before_adoption=True)
  commands=[p,waited((ell,'replay'),replay,folder/f'replay_l{ell}.log')]
 else:
  if (ell,'producer')in observed:p=waited((ell,'producer'),producer,folder/f'head_l{ell}.log')
  else:
   assert not folder.exists(),'Do not overwrite unclassified proof work';folder.mkdir();p=run(producer,folder/f'head_l{ell}.log')
  commands=[p,run(replay,folder/f'replay_l{ell}.log')]
 saved.write_text(json.dumps(commands,indent=2));print('REPLAYED',ell,flush=True);return dict(spin=ell,commands=commands)
with concurrent.futures.ThreadPoolExecutor(max_workers=6)as pool:jobs=list(pool.map(work,rows))
for h in handles.values():
 if h:K.CloseHandle(h)
(BASE/'jobs.json').write_text(json.dumps(jobs,indent=2))
command=[sys.executable,str(OLD/'assemble.py'),'--witness',str(WITNESS),'--base',str(BASE)]
assembly=run(command,BASE/'assembly.log')
campaign=dict(status='COMPLETE_REPLAYED_FINITE_HEAD_ONLY',witness_sha256=EXPECTED,jobs=jobs,assembly=assembly,script_sha256=sha(__file__),elapsed_seconds=time.time()-start,orchestration_change=dict(path=str((BASE/'orchestration_switch.json').relative_to(ROOT)),sha256=sha(BASE/'orchestration_switch.json'),script=str(Path(__file__).relative_to(ROOT))),scope='Every energy, even spins0..40; scheduler width changed, all quality targets and mathematical inputs retained.')
(BASE/'campaign.json').write_text(json.dumps(campaign,indent=2)+'\n');print('COMPLETE',flush=True)
