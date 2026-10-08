"""Continue unchanged proof workers; relax only unstarted quality budgets."""
import os,sys,json,time,hashlib,subprocess,concurrent.futures,ctypes,re
from ctypes import wintypes
from pathlib import Path
from decimal import Decimal
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/peak_certification_20261007/lower263/head'
WITNESS=ROOT/'output/weak_regions_20261007/lower/search/angle263p2747_from261.npz'
EXPECTED='8931658c2a353a2bbf9b9cb074db490ba5d15f044c80837236d3fd9dd801c541'
OLD=ROOT/'research/boundary_length_20261007/head';ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==EXPECTED and not (BASE/'campaign.json').exists()
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.OpenProcess.argtypes=[wintypes.DWORD,wintypes.BOOL,wintypes.DWORD];kernel.OpenProcess.restype=wintypes.HANDLE
kernel.WaitForSingleObject.argtypes=[wintypes.HANDLE,wintypes.DWORD];kernel.WaitForSingleObject.restype=wintypes.DWORD
kernel.GetExitCodeProcess.argtypes=[wintypes.HANDLE,ctypes.POINTER(wintypes.DWORD)];kernel.CloseHandle.argtypes=[wintypes.HANDLE]
switch=json.loads((BASE/'orchestration_switch.json').read_text());handles={};observed={}
for r in switch['active_workers']:
    ell=int(re.search(r'--spin (\d+)',r['CommandLine'])[1]);assert ell in[6,8,10]
    h=kernel.OpenProcess(0x100000|0x1000,False,r['ProcessId']);assert h,('Open worker failed',r,ctypes.get_last_error())
    handles[ell]=h;observed[ell]=r
budget=json.loads((BASE/'budgets.json').read_text());rows=budget['rows'];changes=[]
for r in rows:
    ell=r['spin'];r['max_cells']=7000
    if ell<=10:continue
    assert not (BASE/f'spin{ell}').exists(),'Only unstarted jobs may receive relaxed targets.'
    before=r['target'];r['slack']='.000010';r['target']=str(Decimal(r['numerical'])+Decimal(r['slack']));r['max_cells']=3000
    changes.append(dict(spin=ell,previous_target=before,new_target=r['target'],previous_max_cells=7000,new_max_cells=3000))
budget.update(slack_sum_exact=str(sum((F(r['slack'])for r in rows),F(0))),quality_changes=changes,scope='Quality targets only. All new bounds use actual complete-partition enclosures and their full higher-precision replay. Active original workers were preserved without changed numerical settings.')
(BASE/'budgets_relaxed.json').write_text(json.dumps(budget,indent=2));start=time.time()
def run(command,log):
    tic=time.time()
    with log.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
    entry=dict(command=command,log=str(log.relative_to(ROOT)),exit_code=r.returncode,elapsed_seconds=time.time()-tic)
    assert r.returncode==0,entry
    return entry
def work(row):
    ell=row['spin'];folder=BASE/f'spin{ell}';folder.mkdir(exist_ok=True);saved=folder/'commands.json'
    if saved.exists():
        commands=json.loads(saved.read_text());assert len(commands)==2 and all(c['exit_code']==0 for c in commands)
        return dict(spin=ell,commands=commands)
    command=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'head_l{ell}.json'),'--spin',str(ell),'--model',row['model'],'--target',row['target'],'--max-cells',str(row['max_cells']),'--bits','256']
    if ell in handles:
        h=handles[ell];tic=time.time()
        while kernel.WaitForSingleObject(h,1000)==258:pass
        code=wintypes.DWORD();assert kernel.GetExitCodeProcess(h,ctypes.byref(code));kernel.CloseHandle(h)
        entry=dict(command=command,observed_original_process=observed[ell],exit_code=code.value,log=str((folder/f'head_l{ell}.log').relative_to(ROOT)),observed_wait_seconds=time.time()-tic,adopted_original_worker=True)
        assert code.value==0,entry;commands=[entry]
        d=json.loads((folder/f'head_l{ell}.json').read_text());assert d['final']and d['witness_sha256']==EXPECTED and d['spin']==ell and d['target']==row['target']
    else:commands=[run(command,folder/f'head_l{ell}.log')]
    command=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(folder/f'replay_l{ell}.json'),'--replay',str(folder/f'head_l{ell}.json'),'--bits','384']
    commands.append(run(command,folder/f'replay_l{ell}.log'));saved.write_text(json.dumps(commands,indent=2));print('REPLAYED',ell,flush=True)
    return dict(spin=ell,commands=commands)
with concurrent.futures.ThreadPoolExecutor(max_workers=6)as pool:jobs=list(pool.map(work,rows))
(BASE/'jobs.json').write_text(json.dumps(jobs,indent=2))
command=[sys.executable,str(OLD/'assemble.py'),'--witness',str(WITNESS),'--base',str(BASE)]
assembly=run(command,BASE/'assembly.log');preparation=json.loads((BASE/'preparation.json').read_text())
campaign=dict(status='COMPLETE_REPLAYED_FINITE_HEAD_ONLY',witness_sha256=EXPECTED,jobs=jobs,commands=preparation['commands']+[assembly],elapsed_seconds_after_orchestration_change=time.time()-start,script_sha256=sha(Path(__file__)),orchestration_change=dict(path=str((BASE/'orchestration_switch.json').relative_to(ROOT)),sha256=sha(BASE/'orchestration_switch.json'),old_preparation_and_proof_records_preserved=True,active_workers_adopted=[6,8,10],new_worker_limit=6),budgets_path=str((BASE/'budgets_relaxed.json').relative_to(ROOT)),budgets_sha256=sha(BASE/'budgets_relaxed.json'),scope='Every energy, even spins0..40 only. Same-witness spin42..infinity tail remains separate. Target relaxation affects requested quality, never proof validity.')
(BASE/'campaign.json').write_text(json.dumps(campaign,indent=2));print((BASE/'head_summary.json').read_text(),flush=True)
