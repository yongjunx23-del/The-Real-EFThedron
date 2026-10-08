"""Wait for the preserved ordinary producer, then refine and fully replay spin2."""
from prepare_head import *
def run(command,log):
    tic=time.time()
    with log.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=ENV,stdout=f,stderr=subprocess.STDOUT)
    record=dict(command=command,log=str(log.relative_to(ROOT)),exit_code=r.returncode,elapsed_seconds=time.time()-tic)
    assert r.returncode==0,record
    return record
def main():
    out=BASE/'spin2_refined';assert not(out/'commands.json').exists();start=time.time()
    source=BASE/'spin2/head_l2.json'
    while True:
        if source.exists():
            try:ready=json.loads(source.read_text()).get('final',False)
            except (json.JSONDecodeError,PermissionError):ready=False
            if ready:break
        time.sleep(10)
    command=[sys.executable,str(Path(__file__).with_name('refine_spin2_threshold.py'))]
    rows=[run(command,BASE/'threshold_refinement.log')]
    command=[sys.executable,str(OLD/'certify_head_constant.py'),'--witness',str(WITNESS),'--moments',str(BASE/'moments128'),'--out',str(out/'replay_l2.json'),'--replay',str(out/'head_l2.json'),'--bits','384']
    rows.append(run(command,out/'replay_l2.log'))
    (out/'commands.json').write_text(json.dumps(dict(status='PASS_REFINED_SPIN2_AND_FULL_REPLAY',commands=rows,witness_sha256=EXPECTED,script_sha256=sha(__file__),elapsed_seconds_including_wait=time.time()-start),indent=2))
    print('REFINED_SPIN2_AND_REPLAY_COMPLETE',flush=True)
if __name__=='__main__':main()
