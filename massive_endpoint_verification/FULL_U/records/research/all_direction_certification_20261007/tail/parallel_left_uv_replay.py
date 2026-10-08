"""Fresh directed-interval replay of one frozen left finite-UV partition."""
from concurrent.futures import ProcessPoolExecutor
from finite_uv_bridge220 import *
import finite_uv_bridge220 as uv
from uniform_source_uv import LEFT,FOLDER

def initialize(bits):
    global BANK
    ctx.prec=bits;uv.WITNESS=LEFT;uv.OUT=FOLDER;BANK=uv.FiniteUV()

def replay_spin(record):
    tic=time.time();l=record['spin'];cells=[];cursor=F(0)
    for c in record['cells']:
        a,b=F(c['a']),F(c['b']);assert cursor==a<b
        v=BANK.field(l,a,b)-1
        assert v.is_finite()and v.upper()<-ar(F(1,10000)),(l,a,b,str(v))
        cells.append(dict(a=str(a),b=str(b),upper_exact=str(upper(v))));cursor=b
    assert cursor==min(F(256),F(l*(l+1),45))==F(record['tau_end_exact'])
    out=dict(spin=l,status='PASS',tau_end_exact=str(cursor),cells=cells,evaluations=len(cells),seconds=time.time()-tic)
    print(json.dumps({k:v for k,v in out.items()if k!='cells'}|dict(cell_count=len(cells))),flush=True)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('reference');ap.add_argument('--precision',type=int,default=480);ap.add_argument('--workers',type=int,default=6);a=ap.parse_args();tic=time.time()
    path=Path(a.reference).resolve();d=json.loads(path.read_text())
    assert d['status']=='PASS_LEFT_FINITE_UV'and not d['failed']and d['witness_sha256']==sha(LEFT)and d['objective_exact']==['-1','0']
    assert F(d['kappa_exact'])==F(1,10000)and a.precision>d['precision_bits']
    assert d['spins']==[r['spin']for r in d['records']]
    for p,h in d['sha256_manifest'].items():assert sha(ROOT/p)==h,p
    with ProcessPoolExecutor(max_workers=a.workers,initializer=initialize,initargs=(a.precision,))as pool:
        records=list(pool.map(replay_spin,d['records']))
    manifest=dict(d['sha256_manifest'])
    for p in[Path(__file__),Path(uv.__file__),path]:manifest[str(p.relative_to(ROOT))]=sha(p)
    for p,h in manifest.items():assert sha(ROOT/p)==h,p
    out=dict(status='PASS_LEFT_FINITE_UV',witness_sha256=sha(LEFT),objective_exact=['-1','0'],spins=d['spins'],kappa_exact='1/10000',precision_bits=a.precision,records=records,failed=[],reference_certificate_sha256=sha(path),acceptance_policy='Every original closed interval replayed independently, finite and directed upper<-1/10000. No adaptive changes or midpoint acceptance.',crossing_omitted_kmin=64,crossing_tail_Bernstein_coefficients=d['crossing_tail_Bernstein_coefficients'],seconds=time.time()-tic,sha256_manifest=manifest)
    target=FOLDER/f"left_finite_uv_{min(d['spins'])}_{max(d['spins'])}_{a.precision}.json";assert not target.exists();target.write_text(json.dumps(out,indent=2))
    print(json.dumps(dict(status=out['status'],file=str(target.relative_to(ROOT)),seconds=out['seconds'])),flush=True)

if __name__=='__main__':main()
