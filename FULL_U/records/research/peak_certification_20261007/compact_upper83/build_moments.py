"""Fresh degree-exact contractions for the analytic plan's finite complement."""
from pathlib import Path
import sys,time,json,hashlib,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import ctx
from certified_roy_moments import build as roy,save
from certified_strip_moments import build as strips
WITNESS=ROOT/'output/weak_regions_20261007/upper/angle83_repaired.npz'
SHA='68e39591f3f9a9d3b5cdfe3ade3924fecb227a329de4a9543729f0607e8f663a'
OUT=ROOT/'output/peak_certification_20261007/upper83/compact'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=['roy','strips'],required=True);ap.add_argument('--bits',type=int,default=512)
    a=ap.parse_args();ctx.prec=a.bits;tic=time.time();assert sha(WITNESS)==SHA
    plan_path=OUT/'analytic_piecewise_512.json';plan=json.loads(plan_path.read_text())
    assert plan['witness_sha256']==SHA and plan['status']=='PASS_PIECEWISE_ALLSPIN_COMPACT_ENVELOPES'
    jmax=max(60,max(r['minimum_even_spin']-2 for r in plan['records']))
    folder=OUT/'moments';folder.mkdir(parents=True,exist_ok=True);assert not(folder/f'{a.kind}_manifest.json').exists()
    data=[roy(WITNESS,jmax=jmax,kmax=256,skmax=jmax)]if a.kind=='roy'else strips(WITNESS,jmax=jmax,kmax=256)
    records=[]
    for table in data:
        md=table['metadata'];label='roy'if a.kind=='roy'else f"strip_{md['A']}_{md['B']}"
        path=folder/f'{label}_moments.json';assert not path.exists();save(table,path)
        records.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),metadata=md));print('SAVED',label,time.time()-tic,flush=True)
    deps=[Path(__file__),WITNESS,plan_path,ROOT/'research/certified_roy_moments.py',ROOT/'research/certified_strip_moments.py',ROOT/'research/validated_quadrature.py']
    record=dict(status='PASS_POLYNOMIAL_MOMENT_DEGREE_ENCLOSURES',witness_sha256=SHA,tables=records,seconds=time.time()-tic,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in deps},scope='Degree-exact integration with new witness contractions; separate geometric remainder bounds are still required.')
    (folder/f'{a.kind}_manifest.json').write_text(json.dumps(record,indent=2));print('PASS',a.kind,jmax,time.time()-tic,flush=True)
if __name__=='__main__':main()
