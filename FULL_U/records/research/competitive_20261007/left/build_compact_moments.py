"""Fresh left-witness contractions of certified compact source moments."""
from pathlib import Path
import sys,time,json,hashlib,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import ctx
from certified_roy_moments import build as roy,save
from certified_strip_moments import build as strips
WITNESS=ROOT/'output/full_direction_20261007/left/angle180_final.npz'
EXPECTED='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=['roy','strips'],required=True)
    ap.add_argument('--bits',type=int,default=512);a=ap.parse_args();ctx.prec=a.bits;tic=time.time()
    assert sha(WITNESS)==EXPECTED
    folder=ROOT/'output/competitive_20261007/left/compact/moments';folder.mkdir(parents=True,exist_ok=True)
    data=[roy(WITNESS,jmax=98,kmax=256,skmax=98)]if a.kind=='roy'else strips(WITNESS,jmax=98,kmax=256)
    records=[]
    for table in data:
        md=table['metadata'];label='roy'if a.kind=='roy'else f"strip_{md['A']}_{md['B']}"
        p=folder/f'{label}_moments.json';save(table,p)
        records.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),metadata=md))
        print('saved',label,'seconds',time.time()-tic,flush=True)
    record=dict(status='PASS_POLYNOMIAL_MOMENT_DEGREE_ENCLOSURES',witness_sha256=EXPECTED,tables=records,seconds=time.time()-tic,
      sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),WITNESS,ROOT/'research/certified_roy_moments.py',ROOT/'research/certified_strip_moments.py',ROOT/'research/validated_quadrature.py']},
      scope='Certified degree-exact polynomial integration with fresh left-witness contractions. Geometric truncation requires separate remainder bounds.')
    (folder/f'{a.kind}_manifest.json').write_text(json.dumps(record,indent=2));print('PASS',a.kind,time.time()-tic,flush=True)
if __name__=='__main__':main()
