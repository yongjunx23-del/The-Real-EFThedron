"""Fresh angle85-witness contractions of certified compact source moments."""
from pathlib import Path
import sys,time,json,hashlib,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import ctx
from certified_roy_moments import build as roy,save
from certified_strip_moments import build as strips
WITNESS=ROOT/'output/full_direction_20261007/right/angle85_220.npz'
EXPECTED='58fafda2eafa3a741610175f416a009c84f9a18f309d36e74d0e4007f65fd6d8'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=['roy','strips'],required=True)
    ap.add_argument('--bits',type=int,default=512);a=ap.parse_args();ctx.prec=a.bits;tic=time.time()
    assert sha(WITNESS)==EXPECTED
    folder=ROOT/'output/boundary_length_20261007/compact/angle85/moments';folder.mkdir(parents=True,exist_ok=True)
    data=[roy(WITNESS,jmax=110,kmax=256,skmax=110)]if a.kind=='roy'else strips(WITNESS,jmax=110,kmax=256)
    records=[]
    for table in data:
        md=table['metadata'];label='roy'if a.kind=='roy'else f"strip_{md['A']}_{md['B']}"
        p=folder/f'{label}_moments.json';save(table,p)
        records.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),metadata=md))
        print('saved',label,'seconds',time.time()-tic,flush=True)
    record=dict(status='PASS_POLYNOMIAL_MOMENT_DEGREE_ENCLOSURES',witness_sha256=EXPECTED,tables=records,seconds=time.time()-tic,
      sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),WITNESS,ROOT/'research/certified_roy_moments.py',ROOT/'research/certified_strip_moments.py',ROOT/'research/validated_quadrature.py']},
      scope='Certified degree-exact polynomial integration with fresh angle85-witness contractions. Geometric truncation requires separate remainder bounds.')
    (folder/f'{a.kind}_manifest.json').write_text(json.dumps(record,indent=2));print('PASS',a.kind,time.time()-tic,flush=True)
if __name__=='__main__':main()
