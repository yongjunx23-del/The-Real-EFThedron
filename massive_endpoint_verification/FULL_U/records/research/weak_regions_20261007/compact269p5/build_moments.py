"""New exact-witness compact contractions; original builders remain frozen."""
from pathlib import Path
import sys,time,json,hashlib,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import ctx
from certified_roy_moments import build as roy,save
from certified_strip_moments import build as strips
SPECS={'269p5':(60,'f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae')}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--angle',choices=list(SPECS),required=True)
    ap.add_argument('--kind',choices=['roy','strips'],required=True);ap.add_argument('--bits',type=int,default=512)
    a=ap.parse_args();ctx.prec=a.bits;tic=time.time();jmax,expected=SPECS[a.angle]
    witness=ROOT/'output/weak_regions_20261007/lower/search/angle269p5.npz'
    assert sha(witness)==expected
    folder=ROOT/'output/weak_regions_20261007/fresh269p5/compact/moments';folder.mkdir(parents=True,exist_ok=True)
    assert not (folder/f'{a.kind}_manifest.json').exists()
    data=[roy(witness,jmax=jmax,kmax=256,skmax=jmax)]if a.kind=='roy'else strips(witness,jmax=jmax,kmax=256)
    records=[]
    for table in data:
        md=table['metadata'];label='roy'if a.kind=='roy'else f"strip_{md['A']}_{md['B']}"
        path=folder/f'{label}_moments.json';assert not path.exists();save(table,path)
        records.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),metadata=md))
        print('saved',a.angle,label,time.time()-tic,flush=True)
    record=dict(status='PASS_POLYNOMIAL_MOMENT_DEGREE_ENCLOSURES',witness_sha256=expected,tables=records,seconds=time.time()-tic,
      sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),witness,ROOT/'research/certified_roy_moments.py',ROOT/'research/certified_strip_moments.py',ROOT/'research/validated_quadrature.py']},
      scope='Certified degree-exact polynomial integration with fresh witness contractions. Geometric truncation requires separate remainder bounds.')
    (folder/f'{a.kind}_manifest.json').write_text(json.dumps(record,indent=2));print('PASS',a.angle,a.kind,time.time()-tic,flush=True)
if __name__=='__main__':main()
