"""Build the higher-degree compact source tables, without UV tables."""
from pathlib import Path
import sys,time,json,hashlib,argparse
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import ctx
from certified_roy_moments import build as roy,save
from certified_strip_moments import build as strips


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--kind',choices=['roy','strips'],required=True)
    p.add_argument('--bits',type=int,default=512);a=p.parse_args();ctx.prec=a.bits;tic=time.time()
    witness=ROOT/'output/competitive_20261007/robust/g2_robust_stage1.npz'
    folder=ROOT/'output/competitive_20261007/compact/moments';folder.mkdir(parents=True,exist_ok=True)
    data=[roy(witness,jmax=108,kmax=256,skmax=108)] if a.kind=='roy' else strips(witness,jmax=108,kmax=256)
    records=[]
    for table in data:
        md=table['metadata'];label='roy' if a.kind=='roy' else f"strip_{md['A']}_{md['B']}"
        path=folder/f'{label}_moments.json';save(table,path)
        records.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),metadata=md))
        print('saved',label,'seconds',time.time()-tic,flush=True)
    record=dict(status='PASS_POLYNOMIAL_MOMENT_DEGREE_ENCLOSURES',witness_sha256=sha(witness),tables=records,seconds=time.time()-tic,
        sha256_manifest={str(x.relative_to(ROOT)):sha(x) for x in [Path(__file__),ROOT/'research/certified_roy_moments.py',ROOT/'research/certified_strip_moments.py',ROOT/'research/validated_quadrature.py']},
        scope='Certified Gauss nodes and weights with asserted polynomial exactness degree. Geometric expansions using these moments require separate remainders.')
    (folder/f'{a.kind}_manifest.json').write_text(json.dumps(record,indent=2))
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
