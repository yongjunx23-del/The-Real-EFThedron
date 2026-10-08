"""Rigorous polynomial source tables for a specified fixed witness."""
import argparse,hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import ctx
from certified_roy_moments import build as roy,save,load
from certified_strip_moments import build as strips
from certified_uv_moments import build as uv,save as save_uv
ap=argparse.ArgumentParser();ap.add_argument('--witness',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
ap.add_argument('--lmax',type=int,default=40);ap.add_argument('--N',type=int,default=64);ap.add_argument('--Nuv',type=int,default=32);ap.add_argument('--bits',type=int,default=256)
a=ap.parse_args();ctx.prec=a.bits;a.out.mkdir(parents=True,exist_ok=True);tic=time.time();rec=[]
jmax=a.lmax+a.Nuv+2
for label,d in [('roy',roy(a.witness,jmax,a.N,a.N))]+[(f"strip_{x['metadata']['A']}_{x['metadata']['B']}",x) for x in strips(a.witness,jmax,a.N)]:
    path=a.out/(label+'_moments.json');save(d,path)
    ud=uv(d,a.lmax,a.Nuv);ud['metadata']['source_table_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    up=a.out/(label+'_uv_moments.json');save_uv(ud,up)
    rec.append(dict(label=label,moment_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),uv_sha256=hashlib.sha256(up.read_bytes()).hexdigest(),moment_metadata=d['metadata'],uv_metadata=ud['metadata']))
    print('DONE',label,time.time()-tic,flush=True)
r=dict(witness_sha256=hashlib.sha256(a.witness.read_bytes()).hexdigest(),seconds=time.time()-tic,tables=rec)
(a.out/'manifest.json').write_text(json.dumps(r,indent=2))
