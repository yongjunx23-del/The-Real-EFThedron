import os,sys,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/weak_regions_20261007/upper';sys.path[:0]=[str(ROOT/'research'),str(ROOT/'research/lower_g0_upper_g2_20261006')]
from matrix_crossing_series import coefficient_polynomial,shifted_certificate,replay_identity
ref=OUT/'angle83_refined.npz';base=OUT/'angle83p1759644845.npz';x=np.load(ref);y=np.load(base);assert np.array_equal(x['direction'],y['direction'])and np.array_equal(x['names'],y['names']);records=[]
for w in [1e-8,1e-7,1e-6,1e-5,1e-4,.001,.01,.05,.1]:
    path=OUT/f'angle83_repair_w{w:g}.npz';lam=(1-w)*x['multipliers']+w*y['multipliers'];np.savez_compressed(path,names=x['names'],multipliers=lam,direction=x['direction']);poly=coefficient_polynomial(path);cert=shifted_certificate(poly,10,1);rec={'old83_weight':w,'witness':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'crossing_pass':cert is not None};records.append(rec);print(rec,flush=True)
    if cert is not None:
        rec['certificate']=cert;rec['identity_replay']=replay_identity(poly,path);break
(OUT/'repair83_search.json').write_text(json.dumps({'records':records,'parents_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in[ref,base]},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'claim':'Exact crossing component repair only; other components must be freshly audited.'},indent=2))
