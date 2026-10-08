"""Frozen parameters and reviewed dependencies for the angleupper95 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/six_targets_20261007/upper95_direct/witness.npz'
SHA='767b9d4aa9247f505004c392031f60ed82a3aedc7b1354e7a21ee9657e7a24ec'
OUT=ROOT/'output/six_targets_20261007/upper95_direct/tail'
AUDITS=ROOT/'output/six_targets_20261007/upper95_direct/audits'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(1203580983497459,9007199254740992),F(4463211602388183,4503599627370496))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
