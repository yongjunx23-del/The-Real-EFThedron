"""Frozen parameters and reviewed dependencies for the angle259 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/boundary_length_20261007/search/angle259p2941_from255.npz'
SHA='321695c10783e2e0f6e6e28589f132fdb7a68497314b1cf92eb71a07612c704f'
OUT=ROOT/'output/boundary_length_20261007/tail/angle259'
AUDITS=ROOT/'output/boundary_length_20261007/fresh_audit_259/angle259p2941_from255'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(-6692978864176391,36028797018963968),F(-4425208767245823,4503599627370496))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
