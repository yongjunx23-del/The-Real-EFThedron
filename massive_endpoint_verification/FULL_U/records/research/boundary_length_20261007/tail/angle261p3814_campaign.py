"""Frozen parameters and reviewed dependencies for the angle261p3814 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/boundary_length_20261007/search/angle261p3814.npz'
SHA='5bd1b99e5fd4e237fe7e98583aa683ab076d4d76319b05b06c023e583f48383e'
OUT=ROOT/'output/boundary_length_20261007/tail/angle261p3814'
AUDITS=ROOT/'output/boundary_length_20261007/fresh_audit/angle261p3814'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(-5399142818899535,36028797018963968),F(-8905488025140847,9007199254740992))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
