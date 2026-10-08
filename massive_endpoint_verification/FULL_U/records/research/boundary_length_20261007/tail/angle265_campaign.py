"""Frozen parameters and reviewed dependencies for the angle265 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/boundary_length_20261007/search/angle265.npz'
SHA='5aae4902c2db47f5e3a5f5a37d87c136cd7ec3cb191752d8f3b2f7ee04d6d8c3'
OUT=ROOT/'output/boundary_length_20261007/tail/angle265'
AUDITS=ROOT/'output/boundary_length_20261007/fresh_audit/angle265'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(-785029141123105,9007199254740992),F(-4486462071114449,4503599627370496))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
