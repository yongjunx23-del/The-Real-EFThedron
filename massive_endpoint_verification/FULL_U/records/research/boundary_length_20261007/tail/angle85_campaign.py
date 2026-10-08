"""Frozen parameters and reviewed dependencies for the angle85 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/full_direction_20261007/right/angle85_220.npz'
SHA='58fafda2eafa3a741610175f416a009c84f9a18f309d36e74d0e4007f65fd6d8'
OUT=ROOT/'output/boundary_length_20261007/tail/angle85'
AUDITS=ROOT/'output/full_direction_20261007/right/fresh_candidate_audits/angle85_220'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(24532160660097,281474976710656),F(4486462071114449,4503599627370496))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
