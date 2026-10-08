"""Frozen parameters and reviewed dependencies for the angle255 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/full_direction_20261007/left/angle255.npz'
SHA='5b7febef7cb467ff224cd7162a82d76a3534fd3439102f459733e513b4f31615'
OUT=ROOT/'output/boundary_length_20261007/tail/angle255'
AUDITS=ROOT/'output/full_direction_20261007/left/batch_audit/angle255'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(-2331234710160199,9007199254740992),F(-8700286382685973,9007199254740992))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}

CROSSING=ROOT/'output/full_direction_20261007/left/components/angle255_crossing.json'
