"""Frozen parameters and reviewed dependencies for the anglelower263 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/weak_regions_20261007/lower/search/angle263p2747_from261.npz'
SHA='8931658c2a353a2bbf9b9cb074db490ba5d15f044c80837236d3fd9dd801c541'
OUT=ROOT/'output/peak_certification_20261007/lower263/tail'
AUDITS=ROOT/'output/weak_regions_20261007/lower/audit/angle263p2747_from261'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(-2109641626282615,18014398509481984),F(-8945221711437375,9007199254740992))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
