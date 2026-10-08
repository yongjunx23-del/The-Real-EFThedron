"""Frozen parameters and reviewed dependencies for the angle269p5 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/weak_regions_20261007/lower/search/angle269p5.npz'
SHA='f7f5f929ee8885d8fa9f795912405c063c31bf5a6df7d09a71be5562569f46ae'
OUT=ROOT/'output/weak_regions_20261007/fresh269p5/tail'
AUDITS=ROOT/'output/weak_regions_20261007/lower/audit/angle269p5'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(-2515252609197481,288230376151711744),F(-9006856288192973,9007199254740992))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
