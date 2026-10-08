"""Frozen parameters and reviewed dependencies for the angleupper83 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/weak_regions_20261007/upper/angle83_repaired.npz'
SHA='68e39591f3f9a9d3b5cdfe3ade3924fecb227a329de4a9543729f0607e8f663a'
OUT=ROOT/'output/peak_certification_20261007/upper83/tail'
AUDITS=ROOT/'output/weak_regions_20261007/upper/repaired83_audit/audits/angle83_repaired'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(535119985632953,4503599627370496),F(8943390007067311,9007199254740992))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
