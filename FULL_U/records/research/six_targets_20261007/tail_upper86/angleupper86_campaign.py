"""Frozen parameters and reviewed dependencies for the angleupper86 tail campaign."""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json,numpy as np
ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'research/all_direction_certification_20261007/tail'
sys.path[:0]=[str(OLD),str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research/tail_certification'),str(ROOT/'research')]
WITNESS=ROOT/'output/weak_regions_20261007/upper/angle86_refined.npz'
SHA='14e2b30df27e95ede876a8ad75ca190904547fb4ace535c14c6a11725d2df97f'
OUT=ROOT/'output/six_targets_20261007/upper86/tail'
AUDITS=ROOT/'output/weak_regions_20261007/upper/audits/angle86_refined'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(WITNESS)==SHA
OBJECTIVE=tuple(F(float(x))for x in np.load(WITNESS)['direction'])
assert OBJECTIVE==(F(8816243287893001,144115188075855872),F(8990329285288823,9007199254740992))

def source_manifest():
 p=OUT/'prepared_components.json';d=json.loads(p.read_text());assert d['witness_sha256']==SHA and d['objective_exact']==list(map(str,OBJECTIVE))
 for q,h in d['sha256_manifest'].items():assert sha(ROOT/q)==h,q
 return dict(d['sha256_manifest'])|{str(p.relative_to(ROOT)):sha(p)}
