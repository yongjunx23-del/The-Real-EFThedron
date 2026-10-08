"""Independent identity/coverage/accounting review of the threshold-cap refinement."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'output/peak_certification_20261007/upper83/head';hashes={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    d=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
    for key in ['dependencies','dependency_sha256']:
        for rel,h in d.get(key,{}).items():assert sha(ROOT/rel)==h;hashes[rel]=h
    return d
def dy(v):return F(int(v[0]))*F(2)**int(v[1])
old=read(BASE/'spin2/head_l2.json');new=read(BASE/'spin2_refined/head_l2.json')
assert old['final']and new['final']and old['model']=='centered'and new['model']=='constant_centered'
for key in ['witness_sha256','direction_exact','spin','N','precision_bits','source_manifest_sha256']:assert old[key]==new[key]
assert new['threshold_refinement_source_sha256']==sha(BASE/'spin2/head_l2.json')
assert new['spin']==2 and new['precision_bits']==256 and new['new_threshold_cells']==64
oldgaps=[tuple(map(F,x))for x in old['gap_beta_intervals']];newgaps=[tuple(map(F,x))for x in new['gap_beta_intervals']]
assert newgaps[1:]==oldgaps[1:]and newgaps[0]==(F(0),F(1,1024))and oldgaps[0]==(F(0),F(1,256))
oldcells={(c['a'],c['b']):c for c in old['cells_data']};newcells={(c['a'],c['b']):c for c in new['cells_data']}
assert len(oldcells)==len(old['cells_data'])and len(newcells)==len(new['cells_data'])
for key,c in oldcells.items():assert newcells[key]==c
added=[c for key,c in newcells.items()if key not in oldcells];assert len(added)==64
added.sort(key=lambda c:F(c['a']));assert F(added[0]['a'])==F(1,1024)and F(added[-1]['b'])==F(1,256)
parts=sorted([(F(c['a']),F(c['b']))for c in new['cells_data']]+newgaps)
assert parts[0][0]==0 and parts[-1][1]==1 and all(a<b for a,b in parts)and all(x[1]==y[0]for x,y in zip(parts,parts[1:]))
assert dy(new['gap_upper_dyadic'])+sum((dy(c['upper_dyadic'])for c in new['cells_data']),F(0))<=dy(new['upper_dyadic'])
assert all(dy(c['upper_dyadic'])>=0 for c in new['cells_data'])and dy(new['gap_upper_dyadic'])>=0
hashes[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
report=dict(status='PASS_THRESHOLD_REFINEMENT_IDENTITY_COVERAGE_AND_ACCOUNTING_REVIEW',witness_sha256=new['witness_sha256'],old_upper_exact=str(dy(old['upper_dyadic'])),new_upper_exact=str(dy(new['upper_dyadic'])),retained_old_cells=len(oldcells),new_cells=len(added),new_threshold_beta_upper_exact='1/1024',old_threshold_beta_upper_exact='1/256',sha256_manifest=hashes,source_review='The wrapper changes the existing threshold cap only. ConstantCentered inherits every ordinary-cell method from Centered and delegates its cap to the unchanged FiniteSpin constant-margin proof. The newly exposed interval is integrated, old ordinary upper enclosures copied exactly, other gaps preserved, and all components summed outward. No widened normalized source branch is introduced.',remaining='The final complete head must use the larger of this producer and its mandatory full384-bit replay, instead of the old spin2 pair. This review alone is not a complete head or directional certificate.')
(ROOT/'output/peak_certification_20261007/upper83/cap_refinement_peer_review.json').write_text(json.dumps(report,indent=2));print(report['status'],float(dy(old['upper_dyadic'])),float(dy(new['upper_dyadic'])))
