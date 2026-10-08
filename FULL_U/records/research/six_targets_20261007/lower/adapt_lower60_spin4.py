"""Separate lower60 spin4 adapter with conservative threshold-gap union."""
from pathlib import Path
import json,hashlib
BASE=Path(__file__).parent;ROOT=BASE.parents[2]
old=BASE/'mixture_head_spin4_v2.py';new=BASE/'mixture_head_lower60_spin4.py'
s=old.read_text()
replacements=[
 ("choices=['lower15','lower45']","choices=['lower60']"),
 ("    assert records[0]['gap_beta_intervals']==records[1]['gap_beta_intervals']\n    gaps=[tuple(map(F,x))for x in records[0]['gap_beta_intervals']]", "    original_gaps=[[tuple(map(F,x))for x in r['gap_beta_intervals']]for r in records]\n    assert original_gaps[0][1:]==original_gaps[1][1:]\n    assert original_gaps[0][0][0]==original_gaps[1][0][0]==0\n    gaps=[(F(0),max(g[0][1]for g in original_gaps))]+original_gaps[0][1:]"),
 ("    gap=sum(w*g for w,g in zip(weights,gap_components));total=gap\n    oldcells=[[ (F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in r['cells_data']]for r in records]", "    original_gap_components=list(gap_components)\n    oldcells=[[ (F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in r['cells_data']]for r in records]\n    charged_cells=[]\n    for j in[0,1]:\n        charged=[(lo,hi,v)for lo,hi,v in oldcells[j]if any(max(lo,ga)<min(hi,gb)for ga,gb in gaps)]\n        gap_components[j]+=sum((v for lo,hi,v in charged),F(0))\n        charged_cells.append([dict(a=str(lo),b=str(hi),upper_exact=str(v))for lo,hi,v in charged])\n    gap=sum(w*g for w,g in zip(weights,gap_components));total=gap"),
 ("gap_intervals=records[0]['gap_beta_intervals'],gap_components_exact", "gap_intervals=[[str(lo),str(hi)]for lo,hi in gaps],gap_coarsening=dict(original_gap_intervals=[r['gap_beta_intervals']for r in records],original_gap_components_exact=list(map(str,original_gap_components)),charged_whole_cells=charged_cells,rule='Original analytic gaps plus full nonnegative cost of every ordinary base cell whose interior intersects the enlarged union; whole cells are never scaled by overlap.'),gap_components_exact")]
for a,b in replacements:assert s.count(a)==1;s=s.replace(a,b)
assert not new.exists();new.write_text(s)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=ROOT/'output/six_targets_20261007/lower/lower60_spin4';out.mkdir(exist_ok=True)
(out/'adaptation.json').write_text(json.dumps(dict(original_source=str(old.relative_to(ROOT)),original_sha256=sha(old),source=str(new.relative_to(ROOT)),source_sha256=sha(new),replacements=replacements,scope='All field formulas unchanged. Only lower60 case selection and exact threshold-gap coarsening change; retained ordinary geometry is the complement of the union.'),indent=2))
print(sha(new))
