"""Peer checks of upper95 units, normalized window and unequal-gap accounting."""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/six_targets_20261007/upper'))
import upper95_head as h
from upper95_common import canonical_records,union_gap_data,source_paths,OUT,WITNESSES
F=h.F;ctx=h.ctx;ctx.prec=320
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
checks=[];gaps_report=[]
for ell in[0,2]:
 raw=[json.loads((ROOT/p).read_text())for p in source_paths(ell)]
 cells,gaps,oldbounds=canonical_records(raw,ell);ordinary,union,expanded,ledger=union_gap_data(cells,gaps,oldbounds)
 # Independently classify the common refinement by midpoint membership.
 knots=sorted({F(0),F(1)}|{v for gs in gaps for pair in gs for v in pair}|{v for cc in cells for a,b,u in cc for v in[a,b]})
 independent=[(a,b)for a,b in zip(knots,knots[1:])if not any(lo<(a+b)/2<hi for gs in gaps for lo,hi in gs)]
 assert independent==ordinary
 for j in[0,1]:
  used=[k for k,(a,b,u)in enumerate(cells[j])if any(max(a,lo)<min(b,hi)for lo,hi in union)]
  assert used==ledger[j]['intersecting_cell_indices']
  assert expanded[j]==oldbounds[j]+sum((cells[j][k][2]for k in used),F(0))
 first=h.FiniteCentered(ROOT/WITNESSES[0],ell,ROOT/'output/peak_certification_20261007/upper83/head/moments128',N=128)
 old=h.SpinZero(ROOT/WITNESSES[1])if ell==0 else h.LowSpin(ROOT/WITNESSES[1],ell)
 second=(h.DirectionalZero if ell==0 else h.DirectionalLow)(old,F(6),F(107,5))
 if ell==2:
  for model in[first,second]:
   assert all(s['A']==4 and s['B']==F(119,10)for s in model.sources if s['kind']=='roy')
   for beta in[F(1,10),F(1,2),F(7,10),F(4,5),F(81,100)]:
    u,v,_=model.fields(beta,beta);nu,nv,pinv,e=model.normalized_interval(beta,beta)
    assert u.is_finite()and v.is_finite()and nu.is_finite()and nv.is_finite()
    assert u.overlaps(nu)and(v*pinv*e*e).overlaps(nv)
    checks.append(dict(spin=ell,model=type(model).__name__,beta=str(beta),U_units_overlap=True,V_scaled_units_overlap=True))
 for beta in[F(1,10),F(7,10),F(9,10),F(99,100)]:
  u,v,_=second.fields(beta,beta);ou,ov,_=old.fields(beta,beta)
  assert u.overlaps(ou)and(v-ov).overlaps(second.addition(beta,beta))
  checks.append(dict(spin=ell,model=type(second).__name__,beta=str(beta),exact_directional_addition_overlap=True))
 gaps_report.append(dict(spin=ell,ordinary_cells=len(ordinary),union_gaps=len(union),whole_cells_added=[len(row['intersecting_cell_indices'])for row in ledger],independent_gap_complement_and_cost_check=True))
paths=[Path(__file__),ROOT/'research/six_targets_20261007/upper/upper95_head.py',ROOT/'research/six_targets_20261007/upper/upper95_common.py',ROOT/'research/full_direction_20261007/head_family_extended.py',ROOT/'research/all_direction_certification_20261007/head/finite_centered.py',ROOT/'research/head_certification_lowspin.py',ROOT/'research/head_certification_spin0.py']
r=dict(status='PASS_UPPER95_PEER_GAP_AND_FIELD_UNIT_CHECKS',checks=checks,gap_checks=gaps_report,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in paths},scope='Independent exact gap complement/accounting reconstruction and numerical interval overlap diagnostics at selected points supporting source-level unit review. These point diagnostics do not establish continuum validity on their own; complete producer/replay bounds remain required.')
p=OUT/'peer_field_and_gap_checks.json';assert not p.exists();p.write_text(json.dumps(r,indent=2));print(r['status'],len(checks),flush=True)
