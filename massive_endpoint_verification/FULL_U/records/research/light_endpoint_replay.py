"""Small adapters around unchanged original proof code; no optimization."""
from pathlib import Path
from fractions import Fraction as F
import argparse,importlib.util,json,sys
ROOT=Path(__file__).resolve().parents[1]
def load(relative,name):
 p=ROOT/relative;sys.path.insert(0,str(p.parent));s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def main():
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['c08-finite','c21-sign']);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 if a.mode=='c08-finite':
  module=load('research/full_direction_20261007/combine_dense_family.py','endpoint_dense')
  module.CASES=[(F(0),F(0))]
  module.main()
  result=json.loads((ROOT/'output/full_direction_20261007/certification/dense_finite_family.json').read_text())
 else:
  module=load('research/all_direction_certification_20261007/real220_local_coercivity.py','endpoint_local')
  co=module.Coercivity(ROOT/'output/competitive_20261007/left/real/coercivity_descriptor.json')
  assert co.compact==10 and co.large==3000 and sorted(co.cells)==list(range(42,449,2))
  result={'status':'PASS_SAVED_FINITE_UV_SIGN_PARTITIONS','even_spins':len(co.cells),'spin_range':[min(co.cells),max(co.cells)],'cells':sum(map(len,co.cells.values())),'compact_kappa_exact':str(co.compact),'large_tau_kappa_exact':str(co.large),'scope':'Original Coercivity constructor: paired precision records, exact partitions, negative upper enclosures and all finite spins. Interval fields are retained saved enclosures, not freshly recomputed.'}
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
