"""Apply the reviewed unbounded-tau theorem to the fixed left witness."""
import json,sys
from pathlib import Path
import large_tau_highspin220 as base
ROOT=base.ROOT;OUT=ROOT/'output/all_direction_certification_20261007/tail/left/refined'
WITNESS=ROOT/'output/full_direction_20261007/left/angle180_final.npz'
def main():
 base.OUT=OUT;base.WITNESS=WITNESS;base.main()
 bits=int(sys.argv[sys.argv.index('--precision')+1])if'--precision'in sys.argv else 320
 L=int(sys.argv[sys.argv.index('--lmin')+1])if'--lmin'in sys.argv else 450
 p=OUT/f'large_tau_L{L}_{bits}.json';d=json.loads(p.read_text());d['objective_exact']=['-1','0'];d['conservative_objective_offset_exact']='1';d['note']='The reused absolute forward estimate adds1 while the left objective b=0, so this is a conservative upper envelope.';d['sha256_manifest'][str(Path(__file__).relative_to(ROOT))]=base.sha(__file__);p.write_text(json.dumps(d,indent=2))
if __name__=='__main__':main()
