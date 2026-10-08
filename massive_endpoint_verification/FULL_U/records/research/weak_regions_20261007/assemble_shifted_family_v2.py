"""Retain unchanged assembly checks, substituting validated repaired spin2."""
from pathlib import Path
import json
import assemble_shifted_family as original
from repair_shifted_spin2 import corrected_record,sha,ROOT
CHECK=original.check_head
class Inputs(original.Inputs):
    def __init__(self):
        super().__init__()
        self.bind(Path(__file__))
        self.bind(ROOT/'research/weak_regions_20261007/repair_shifted_spin2.py')
        self.read(ROOT/'output/weak_regions_20261007/head_corrected/CORRECTION_NOTICE.json')
original.Inputs=Inputs
def check_head(inputs,path,ell,bits,wp,base,direction,delta):
    path=Path(path).resolve()
    if ell==2:
        raw=inputs.read(path);source=inputs.read(raw['source_record'],raw['source_record_sha256'])
        repaired=ROOT/'output/weak_regions_20261007/head_corrected'/f'head_l2_{bits}.json'
        d=inputs.read(repaired)
        assert d==corrected_record(raw,source,path),'Corrected record must equal deterministic reconstruction.'
        assert d['correction']['status']=='PASS_DETERMINISTIC_MONOTONICITY_REPAIR'
        inputs.bind(d['correction']['repair_script'],d['correction']['repair_script_sha256'])
        return CHECK(inputs,repaired,ell,bits,wp,base,direction,delta)
    return CHECK(inputs,path,ell,bits,wp,base,direction,delta)
original.check_head=check_head
if __name__=='__main__':original.main()
