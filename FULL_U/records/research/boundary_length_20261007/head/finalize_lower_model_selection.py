"""Distinguish a valid-but-loose original cap from an original sign failure."""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'output/boundary_length_20261007/head/angle268'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    path=BASE/'model_selection_final.json'
    saved=BASE/'model_selection_final_initial_metadata.json'
    assert not saved.exists()
    shutil.copyfile(path,saved)
    packet=json.loads(path.read_text()); preflight=json.loads((BASE/'preflight.json').read_text())
    assert preflight['status']=='PASS_PREFLIGHT'
    packet['original_failed_spins']=[]
    packet['replaced_for_quality_spins']=[2]
    packet['replacement_reason']='The original quadratic cap passed sign tests but its support cost exceeded the spin quality target.'
    packet['initial_metadata_record_sha256']=sha(saved)
    packet['dependencies_sha256'][str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    path.write_text(json.dumps(packet,indent=2))
    print(json.dumps(dict(status=packet['status'],final_path=str(path.relative_to(ROOT)),sha256=sha(path))))

if __name__=='__main__':main()
