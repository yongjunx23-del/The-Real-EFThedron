"""Rebuild every rational moment and remainder bound without modifying cache."""
from finite_source_corrections import *

def main():
 tic=time.time();p=OUT/'finite_source_N12_K80.pkl';old=pickle.loads(p.read_bytes());fresh=build(old['N'],old['K'])
 keys=['witness_sha256','N','K','moments','remainder_coefficients','remainder_by_source','definition','denominator_lower','proof','builder_sha256']
 for key in keys:assert old[key]==fresh[key],key
 files=[Path(__file__),Path(build.__code__.co_filename),p,p.with_suffix('.json'),WITNESS]
 report=dict(status='PASS_EXACT_RATIONAL_REBUILD',moment_count=len(old['moments'])*len(old['moments'][0]),remainder_coefficient_count=sum(len(g['coefficient_abs_upper_exact'])for g in old['remainder_by_source']),all_entries_equal=True,verified_keys=keys,seconds=time.time()-tic,sha256_manifest={str(q.relative_to(ROOT)):sha(q)for q in files})
 (OUT/'finite_source_exact_rebuild.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

if __name__=='__main__':main()
