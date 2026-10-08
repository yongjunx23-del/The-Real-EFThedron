"""Exact saved-cell accounting for supported original all-energy head formats.

No kernels are reevaluated. Unsupported head formats are enumerated explicitly.
This checks every included record with the specific cells_data/gap schema, which
includes authenticated ancestors as well as selected records; it is not a claim
to verify all proof formats or all physical assumptions.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,time

def dy(x):return F(int(x[0]))*F(2)**int(x[1])

def check(package):
 start=time.time();m=json.loads((package/'DEPENDENCIES.json').read_text());rows=[];unsupported=[];notfinal=[]
 for row in m['files']:
  name=row['path'];base=Path(name).name
  if not name.startswith('output/')or not name.endswith('.json'):continue
  if not any(s in base for s in ['head_l','head_spin','refined_spin','wide_head_spin','half_head_spin']):continue
  d=json.loads((package/'records'/name).read_text(encoding='utf-8-sig'))
  needed={'cells_data','gap_beta_intervals','gap_upper_dyadic','upper_dyadic'}
  if not needed<=d.keys():
   unsupported.append({'path':name,'status':d.get('status'),'reason':'Different record schema; this checker does not validate it.'});continue
  if d.get('final')is False:
   notfinal.append(name);continue
  cells=d['cells_data'];gaps=[tuple(map(F,p))for p in d['gap_beta_intervals']]
  if not all({'a','b','upper_dyadic'}<=c.keys()for c in cells):
   unsupported.append({'path':name,'status':d.get('status'),'reason':'Unsupported cell schema.'});continue
  ordinary=[(F(c['a']),F(c['b']))for c in cells];parts=sorted(ordinary+gaps)
  assert parts and parts[0][0]==0 and parts[-1][1]==1,name
  assert all(a<b for a,b in parts),name
  assert all(a[1]==b[0]for a,b in zip(parts,parts[1:])),name
  assert d.get('cells',len(cells))==len(cells),name
  total=dy(d['gap_upper_dyadic'])+sum((dy(c['upper_dyadic'])for c in cells),F())
  assert total<=dy(d['upper_dyadic']),name
  rows.append({'path':name,'spin':d.get('spin'),'ordinary_cells':len(cells),'gap_intervals':len(gaps),'exact_partition_0_1':True,'saved_cell_and_cap_sum_le_declared_upper':True,'slack_exact':str(dy(d['upper_dyadic'])-total),'source_sha256':row['sha256']})
 return {'status':'PASS_SUPPORTED_SAVED_HEAD_PARTITIONS_AND_EXACT_SUMS','checked_records':len(rows),'ordinary_cells':sum(r['ordinary_cells']for r in rows),'gap_intervals':sum(r['gap_intervals']for r in rows),'records':rows,'unsupported_records':unsupported,'nonfinal_records_not_checked':notfinal,'seconds':time.time()-start,'scope':'Exact partitions and sums for the explicitly supported saved-record schema only. Kernel values, derivative bounds, caps and omitted-spin tail inequalities are not freshly recomputed. Unsupported records are listed and must use their original specialized verifier.'}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--report',type=Path);a=ap.parse_args();r=check(a.package.resolve())
 if a.report:a.report.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps({k:v for k,v in r.items()if k not in ['records','unsupported_records','nonfinal_records_not_checked']},indent=2));print('UNSUPPORTED',len(r['unsupported_records']),'NONFINAL',len(r['nonfinal_records_not_checked']))
