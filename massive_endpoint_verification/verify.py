#!/usr/bin/env python3
"""Verify three saved continuum certificates, without running an optimizer.
Default checks include exact witnesses, hash dependencies, saved partitions,
finite-spin sign partitions, directed finite-block merging and original
certificate assembly. --fresh-local additionally recomputes C21's finite-spin
real-norm integral at two precisions. This is not a fresh replay of all integrals.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time
ROOT=Path(__file__).resolve().parent
PACKAGE=ROOT/'FULL_U';RECORDS=PACKAGE/'records';IDS={'C08','C20','C21'}
def require(ok,message):
 if not ok:raise RuntimeError(message)
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for block in iter(lambda:f.read(4194304),b''):h.update(block)
 return h.hexdigest()
def dyad(x):return F(int(x[0]))*F(2)**int(x[1])
def relpath(name):
 p=(RECORDS/name.replace('\\','/')).resolve();require(p.is_relative_to(RECORDS),'Unsafe dependency path');return p
def check_hashes(d):
 checked={}
 def visit(v):
  if isinstance(v,dict):
   for k,item in v.items():
    if isinstance(k,str)and isinstance(item,str)and len(item)==64 and all(c in '0123456789abcdef'for c in item)and k.replace('\\','/').startswith(('output/','research/')):
     p=relpath(k);require(p.is_file()and sha(p)==item,'Dependency mismatch: '+k);checked[k.replace('\\','/')]=item
    visit(item)
  elif isinstance(v,list):
   for x in v:visit(x)
 visit(d);return checked
def partition(parts):
 parts=sorted(parts);require(parts and parts[0][0]==0 and parts[-1][1]==1,'Partition endpoints');require(all(a<b for a,b in parts)and all(x[1]==y[0]for x,y in zip(parts,parts[1:])),'Partition gap or overlap')
def exact_witnesses(registry):
 import numpy as np
 data=read(PACKAGE/'data/multipliers_exact.json')['certificates'];require({d['cut_id']for d in data}==IDS,'Wrong multiplier targets');count=0
 for d in data:
  r=registry[d['cut_id']];recipe=r['witness_recipe'];p=relpath(recipe.get('witness',recipe.get('base_witness','output/angles_resummed/angle_000p000.npz')))
  expected=recipe.get('witness_sha256',recipe.get('base_sha256'));require(expected and sha(p)==expected,'Raw witness hash')
  n=np.load(p,allow_pickle=False);require(len(n['multipliers'])==r['bank_row_count']==len(d['rows']),'Multiplier count')
  require([F(float(v))for v in n['direction']]==list(map(F,d['objective_exact'])),'Witness direction')
  for name,value,row in zip(n['names'],n['multipliers'],d['rows']):
   require(str(name)==row['name']and F(float(value))==F(row['exact']),'Exact multiplier mismatch');count+=1
 require(count==590,'Expected 150+220+220 exact multipliers');return {'status':'PASS_EXACT_WITNESSES','multipliers':count}
def saved_partitions():
 reports=[];unsupported=[]
 for p in sorted(RECORDS.rglob('*.json')):
  if not any(x in p.name for x in ['head_l','refined_spin','head_spin']):continue
  d=read(p)
  if d.get('final')is False:continue
  if {'cells_data','gap_beta_intervals','upper_dyadic'}<=d.keys() and ('gap_upper_dyadic'in d or 'gap_upper'in d):
   cells=d['cells_data'];gaps=[tuple(map(F,x))for x in d['gap_beta_intervals']]
   if not all({'a','b','upper_dyadic'}<=c.keys()for c in cells):
    unsupported.append({'path':p.relative_to(RECORDS).as_posix(),'reason':'Different cells_data schema'});continue
   partition([(F(c['a']),F(c['b']))for c in cells]+gaps)
   if 'gap_upper_dyadic'in d:cap=dyad(d['gap_upper_dyadic'])
   else:
    from flint import arb,ctx
    ctx.prec=1024;cap=dyad(arb(d['gap_upper']).upper().man_exp())
   require(cap>=0 and all(dyad(c['upper_dyadic'])>=0 for c in cells),'Negative head cell or cap cost')
   value=cap+sum((dyad(c['upper_dyadic'])for c in cells),F())
   require(value<=dyad(d['upper_dyadic']),'Saved cell/cap sum exceeds record')
   reports.append({'path':p.relative_to(RECORDS).as_posix(),'cells':len(cells),'caps':len(gaps)})
  else:unsupported.append({'path':p.relative_to(RECORDS).as_posix(),'reason':'Different schema; selected C08 direction records use c08_heads and remaining records use original specialized checkers.'})
 require(len(reports)>0,'No supported head partitions checked')
 return {'unsupported_records':unsupported,'status':'PASS_SUPPORTED_HEAD_PARTITIONS','checked_records':len(reports),'ordinary_cells':sum(r['cells']for r in reports),'caps':sum(r['caps']for r in reports),'records':reports,'scope':'Supported cells_data schema only; C08 direction records are checked separately. No interval integrand is reevaluated here.'}
def c08_heads():
 values={}
 for ell in [0,2,4,6]:
  candidates=[];parts=[]
  for suffix in ['','_replay']:
   p=RECORDS/f'output/directional_certification/refined{suffix}_spin{ell}.json';d=read(p);check_hashes(d)
   require(d['status']=='PASS'and d['spin']==ell and len(d['directions'])==1,'C08 head schema')
   require(d['witness_sha256']==sha(RECORDS/'output/angles_resummed/angle_000p000.npz'),'Head witness')
   row=d['directions'][0];require(F(row['slope_exact'])==F(row['crossing_increment_exact'])==0,'Head slope')
   cells=row['cells']+row['gaps'];part=sorted((F(c['a']),F(c['b']))for c in cells);partition(part)
   value=sum((F(c['upper_exact'])for c in cells),F());require(value==F(row['head_upper_exact'])and all(F(c['upper_exact'])>=0 for c in cells),'C08 head cell sums')
   if row['threshold']is not None:require(dyad(row['threshold']['new_V_over_p_upper_dyadic'])<0,'C08 threshold sign')
   candidates.append(value);parts.append(part)
  require(parts[0]==parts[1],'C08 replay partition');values[str(ell)]=max(candidates)
 return values

def release_integrity():
 candidates=[ROOT/'MANIFEST.json',ROOT.parent/'MANIFEST.json']
 manifests=[p for p in candidates if p.is_file()]
 if not manifests:return {'status':'UNCHECKED_RELEASE_MANIFEST','reason':'No release manifest in staging. Delivered archive should include the enclosing manifest.'}
 reports=[]
 for path in manifests:
  doc=read(path);seen=set()
  for row in doc['files']:
   name=row['path'];file=(path.parent/name).resolve()
   require(file.is_relative_to(path.parent)and file!=path.parent and name not in seen,'Unsafe release manifest path');seen.add(name)
   require(file.is_file()and sha(file)==row['sha256'],'Release manifest mismatch: '+name)
   if 'bytes'in row:require(file.stat().st_size==row['bytes'],'Release byte count')
  if 'file_count'in doc:require(len(seen)==doc['file_count'],'Release count mismatch')
  reports.append({'manifest':str(path),'files':len(seen)})
 return {'status':'PASS_RELEASE_MANIFEST','manifests':reports}

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True,help='New report directory outside this extracted package');ap.add_argument('--fresh-local',action='store_true',help='Also recompute C21 finite42..448 weighted norm at256/384bits');a=ap.parse_args()
 require(not sys.flags.optimize,'Do not run Python with -O: original code uses assertions.');out=a.out.resolve();require(not out.exists()and not out.is_relative_to(ROOT),'Output must be a new directory outside verification/')
 import flint,numpy
 out.mkdir(parents=True);start=time.time();result={'status':'RUNNING','scope':'Saved certificate checks and assembly; no optimization; not a fresh complete interval replay.','python':sys.version,'python_flint':getattr(flint,'__version__','?')}
 registry={r['cut_id']:r for r in read(PACKAGE/'data/selection_registry.json')['records']};require(set(registry)==IDS,'Endpoint registry')
 result['release_integrity']=release_integrity()
 result['exact_witnesses']=exact_witnesses(registry);print('PASS exact590 multipliers',flush=True)
 # Verify original byte identities for every retained record in the dependency manifest.
 rows=read(PACKAGE/'DEPENDENCIES.json')['files'];checked=0
 for row in rows:
  p=relpath(row['path']);require(p.is_file()and sha(p)==row['sha256'],'Original dependency bytes: '+row['path']);checked+=1
 result['dependency_files']=checked;result['saved_partitions']=saved_partitions();write(out/'saved_partitions.json',result['saved_partitions']);print('PASS saved head partitions',flush=True)
 heads=c08_heads();result['c08_heads']={k:str(v)for k,v in heads.items()}
 with tempfile.TemporaryDirectory(prefix='endpoint-work-',dir=out)as temp:
  work=Path(temp)/'records';shutil.copytree(RECORDS,work)
  def run(script,args,log):
   command=[sys.executable,'-B',str(PACKAGE/'bundle_run.py'),'--package',str(PACKAGE),'--workdir',str(work),script]+list(map(str,args))
   env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
   with (out/log).open('w')as stream:ret=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,env=env)
   require(ret.returncode==0,'Original checker failed: '+log);return command
  # Original exact finite-block merge, restricting its CASES list to b=0.
  run('research/light_endpoint_replay.py',['c08-finite','--out','endpoint_C08_finite.json'],'C08_finite.log')
  finite=read(work/'endpoint_C08_finite.json');frow=finite['results'][0];require(frow['b']=='0'and frow['delta_X4']=='0'and frow['status']=='PASS','C08 finite direction')
  old=read(RECORDS/'output/full_direction_20261007/certification/dense_finite_family.json');orow=next(r for r in old['results']if r['b']=='0');require(F(frow['finite_support_upper_exact'])==F(orow['finite_support_upper_exact']),'C08 directed finite merge differs')
  shutil.copy2(work/'endpoint_C08_finite.json',out/'C08_finite_merge.json')
  # Original all-energy and infinite-spin tail assembler checks its coercivity margin.
  tailpath=work/'output/certification_tail/fixed_tail_certificate.json';original=tailpath.read_bytes()
  run('research/tail_certification/verify_fixed_tail.py',['--assemble-only'],'C08_tail.log')
  tail=read(tailpath);shutil.copy2(tailpath,out/'C08_tail_assembly.json');tailpath.write_bytes(original)
  saved=read(RECORDS/'output/directional_certification/tail_family_4b.json');check_hashes(saved);trow=next(r for r in saved['results']if r['b']=='0')
  require(F(tail['disk_support_tail_upper_exact'])==F(trow['tail_support_upper_exact']),'C08 tail mismatch')
  total=sum(heads.values(),F())+F(frow['finite_support_upper_exact'])+F(trow['tail_support_upper_exact'])
  results={'C08':{'objective_exact':['1','0'],'support_upper_exact':str(total)}}
  run('research/light_endpoint_replay.py',['c21-sign','--out','endpoint_C21_sign.json'],'C21_sign.log')
  shutil.copy2(work/'endpoint_C21_sign.json',out/'C21_saved_sign_partitions.json');print('PASS C21 saved finite-spin sign partitions',flush=True)
  # Recheck C20 saved sign leaves and continuum domains with original assemblers.
  # Restore their output bytes before subsequent manifest checks.
  sign_routes=[
   ('research/competitive_20261007/robust/assemble_compact_coercivity.py','output/competitive_20261007/compact/complete_compact_coercivity.json','C20_compact_sign'),
   ('research/all_direction_certification_20261007/tail/assemble_robust_uv.py','output/all_direction_certification_20261007/tail/g2_robust_stage1/complete_uv_coercivity.json','C20_base_UV_sign'),
   ('research/all_direction_certification_20261007/tail/assemble_strong_uv.py','output/all_direction_certification_20261007/tail/g2_robust_stage1/strong_uv_coercivity.json','C20_strong_UV_sign')]
  for script,name,label in sign_routes:
   target=work/name;frozen=target.read_bytes()
   run(script,[],label+'.log');shutil.copy2(target,out/(label+'.json'));target.write_bytes(frozen)
  print('PASS C20 compact and UV saved sign partitions',flush=True)
  run('research/all_direction_certification_20261007/assemble_competitive_g2.py',['--strong-uv','--refined-head','--out','endpoint_C20.json'],'C20_assembly.log');results['C20']=read(work/'endpoint_C20.json')
  run('research/all_direction_certification_20261007/assemble_competitive_left.py',['--head','output/all_direction_certification_20261007/head/left/head_refined_summary.json','--finite-cost','output/competitive_20261007/left/real/local_cost2048.json'],'C21_assembly.log');results['C21']=read(work/'output/competitive_20261007/left/complete_left_certificate.json')
  for cut,d in results.items():
   require(list(map(F,d['objective_exact']))==list(map(F,registry[cut]['objective_exact'])),'Objective mismatch')
   require(F(d['support_upper_exact'])==F(registry[cut]['exact_certified_support']),'Exact support mismatch for '+cut)
   require(F(d['support_upper_exact'])<=F(registry[cut]['plotted_support_exact']),'Incorrect outward rounding')
   write(out/(cut+'_assembly.json'),d);print('PASS',cut,{'C08':'g0 <= 2.705955','C20':'g2 <= 1.052351','C21':'g0 >= -7.751502'}[cut],flush=True)
  if a.fresh_local:
   print('Recomputing C21 finite42..448 real norm at256/384bits',flush=True)
   run('research/all_direction_certification_20261007/real220_local_coercivity.py',['--descriptor','output/competitive_20261007/left/real/coercivity_descriptor.json','--cells','2048','--bits','256','--replay-bits','384','--out','endpoint_C21_fresh_local.json'],'C21_fresh_local.log')
   fresh=read(work/'endpoint_C21_fresh_local.json');old=read(RECORDS/'output/competitive_20261007/left/real/local_cost2048.json')
   require(fresh['status']=='PASS_LOCAL_COERCIVITY_WEIGHTED_REAL_COST','Fresh local status')
   require(F(fresh['cost_upper_exact'])==F(old['cost_upper_exact']),'Fresh local upper differs; inspect and reassemble explicitly')
   shutil.copy2(work/'endpoint_C21_fresh_local.json',out/'C21_fresh_local.json');result['fresh_local']={'status':'PASS_FRESH_LOCAL_EXACT_MATCH','scope':'Only C21 finite42..448 real norm recomputed, using saved coercivity proofs.'}
 result['status']='PASS_THREE_ENDPOINT_CERTIFICATES';result['bounds']={c:{'objective_exact':registry[c]['objective_exact'],'exact_support':registry[c]['exact_certified_support'],'outward_support':registry[c]['plotted_support_exact']}for c in sorted(IDS)};result['elapsed_seconds']=time.time()-start;write(out/'result.json',result);print('Reports:',out)
if __name__=='__main__':main()
