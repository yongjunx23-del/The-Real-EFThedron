"""Run one selected original assembler in a disposable copy and compare exactly.

Before moving any existing derived output, record its SHA256 and preserve it as
.before_reassembly. Input records are never deleted. Use a fresh copy per route;
several C indices share the same family producer. This is saved-component
assembly, NOT fresh interval reevaluation of all head and tail kernels.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,subprocess,sys,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def safe(root,name):
 p=(root/name).resolve();assert p.is_relative_to(root)and p!=root;return p
def objective(d):return list(map(F,d.get('objective_exact',d.get('direction_exact',['1',d.get('slope_exact','0')]))))
def support(d):return F(d.get('support_upper_exact',d.get('certified_upper_exact','0')))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--workdir',type=Path,required=True);ap.add_argument('--cut',required=True);a=ap.parse_args()
 package=a.package.resolve();work=a.workdir.resolve();assert work.is_dir()and work!=package/'records'and not work.is_relative_to(package)
 routes=json.loads((package/'data/proof_replay_routes.json').read_text())['routes'];r=next(x for x in routes if x['cut_id']==a.cut)
 if r['original_assembly']is None:
  assert a.cut=='C28';print('C28 is the analytic inequality g2>=0; no numerical assembly.');return
 route=r['original_assembly'];registry=json.loads((package/'data/selection_registry.json').read_text());expected=next(x for x in registry['records']if x['cut_id']==a.cut)
 transcript=work/('reassembly_'+a.cut+'.json');assert not transcript.exists();moves=[]
 for name in route['stash_existing_derived_outputs']:
  p=safe(work,name)
  if not p.exists():moves.append({'path':name,'existed':False});continue
  dest=p.with_name(p.name+'.before_reassembly');assert p.is_file()and not dest.exists()
  h=sha(p);p.rename(dest);assert sha(dest)==h;moves.append({'path':name,'existed':True,'sha256_before_move':h,'preserved_as':dest.relative_to(work).as_posix()})
 command=[sys.executable,str(package/'bundle_run.py'),'--package',str(package),'--workdir',str(work),route['script']]+route['arguments']
 start=time.time();p=subprocess.run(command,capture_output=True,text=True)
 (work/('reassembly_'+a.cut+'.log')).write_text(p.stdout+p.stderr)
 report={'cut_id':a.cut,'command':command,'moves':moves,'returncode':p.returncode,'seconds':time.time()-start,'scope':route['scope']}
 transcript.write_text(json.dumps(report,indent=2)+'\n');assert p.returncode==0,p.stderr
 produced=safe(work,route['produced_record']);d=json.loads(produced.read_text());obj=list(map(F,expected['objective_exact']))
 if 'results'in d:
  choices=[v for v in d['results']if objective(v)==obj];assert len(choices)==1;d=choices[0]
 assert objective(d)==obj
 assert support(d)==F(expected['exact_certified_support']),(support(d),expected['exact_certified_support'])
 checked=[]
 for peer in routes:
  pr=peer['original_assembly']
  if not pr or pr['script']!=route['script']or pr['arguments']!=route['arguments']:continue
  e=next(x for x in registry['records']if x['cut_id']==peer['cut_id']);pd=json.loads(safe(work,pr['produced_record']).read_text());po=list(map(F,e['objective_exact']))
  if 'results'in pd:
   vals=[v for v in pd['results']if objective(v)==po];assert len(vals)==1;pd=vals[0]
  assert objective(pd)==po and support(pd)==F(e['exact_certified_support']),peer['cut_id'];checked.append(peer['cut_id'])
 report.update(status='PASS_ORIGINAL_ASSEMBLY_EXACT_SELECTED_SUPPORT',produced_record=route['produced_record'],produced_sha256=sha(produced),objective_exact=expected['objective_exact'],support_exact=expected['exact_certified_support'],selected_cuts_checked_for_this_shared_invocation=checked)
 transcript.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
