"""Select the last figure's proof dependency closure; never modify source records.

JSON bindings and file references are followed recursively. Python imports and
literal file inputs are added. Original records stay byte-identical in records/.
The separate checker verifies integrity and declared closure, not the mathematics.
"""
from pathlib import Path, PureWindowsPath
import argparse, ast, collections, functools, hashlib, importlib.metadata, json, platform, re, shutil, sys, zipfile

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/FIGURE_REPRO_20261008'
HEX=re.compile(r'^[0-9a-fA-F]{64}$')
SUFFIXES={'.json','.jsonl','.log','.npz','.npy','.csv','.txt','.md','.py','.gz','.png','.pdf','.svg','.toml'}
EXCLUDE_PARTS={'__pycache__','.git','FIGURE_REPRO_20261008','figure_repro_20261008'}
ACTIVE=[0,1,2,3,4,5,6,8,15,16,17,20,21,28,29,31,33,34,36,37,39,41,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59]

@functools.lru_cache(None)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
 return h.hexdigest()

def normal(value):
 v=str(value).replace('\\','/')
 prefix=ROOT.as_posix().rstrip('/')+'/'
 if v.lower().startswith(prefix.lower()):v=v[len(prefix):]
 if PureWindowsPath(v).is_absolute() or v.startswith('/') or '..'in Path(v).parts:return None
 return v

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seal',action='store_true');ap.add_argument('--extend',action='store_true');a=ap.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 if a.seal:
  seal();return
 files={};queue=collections.deque();bindings=[];refs=[];missing=[];imports=[]
 processed=set()
 if a.extend:
  prior=json.loads((OUT/'DEPENDENCIES.json').read_text())
  files={r['path']:{'parents':set(r['parents']),'reasons':set(r['reasons'])}for r in prior['files']}
  processed=set(files);bindings=prior['record_bindings'];refs=prior['record_references'];missing=prior['missing_bound_references'];imports=prior['local_imports']
  capture='output/competitive_20261007/robust/g2_robust_stage1_capture.json'
  bindings=[b for b in bindings if not(b['parent']==capture and b['field']=='/sha256')]
 module_index=collections.defaultdict(list)
 for p in (ROOT/'research').rglob('*.py'):
  if not EXCLUDE_PARTS.intersection(p.parts):module_index[p.stem].append(p)
 def add(p,parent,reason):
  p=Path(p)
  if not p.is_file() or not p.is_relative_to(ROOT) or (EXCLUDE_PARTS.intersection(p.relative_to(ROOT).parts)and reason!='operational-appendix-source') or p.suffix in {'.zip','.pyc'}:return
  rel=p.relative_to(ROOT).as_posix()
  if rel not in files:files[rel]={'parents':set(),'reasons':set()};queue.append(rel)
  files[rel]['parents'].add(parent);files[rel]['reasons'].add(reason)
 def resolve(v,parent,required=False):
  if not isinstance(v,str)or len(v)>1000:return None
  n=normal(v)
  if n is None:return None
  candidates=[ROOT/n,(ROOT/parent).parent/n,ROOT/'research'/n]
  for p in candidates:
   if p.is_file():return p
  if required:missing.append({'parent':parent,'path':v,'kind':'bound-file-reference'})
  return None
 def ref(v,parent,reason):
  p=resolve(v,parent)
  if p is not None:
   add(p,parent,reason);refs.append({'parent':parent,'path':p.relative_to(ROOT).as_posix(),'kind':reason});return p
  return None
 def walk(o,parent,context=''):
  if isinstance(o,dict):
   for k,v in o.items():
    if isinstance(v,str)and HEX.fullmatch(v) and isinstance(k,str):
     if '/'in k or '\\'in k or Path(k).suffix in SUFFIXES:
      p=resolve(k,parent,True)
      if p is not None:
       add(p,parent,'record-sha-binding');bindings.append({'parent':parent,'path':p.relative_to(ROOT).as_posix(),'sha256':v.lower(),'field':context+'/'+k})
    if isinstance(v,str) and HEX.fullmatch(v):
     bases=[]
     if k=='sha256':bases=['path','file','source','certificate','witness']
     elif k.endswith('_sha256'):bases=[k[:-7],k[:-7]+'_path',k[:-7]+'_file']
     for b in bases:
      if isinstance(o.get(b),str):
       p=resolve(o[b],parent)
       if p is not None:
        add(p,parent,'record-sha-binding');bindings.append({'parent':parent,'path':p.relative_to(ROOT).as_posix(),'sha256':v.lower(),'field':context+'/'+k});break
    walk(v,parent,context+'/'+str(k))
  elif isinstance(o,list):
   for i,v in enumerate(o):walk(v,parent,context+'/'+str(i))
  elif isinstance(o,str):
   if '.'not in o or len(o)>1000:return
   v=o.replace('\\','/')
   if Path(v).suffix.lower()in SUFFIXES:ref(o,parent,'record-file-reference')
 def source(rel):
  p=ROOT/rel
  try:tree=ast.parse(p.read_text(encoding='utf-8-sig'))
  except (SyntaxError,UnicodeError):return
  for node in ast.walk(tree):
   if isinstance(node,(ast.Import,ast.ImportFrom)):
    names=[x.name for x in node.names]if isinstance(node,ast.Import)else [node.module or '']
    for name in names:
     stem=name.split('.')[-1];choices=module_index.get(stem,[])
     if choices:
      same=[q for q in choices if q.parent==p.parent]
      nearby=[q for q in choices if q.parent in p.parents]
      chosen=same or nearby or choices
      for q in chosen:add(q,rel,'local-python-import')
      imports.append({'source':rel,'module':name,'resolved':[q.relative_to(ROOT).as_posix()for q in chosen]})
   elif isinstance(node,ast.Constant)and isinstance(node.value,str):
    v=node.value
    if len(v)<1000 and '\n'not in v and Path(v.replace('\\','/')).suffix.lower()in SUFFIXES:ref(v,rel,'source-literal-existing-file')
 seeds=json.loads((ROOT/'output/six_targets_20261007/comparison/certified_polygon.json').read_text())
 for i in ACTIVE:
  if i==28:continue
  ref(seeds['cuts'][i]['source'],'selection','selected-cut-certificate')
 if (OUT/'data/proof_seeds.json').is_file():
  extra=json.loads((OUT/'data/proof_seeds.json').read_text())
  assert extra['selected_cut_indices']==ACTIVE
  for name in extra['paths']:ref(name,'selection','inventory-explicit-seed')
 opt=ROOT/'research/figure_repro_20261008/optimization_seed_manifest.json'
 if opt.is_file():
  extra=json.loads(opt.read_text())
  rows=extra['files'] if 'files'in extra else extra['seeds']
  for row in rows:ref(row['path'],'selection','optimization-explicit-seed')
 # Preserve human proof-method reviews and recorded invocations adjacent to
 # included modules/witnesses. These are inputs to reproduction, not new proofs.
 for rel in list(files):
  p=ROOT/rel
  if p.suffix=='.py':
   for pattern in ['*.md','*.txt']:
    for q in p.parent.glob(pattern):add(q,rel,'adjacent-method-documentation')
  if p.suffix=='.npz' or 'complete_'in p.name:
   for pattern in ['*command*.json','*command*.jsonl','*campaign*.json']:
    for q in p.parent.glob(pattern):add(q,rel,'recorded-producer-invocation')
 for name in ['requirements.txt','requirements-handoff.txt','requirements-certification.txt']:
  add(ROOT/name,'selection','runtime-requirements')
 # Reproduction code for the figure and its precise inputs are part of the package.
 for p in (ROOT/'research/primal_dual_no_dagger_20261008').glob('*'):add(p,'selection','figure-source')
 for p in (ROOT/'output/primal_dual_no_dagger_20261008').glob('*'):add(p,'selection','figure-artifact')
 finalplot=json.loads((ROOT/'output/primal_dual_no_dagger_20261008/plot_data.json').read_text())
 bindings.append({'parent':'output/primal_dual_no_dagger_20261008/plot_data.json','path':'research/primal_dual_no_dagger_20261008/make_figure.py','sha256':finalplot['script_sha256'],'field':'/script_sha256 (matching named figure renderer)'})
 # Omit previous comparison registries' broad provenance: the selected-cut registry
 # authenticates the original full polygon separately, without adopting unused cuts.
 add(ROOT/'output/six_targets_20261007/comparison/certified_polygon.json','selection','figure-exact-polygon')
 while queue:
  rel=queue.popleft()
  if rel in processed:continue
  processed.add(rel);p=ROOT/rel
  if p.suffix=='.json':
   if rel=='output/six_targets_20261007/comparison/certified_polygon.json':continue
   try:walk(json.loads(p.read_text(encoding='utf-8-sig')),rel)
   except (UnicodeError,json.JSONDecodeError):pass
  elif p.suffix=='.py':source(rel)
  # Original NPZ witnesses store settings and accompanying optimization metadata
  # in same-stem JSON records; cache files can be named therein.
  elif p.suffix=='.npz':
   if p.with_suffix('.json').is_file():add(p.with_suffix('.json'),rel,'witness-numerical-settings')
  if len(processed)%250==0:print('CLOSURE',len(processed),len(files),flush=True)
 # Small local binary runtime used in the original proof; source inputs remain
 # portable if python-flint is installed instead on a non-Windows host.
 for p in (ROOT/'tmp/certification_deps').rglob('*'):
  if p.is_file() and '__pycache__'not in p.parts and p.suffix!='.pyc':add(p,'selection','original-python-flint-runtime')
 for p in (ROOT/'research/figure_repro_20261008').iterdir():
  if p.is_file()and p.suffix in ['.py','.tex','.md','.json']:add(p,'selection','operational-appendix-source')
 runtime={'python':sys.version,'platform':platform.platform(),'packages':{n:importlib.metadata.version(n)for n in ['numpy','scipy','sympy','cvxpy','clarabel','matplotlib','mpmath']},'scope':'Version metadata only; no credentials or environment-variable dump.'}
 (OUT/'environment.json').write_text(json.dumps(runtime,indent=2)+'\n')
 mismatches=[]
 capture='output/competitive_20261007/robust/g2_robust_stage1_capture.json'
 frozen='output/competitive_20261007/robust/g2_robust_stage1.npz'
 if capture in files:
  d=json.loads((ROOT/capture).read_text());assert sha(ROOT/frozen)==d['sha256']
  bindings=[b for b in bindings if not(b['parent']==capture and b['field']=='/sha256')]
  bindings.append({'parent':capture,'path':frozen,'sha256':d['sha256'],'field':'/sha256 (frozen snapshot; source checkpoint is historical mutable provenance)'})
 for b in bindings:
  if sha(ROOT/b['path'])!=b['sha256']:mismatches.append(b)
 rows=[];total=0
 for rel,meta in sorted(files.items()):
  p=ROOT/rel;dest=OUT/'records'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
  if not dest.exists()or sha(dest)!=sha(p):shutil.copy2(p,dest)
  n=p.stat().st_size;total+=n
  rows.append({'path':rel,'bytes':n,'sha256':sha(p),'reasons':sorted(meta['reasons']),'parents':sorted(meta['parents'])})
 result={'format':'FIGURE_SCOPED_DEPENDENCY_CLOSURE_V1','original_workspace':str(ROOT),'selected_cut_indices':ACTIVE,'files':rows,'file_count':len(rows),'uncompressed_bytes':total,'record_bindings':bindings,'record_references':refs,'local_imports':imports,'missing_bound_references':missing,'original_record_hash_mismatches':mismatches,'binding_resolution_notes':[{'record':capture,'field':'sha256','meaning':'Hash authenticates the captured fixed witness g2_robust_stage1.npz. The source field names a running checkpoint that advanced after capture; its later bytes are numerical provenance, not the selected proof witness.'}],'scope':'File integrity and explicit recursively discovered record/source dependency closure. Does not establish mathematical correctness; no fresh interval arithmetic replay is implied. Original aggregate records may contain dependencies for additional family members.'}
 (OUT/'DEPENDENCIES.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'files':len(rows),'bytes':total,'bindings':len(bindings),'missing':len(missing),'mismatches':len(mismatches)},indent=2))

def seal():
 archive=OUT.with_suffix('.zip');assert not archive.exists(),'Immutable package already exists; use a new version.'
 rows=[]
 for p in sorted(OUT.rglob('*')):
  if p.is_file()and p.name!='MANIFEST.json'and '__pycache__'not in p.parts:
   rows.append({'path':p.relative_to(OUT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
 manifest={'format':'FIGURE_REPRO_PACKAGE_V1','files':rows,'file_count':len(rows),'scope':'Integrity manifest. Mathematical checking is described separately in the appendix and verification guide.'}
 (OUT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
 with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True)as z:
  for p in sorted(OUT.rglob('*')):
   if p.is_file()and '__pycache__'not in p.parts:z.write(p,p.relative_to(OUT).as_posix())
 with zipfile.ZipFile(archive)as z:
  for r in rows:
   h=hashlib.sha256();n=0
   with z.open(r['path'])as f:
    for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block);n+=len(block)
   assert n==r['bytes']and h.hexdigest()==r['sha256']
 digest=sha(archive);archive.with_suffix('.zip.sha256').write_text(digest+'  '+archive.name+'\n')
 print(json.dumps({'archive':str(archive),'sha256':digest,'bytes':archive.stat().st_size,'files':len(rows),'readback':'PASS'},indent=2))

if __name__=='__main__':main()
