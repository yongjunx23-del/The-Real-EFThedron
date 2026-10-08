"""Authenticate literal relocations and expose every extra mathematical edit."""
from pathlib import Path
import json,hashlib,difflib
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/six_targets_20261007/source_review'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return p.read_text(encoding='utf-8')
records=[]
for rel,kind in [('research/six_targets_20261007/upper95_direct/adapter_manifest.json','head_compact'),('research/six_targets_20261007/tail_upper95/initial_adaptation_manifest.json','tail')]:
 d=json.loads(read(ROOT/rel))
 for r in d.get('files',d.get('records')):
  old,new=ROOT/r['source'],ROOT/r['output'];assert sha(old)==r['source_sha256']
  expected=read(old)
  for a,b in d['replacements']:expected=expected.replace(a,b)
  for a,b in r.get('additional_replacements',[]):
   assert a in expected;expected=expected.replace(a,b)
  diff=''.join(difflib.unified_diff(expected.splitlines(True),read(new).splitlines(True),fromfile='declared transformation',tofile=str(new.relative_to(ROOT))))
  if kind=='head_compact':assert not diff
  if diff:print(diff)
  records.append(dict(path=str(new.relative_to(ROOT)),sha256=sha(new),matches_declared_transformation=not diff,diff=diff))
(OUT/'upper95_adapter_diffs.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print('CHECKED',len(records),'files')
d=json.loads(read(ROOT/'research/six_targets_20261007/upper95_direct/assembly_adapter_manifest.json'))
assembly=[]
for r in d['files']:
 old,new=ROOT/r['source'],ROOT/r['result'];assert sha(old)==r['source_sha256']
 expected=read(old)
 for a,b in d['replacements']:expected=expected.replace(a,b)
 assert expected==read(new)and sha(new)==r['result_sha256']
 assembly.append(dict(path=str(new.relative_to(ROOT)),sha256=sha(new),literal_relocation_equal=True))
(OUT/'upper95_assembly_diffs.json').write_text(json.dumps(assembly,indent=2),encoding='utf-8')
print('CHECKED',len(assembly),'literal assemblies')
