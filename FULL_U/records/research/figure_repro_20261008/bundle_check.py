"""Check package bytes, recorded SHA bindings, and the declared dependency closure.

This deliberately does not label a hash check a mathematical proof replay.
Run from an extracted package: python bundle_check.py --package .
"""
from pathlib import Path, PurePosixPath
import argparse, hashlib, json, time

def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()

def safe(s):
 p=PurePosixPath(s)
 if p.is_absolute()or '..'in p.parts or '\\'in s or ':'in s or p.as_posix()!=s:raise ValueError('Unsafe member: '+s)
 return s

def check(package):
 start=time.time();d=json.loads((package/'DEPENDENCIES.json').read_text())
 rows=d['files'];idx={safe(r['path']):r for r in rows};assert len(idx)==len(rows)==d['file_count']
 assert len({s.casefold()for s in idx})==len(idx)
 digests={};total=0
 for s,r in idx.items():
  p=package/'records'/s;assert p.is_file(),s
  assert p.stat().st_size==r['bytes'],s
  h=sha(p);assert h==r['sha256'],s;digests[s]=h;total+=r['bytes']
 assert total==d['uncompressed_bytes']
 assert not d['missing_bound_references'],d['missing_bound_references']
 assert not d['original_record_hash_mismatches'],d['original_record_hash_mismatches']
 for b in d['record_bindings']:
  assert b['parent']in idx and b['path']in idx,b
  assert digests[b['path']]==b['sha256'],b
 for r in d['record_references']:
  assert r['path']in idx,r
  assert r['parent']=='selection'or r['parent']in idx,r
 for m in d['local_imports']:
  assert m['source']in idx,m
  assert all(v in idx for v in m['resolved']),m
 full=package/'MANIFEST.json';final_count=None
 if full.exists():
  m=json.loads(full.read_text());final_count=len(m['files']);assert final_count==m['file_count']
  for r in m['files']:
   p=package/safe(r['path']);assert p.is_file()and p.stat().st_size==r['bytes']and sha(p)==r['sha256'],r['path']
 return {'status':'PASS_PACKAGE_BYTE_INTEGRITY_AND_DECLARED_CLOSURE','record_files':len(idx),'record_bytes':total,'record_hash_bindings':len(d['record_bindings']),'plain_file_references':len(d['record_references']),'local_import_records':len(d['local_imports']),'final_manifest_files':final_count,'seconds':time.time()-start,'scope':'Integrity and declared dependency closure only. This is not fresh interval arithmetic, an independent kernel derivation, or formal proof verification.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--report',type=Path);a=ap.parse_args()
 result=check(a.package.resolve());print(json.dumps(result,indent=2))
 if a.report:a.report.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
