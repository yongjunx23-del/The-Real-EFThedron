"""Run preserved research code with narrowly scoped path relocation.

Use a disposable COPY of records/ as --workdir for scripts that write results.
Original scientific records and code bytes remain unchanged. Only construction
of pathlib paths is adapted: backslashes become separators and the one original
workspace prefix is mapped to the supplied copy. Other absolute paths are refused
when supplied as the script name; this is not a filesystem security sandbox.
"""
import argparse, json, os, pathlib, runpy, sys

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--package',type=pathlib.Path,default=pathlib.Path(__file__).resolve().parent);ap.add_argument('--workdir',type=pathlib.Path);ap.add_argument('script');ap.add_argument('args',nargs=argparse.REMAINDER);a=ap.parse_args()
 package=a.package.resolve();work=(a.workdir or package/'records').resolve()
 manifest=json.loads((package/'DEPENDENCIES.json').read_text());original=manifest['original_workspace'].replace('\\','/').rstrip('/')
 script=a.script.replace('\\','/');assert not pathlib.PureWindowsPath(script).is_absolute()and '..'not in pathlib.PurePosixPath(script).parts
 target=work/script;assert target.is_file(),target
 old_new=pathlib.Path.__new__;old_init=pathlib.Path.__init__
 destination=str(work).replace('\\','/').rstrip('/')
 def convert(parts):
  values=[]
  for part in parts:
   if isinstance(part,(str,os.PathLike)):
    s=os.fspath(part).replace('\\','/')
    if s.casefold()==destination.casefold()or s.casefold().startswith(destination.casefold()+'/'):pass
    elif s.casefold()==original.casefold():s=str(work)
    elif s.casefold().startswith(original.casefold()+'/'):s=str(work)+s[len(original):]
    part=s
   values.append(part)
  return values
 def relocated_new(cls,*parts,**kwargs):return old_new(cls,*convert(parts),**kwargs)
 def relocated_init(self,*parts,**kwargs):return old_init(self,*convert(parts),**kwargs)
 pathlib.Path.__new__=relocated_new
 pathlib.Path.__init__=relocated_init
 os.chdir(work)
 sys.path[:0]=[str(target.parent),str(work/'research')]
 if sys.platform=='win32':sys.path.insert(0,str(work/'tmp/certification_deps'))
 sys.argv=[str(target)]+a.args
 runpy.run_path(str(target),run_name='__main__')

if __name__=='__main__':main()
