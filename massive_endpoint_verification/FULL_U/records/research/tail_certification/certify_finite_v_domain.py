import sys,time,json,hashlib,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from certified_finite_v import *

def run(kappa=10,maxcells=100000,replay=None):
 ctx.prec=260;bank=FiniteV();cuts=[F(4),F(5),F(8),F(119,10),F(12),F(24),F(48),F(49),F(100),F(101)]
 stack=[('mu',a,b,0) for a,b in zip(cuts,cuts[1:])]+[('w',F(0),F(1,97),0)]
 if replay:
  old=json.loads(Path(replay).read_text());stack=[(c['variable'],F(c['a']),F(c['b']),0) for c in old['leaves']]
  for var,start,end in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
   ranges=sorted((a,b) for k,a,b,d in stack if k==var);assert ranges[0][0]==start and ranges[-1][1]==end
   assert all(a[1]==b[0] for a,b in zip(ranges,ranges[1:]))
 leaves=[];trials=0;worst=arb(-10**20);start=time.time()
 while stack:
  kind,a,b,depth=stack.pop();trials+=1
  if trials>maxcells:raise RuntimeError(f'Cell budget reached {len(leaves)} leaves, {len(stack)}pending')
  x=arb(ar((a+b)/2),ar((b-a)/2));ok=False
  try:
   vals=bank.uv_fields(x) if kind=='w' else bank.fields(x)
   ok=all(v.is_finite() and v.upper()<-kappa for v in vals.values())
   if all(v.is_finite() for v in vals.values()):
    def key(v):
     m,e=v.man_exp();return F(int(m))*F(2)**int(e)
    upper=max((v.upper() for v in vals.values()),key=key)
   else:upper=arb.pos_inf()
  except (AssertionError,ValueError,ZeroDivisionError):upper=arb.pos_inf()
  if ok:
   if upper>worst:worst=upper
   leaves.append({'variable':kind,'a':str(a),'b':str(b),'upper':str(upper),'upper_dyadic':[str(v) for v in upper.man_exp()]})
  else:
   if depth>=35:raise RuntimeError(f'Unresolved {kind} {a}..{b}: {upper}')
   mid=(a+b)/2;stack.extend([(kind,a,mid,depth+1),(kind,mid,b,depth+1)])
  if trials%1000==0:print('trials',trials,'leaves',len(leaves),'remaining',len(stack),'seconds',time.time()-start,flush=True)
 return dict(witness_sha256=hashlib.sha256((ROOT/'output/angles_resummed/angle_000p000.npz').read_bytes()).hexdigest(),spins=bank.spins,claim=f'V_l(mu)/p_l(mu)<-{kappa} for the listed finite spins and all mu>4',worst_upper=str(worst),kappa=kappa,trials=trials,cells=len(leaves),seconds=time.time()-start,leaves=leaves)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--kappa',type=int,default=10);ap.add_argument('--maxcells',type=int,default=100000);ap.add_argument('--replay',type=Path);ap.add_argument('--output',type=Path,default=ROOT/'output/certification_head/V8_40.json');a=ap.parse_args()
 out=run(a.kappa,a.maxcells,a.replay);a.output.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='leaves'},indent=2))
