"""Old150 b4 compact bridge: spins42,44 and every mu in(4,49]."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/tail_certification'),str(ROOT/'research')]
import certified_finite_v as fv
import certified_roy_moments as rm
import certified_strip_moments as sm
OUT=ROOT/'output/all_direction_certification_20261007/tail/old150_extension'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upper(v):m,e=v.upper().man_exp();return F(int(m))*F(2)**int(e)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=320);ap.add_argument('--replay',type=Path);a=ap.parse_args();fv.ctx.prec=a.precision;tic=time.time()
 mom=OUT/'bridge_moments';mom.mkdir(parents=True,exist_ok=True);w=ROOT/'output/angles_resummed/angle_000p000.npz'
 if not(mom/'roy_moments.json').exists():rm.save(rm.build(w,44,64,44),mom/'roy_moments.json')
 if not(mom/'strip_12_24_moments.json').exists():
  for d in sm.build(w,44,64):rm.save(d,mom/f"strip_{d['metadata']['A']}_{d['metadata']['B']}_moments.json")
 # Redirect only the in-process loader, leaving every frozen source unchanged.
 original=fv.load
 fv.load=lambda path:original(mom/Path(path).name)
 bank=fv.FiniteV([42,44]);bank.lx+=fv.ar(F(72,5))
 cuts=list(map(F,['4','5','8','119/10','12','24','48','49']))
 if a.replay:
  saved=json.loads(a.replay.read_text());assert saved['witness_sha256']==sha(w)
  stack=[(F(c['a']),F(c['b']),0)for c in saved['cells']]
 else:stack=[(x,y,0)for x,y in zip(cuts,cuts[1:])]
 cells=[];trials=0
 while stack:
  lo,hi,depth=stack.pop();trials+=1
  x=fv.ar((lo+hi)/2)+fv.arb(0,fv.ar((hi-lo)/2))
  try:values={l:v+4 for l,v in bank.fields(x).items()}
  except (AssertionError,ZeroDivisionError,ValueError):values={}
  if len(values)==2 and all(v.is_finite()and v.upper()<-1 for v in values.values()):
   cells.append(dict(a=str(lo),b=str(hi),upper_by_spin={str(l):str(upper(v))for l,v in values.items()}));continue
  assert depth<30,(lo,hi,values)
  mid=(lo+hi)/2;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  if trials%1000==0:print('trials',trials,'cells',len(cells),'seconds',time.time()-tic,flush=True)
 cells.sort(key=lambda c:F(c['a']));assert F(cells[0]['a'])==4 and F(cells[-1]['b'])==49
 assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
 files=[w,Path(__file__),Path(fv.__file__),Path(rm.__file__),Path(sm.__file__),ROOT/'research/validated_quadrature.py']+list(mom.glob('*.json'))
 out=dict(status='PASS_COMPACT_BRIDGE',witness_sha256=sha(w),b='4',delta_X4='72/5',spins=[42,44],mu_domain=['4','49'],kappa_exact='1',precision_bits=a.precision,
  cells=cells,trials=trials,seconds=time.time()-tic,claim='V_l(mu)/p_l(mu)<=-1 on both listed spins and the complete compact energy interval, including endpoint limits.',
  source_tables='Fresh exact-polynomial Gauss-Legendre moments, jmax44,N64, built with the existing reviewed moment builders.',
  replay_input_sha256=sha(a.replay)if a.replay else None,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in files})
 target=OUT/f'b4_bridge_{a.precision}.json';target.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['cells','sha256_manifest']},indent=2))
if __name__=='__main__':main()
