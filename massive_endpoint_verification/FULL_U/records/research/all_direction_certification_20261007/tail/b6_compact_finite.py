"""Exact b6 finite bridge: all even42..106, 5<=mu<=8."""
from b4_bridge import *
OUT=ROOT/'output/all_direction_certification_20261007/tail/b6_extension'

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=400);ap.add_argument('--replay');a=ap.parse_args();fv.ctx.prec=a.precision;tic=time.time()
 mom=OUT/'bridge_moments';mom.mkdir(exist_ok=True);w=ROOT/'output/angles_resummed/angle_000p000.npz'
 if not(mom/'roy_moments.json').exists():rm.save(rm.build(w,106,64,106),mom/'roy_moments.json')
 if not(mom/'strip_12_24_moments.json').exists():
  for d in sm.build(w,106,64):rm.save(d,mom/f"strip_{d['metadata']['A']}_{d['metadata']['B']}_moments.json")
 original=fv.load;fv.load=lambda path:original(mom/Path(path).name)
 spins=list(range(42,108,2));bank=fv.FiniteV(spins);bank.lx+=fv.ar(F(107,5))
 print('bank ready',time.time()-tic,flush=True)
 saved=json.loads(Path(a.replay).read_text())if a.replay else None
 stack=[(F(c['a']),F(c['b']),0)for c in saved['cells']]if saved else [(F(i),F(i+1),0)for i in range(5,8)]
 cells=[];trials=0
 while stack:
  lo,hi,depth=stack.pop();trials+=1;x=fv.ar((lo+hi)/2)+fv.arb(0,fv.ar((hi-lo)/2))
  try:values={l:v+6 for l,v in bank.fields(x).items()}
  except(AssertionError,ZeroDivisionError,ValueError):values={}
  if len(values)==len(spins)and all(v.is_finite()and v.upper()<-1 for v in values.values()):
   cells.append(dict(a=str(lo),b=str(hi),upper_by_spin={str(l):str(upper(v))for l,v in values.items()}));continue
  assert not saved and depth<25,(lo,hi,{l:str(v)for l,v in values.items()})
  mid=(lo+hi)/2;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  if trials%200==0:print('trials',trials,'cells',len(cells),'seconds',time.time()-tic,flush=True)
 cells.sort(key=lambda c:F(c['a']));assert F(cells[0]['a'])==5 and F(cells[-1]['b'])==8 and all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
 files=[w,Path(__file__),Path(fv.__file__),Path(rm.__file__),Path(sm.__file__),ROOT/'research/validated_quadrature.py']+list(mom.glob('*.json'))+([Path(a.replay).resolve()]if saved else[])
 out=dict(status='PASS_B6_FINITE_COMPACT_BRIDGE',base_witness_sha256=sha(w),objective_exact=['1','6'],delta_X4_exact='107/5',spins=spins,mu_domain=['5','8'],kappa_exact='1',precision_bits=a.precision,cells=cells,seconds=time.time()-tic,trials=trials,reference_certificate_sha256=sha(a.replay)if saved else None,sha256_manifest={str(p.relative_to(ROOT)):sha(p)for p in files})
 (OUT/f'compact_bridge_{a.precision}.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in['cells','sha256_manifest','spins']},indent=2))

if __name__=='__main__':main()
