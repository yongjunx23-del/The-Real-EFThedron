"""Stronger finite-spin coercivity; preserve the original .001 certificates."""
from finite_uv_bridge220 import *

def main():
 ap=argparse.ArgumentParser();ap.add_argument('reference');ap.add_argument('--precision',type=int,default=400);ap.add_argument('--kappa',default='1/100');ap.add_argument('--replay',action='store_true');a=ap.parse_args()
 ctx.prec=a.precision;tic=time.time();path=Path(a.reference).resolve();ref=json.loads(path.read_text());kap=F(a.kappa)
 assert ref['status']=='PASS_FINITE_UV'and ref['witness_sha256']==sha(WITNESS)
 for p,h in ref['sha256_manifest'].items():assert sha(ROOT/p)==h,p
 if a.replay:assert F(ref['kappa_exact'])==kap
 bank=FiniteUV();records=[];failed=[]
 for old in ref['records']:
  l=old['spin'];t=time.time();stack=[(F(c['a']),F(c['b']),0)for c in old['cells']];cells=[];trials=0
  while stack:
   lo,hi,depth=stack.pop();trials+=1;v=bank.field(l,lo,hi)
   if v.is_finite()and v.upper()<-ar(kap):cells.append(dict(a=str(lo),b=str(hi),upper_exact=str(upper(v))));continue
   if a.replay or depth>=22:failed.append(dict(spin=l,a=str(lo),b=str(hi),value=str(v)));break
   mid=(lo+hi)/2;stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
  cells.sort(key=lambda c:F(c['a']));passed=not failed and not stack
  if passed:
   assert F(cells[0]['a'])==0 and F(cells[-1]['b'])==F(old['tau_end_exact'])
   assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
  records.append(dict(spin=l,status='PASS'if passed else'FAILED_ENVELOPE',tau_end_exact=old['tau_end_exact'],cells=cells,evaluations=trials,seconds=time.time()-t))
  print(json.dumps(dict(spin=l,status=records[-1]['status'],cells=len(cells),evaluations=trials,seconds=time.time()-t)),flush=True)
  if failed:break
 manifest=dict(ref['sha256_manifest'])
 for p in[path,Path(__file__)]:manifest[str(p.relative_to(ROOT))]=sha(p)
 out=dict(status='PASS_FINITE_UV'if not failed else'FAILED_ENVELOPE',witness_sha256=sha(WITNESS),spins=ref['spins'],kappa_exact=str(kap),precision_bits=a.precision,records=records,failed=failed,seconds=time.time()-tic,scope=ref['scope'],acceptance_policy=ref['acceptance_policy'],crossing_omitted_kmin=ref['crossing_omitted_kmin'],crossing_tail_Bernstein_coefficients=ref['crossing_tail_Bernstein_coefficients'],reference_certificate_sha256=sha(path),refinement=not a.replay,sha256_manifest=manifest)
 name=f"finite_uv_strong_{min(ref['spins'])}_{max(ref['spins'])}_{a.precision}.json"
 (OUT/name).write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
