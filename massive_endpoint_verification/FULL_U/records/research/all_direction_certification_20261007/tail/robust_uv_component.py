"""Scoped finite-spin/energy sign proof, with arbitrary objective b included."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,sys,time
ROOT=Path(__file__).resolve().parents[3];sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/lower_g0_upper_g2_20261006'),str(ROOT/'research')]
from tail_finite_uv_region import ctx,ar
from contracted_forward import ContractedRegion as Region
OUT=ROOT/'output/all_direction_certification_20261007/tail'
WITNESSES={'g2_robust_stage1':'output/competitive_20261007/robust/g2_robust_stage1.npz'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(x):
    assert x.is_finite();m,e=x.upper().man_exp();return F(int(m))*F(2)**int(e)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--label',choices=WITNESSES,required=True);ap.add_argument('--lmin',type=int,default=1000)
    ap.add_argument('--hi',default='256');ap.add_argument('--kappa',default='1/1000');ap.add_argument('--precision',type=int,default=280);ap.add_argument('--screen',action='store_true');a=ap.parse_args()
    tic=time.time();ctx.prec=a.precision;kappa=F(a.kappa);assert kappa>=0;folder=OUT/a.label;prep=json.loads((folder/'prepared_components.json').read_text())
    for p,h in prep['sha256_manifest'].items():assert sha(ROOT/p)==h,p
    obj=Region(ROOT/WITNESSES[a.label],folder/'joint_moments.pkl',folder/'source_norms.json',a.lmin)
    b=F(prep['objective_exact'][1]);field=lambda lo,hi:obj.field(lo,hi)+ar(b)
    if a.screen:
        rows=[]
        for t in map(F,['0','1/2','1','3/2','2','9/4','5/2','3','4','8','16','32','64','128','256']):
            if t and 4+F(a.lmin*(a.lmin+1))/t<=48:continue
            v=field(t,t);rows.append(dict(tau=str(t),value=str(v),upper_exact=str(dyad(v))))
        result=dict(status='SCREEN_ONLY',label=a.label,lmin=a.lmin,points=rows,seconds=time.time()-tic)
    else:
        end=F(a.hi);stack=[(F(0),end,0)];leaves=[];failed=[];count=0
        while stack:
            lo,hi,depth=stack.pop();v=field(lo,hi);count+=1
            if v.is_finite() and v.upper()<-ar(kappa):leaves.append(dict(lo=str(lo),hi=str(hi),upper_exact=str(dyad(v))));continue
            mid=(lo+hi)/2
            if depth>=24 or (depth>12 and field(mid,mid).lower()>-ar(kappa)):
                failed.append(dict(lo=str(lo),hi=str(hi),value=str(v)));break
            stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
            if count%2000==0:print(json.dumps(dict(evaluations=count,leaves=len(leaves),seconds=time.time()-tic)),flush=True)
        leaves.sort(key=lambda r:F(r['lo']));passed=not stack and not failed
        if passed:
            assert leaves[0]['lo']=='0' and F(leaves[-1]['hi'])==end
            assert all(F(x['hi'])==F(y['lo']) for x,y in zip(leaves,leaves[1:]))
        result=dict(status='PASS' if passed else 'FAILED_ENVELOPE',label=a.label,lmin=a.lmin,tau_domain=['0',str(end)],leaves=leaves,failed=failed,evaluations=count,seconds=time.time()-tic)
    result.update(kappa_exact=str(kappa),evaluation='Exact real Chebyshev/rational contraction before interval evaluation',witness_sha256=obj.sha,objective_exact=prep['objective_exact'],precision_bits=a.precision,
        scope='All even l>=lmin, with tau=l(l+1)/(mu-4) in the stated domain. Other tau, lower spins and real-kernel support remain uncovered.',
        sha256_manifest={**prep['sha256_manifest'],str(Path(__file__).relative_to(ROOT)):sha(__file__),str((Path(__file__).parent/'contracted_forward.py').relative_to(ROOT)):sha(Path(__file__).parent/'contracted_forward.py'),str((ROOT/'research/lower_g0_upper_g2_20261006/tail_finite_uv_region.py').relative_to(ROOT)):sha(ROOT/'research/lower_g0_upper_g2_20261006/tail_finite_uv_region.py')})
    dest=folder/f'contracted_finite_uv_L{a.lmin}_{"screen" if a.screen else "proof"}_{a.precision}.json';dest.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['leaves','sha256_manifest']},indent=2))
if __name__=='__main__':main()
