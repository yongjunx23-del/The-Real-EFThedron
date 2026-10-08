"""Exact nonpositive objective shifts on an already certified head partition.

The unchanged base model and its authenticated moments evaluate U,V.  A new
NPZ is used solely to define the exact shifted objective: every coefficient
and g0 component must agree bit-for-bit with the complete base witness.
Each saved base cell and cap is a valid monotonicity fallback.  No sampled
value enters an upper bound; new choices use directed Arb arithmetic.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/boundary_length_20261007/head'))
import certify_head_constant as const
from replay_left_family_fast import CellCache,pscaled
from head_certification_spin0 import support_upper
from finite_fast import *

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fracdy(d):return F(int(d[0]))*F(2)**int(d[1])
def upperfrac(x):return fracdy(x.upper().man_exp())
def check_manifest(d):
    for name,h in d.items():assert sha(ROOT/name)==h,name

def bounds(m,cache,a,b,base_upper,shifts):
    width=exact(b-a);ball=arb(exact((a+b)/2),exact((b-a)/2))
    u,v,_=cache.fields(a,b);p=pscaled(m,ball)
    endpoints=[]
    for x in (a,b):
        eu,ev,_=cache.fields(x,x);endpoints.append((eu,ev,pscaled(m,exact(x))))
    try:
        uj,vj,err=cache.series(a,b);pj=pscaled(m,arb_series([ball,1],prec=3))
        jet_ok=all(x.is_finite() for x in (uj[2],vj[2],pj[2],err))
    except (ValueError,ZeroDivisionError):jet_ok=False
    norm=None
    if b<F(1,2) or (m.ell and F(1,2)<=b<F(99,100)):
        try:norm=cache.norm(a,b)
        except (ValueError,ZeroDivisionError):pass
    result=[]
    for shift in shifts:
        delta=exact(shift);options=[(base_upper,'base_monotonicity')]
        rr=m.factor*width*support_upper(u,v+delta*p)
        if rr.is_finite():options.append((upperfrac(rr),'shifted_riemann'))
        if jet_ok:
            vs=vj+delta*pj
            rem=width**3/8*(2*uj[2].abs_upper()+4*vs[2].abs_upper())+4*width*err
            es=sum((support_upper(eu,ev+delta*ep)for eu,ev,ep in endpoints),arb(0))
            tr=m.factor*(width/2*es+rem)
            if tr.is_finite():options.append((upperfrac(tr),'shifted_convex_endpoints'))
        if norm is not None:
            nu,nv,pinv,e=norm;vn=nv+delta
            if vn.upper()<0:
                sm=m.factor*width*nu.abs_upper()**2*pinv.abs_upper()*e.abs_upper()**2/(2*(-vn.upper()))
                if sm.is_finite():options.append((upperfrac(sm),'shifted_negative_V'))
        result.append(min(options))
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--spin',type=int,required=True)
    ap.add_argument('--witness',type=Path,action='append',required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--bits',type=int,default=256)
    ap.add_argument('--skip-below',default='0');a=ap.parse_args();ctx.prec=a.bits;tic=time.time()
    assert not a.out.exists()
    complete_path=ROOT/'output/boundary_length_20261007/complete_angle268_certificate.json'
    complete=json.loads(complete_path.read_text());check_manifest(complete['sha256_manifest'])
    summary_path=ROOT/'output/boundary_length_20261007/head/angle268/head_summary_final.json'
    summary=json.loads(summary_path.read_text());rec=next(x for x in summary['records']if x['spin']==a.spin)
    source=ROOT/rec['path'];r=json.loads(source.read_text());base=ROOT/complete['witness']
    assert sha(source)==rec['sha256'];assert sha(base)==complete['witness_sha256']
    bd=np.load(base);dirs=[];shifts=[];witnesses=[]
    for wp in a.witness:
        w=np.load(wp);assert np.array_equal(w['names'],bd['names'])
        assert np.array_equal(w['multipliers'],bd['multipliers'])
        direction=[F(float(x))for x in w['direction']];base_direction=list(map(F,complete['objective_exact']))
        assert direction[0]==base_direction[0]
        delta=direction[1]-base_direction[1];assert delta<=0
        dirs.append(list(map(str,direction)));shifts.append(delta)
        witnesses.append(dict(path=str(wp.relative_to(ROOT)if wp.is_absolute()else wp),sha256=sha(wp)))
    moments=ROOT/'output/boundary_length_20261007/head/angle268/moments128'
    assert sha(moments/'manifest.json')==r['source_manifest_sha256']
    klass=const.driver.MODELS[r['model']];m=klass(base,a.spin,moments,N=128);cache=CellCache(m)
    assert [str(x)for x in m.base.direction]==complete['objective_exact']
    gaps=[tuple(map(F,x))for x in r['gap_beta_intervals']]
    cells=[(F(c['a']),F(c['b']))for c in r['cells_data']]
    const.driver.exact_partition(cells,gaps)
    gap=fracdy(r['gap_upper_dyadic']);totals=[gap for _ in shifts];data=[];skips=0
    for i,c in enumerate(r['cells_data']):
        lo,hi=F(c['a']),F(c['b']);bu=fracdy(c['upper_dyadic'])
        if bu<=F(a.skip_below):vals=[(bu,'base_monotonicity')for _ in shifts];skips+=1
        else:
            cache.clear();vals=bounds(m,cache,lo,hi,bu,shifts)
        for j,(up,method)in enumerate(vals):assert up>=0;totals[j]+=up
        data.append(dict(a=str(lo),b=str(hi),upper_exact=[str(v[0])for v in vals],methods=[v[1]for v in vals]))
        if (i+1)%500==0:print('SHIFT_HEAD',a.spin,a.bits,i+1,len(cells),time.time()-tic,flush=True)
    deps=const.driver.code_hashes();deps[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    deps['research/all_direction_certification_20261007/head/replay_left_family_fast.py']=sha(ROOT/'research/all_direction_certification_20261007/head/replay_left_family_fast.py')
    result=dict(status='PASS_EXACT_SHIFT_FINITE_HEAD_COMPONENT',spin=a.spin,precision_bits=a.bits,
        base_complete_certificate=str(complete_path.relative_to(ROOT)),base_complete_sha256=sha(complete_path),
        base_witness_sha256=sha(base),witnesses=witnesses,objective_exact=dirs,g2_shift_exact=list(map(str,shifts)),
        source_record=str(source.relative_to(ROOT)),source_record_sha256=sha(source),
        source_manifest_sha256=r['source_manifest_sha256'],base_gap_upper_exact=str(gap),
        gap_intervals=r['gap_beta_intervals'],partition_covers_beta_zero_to_one_exactly=True,
        upper_exact=list(map(str,totals)),upper_float=list(map(float,totals)),cells=data,
        cells_reused_without_evaluation=skips,skip_below=a.skip_below,elapsed_seconds=time.time()-tic,
        dependencies=deps,scope='Fixed even spin, all energies; exact U unchanged and V shifted by delta*p, delta<=0. Tail separate.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2))
    print('COMPLETE',a.spin,a.bits,list(map(float,totals)),time.time()-tic,flush=True)
if __name__=='__main__':main()
