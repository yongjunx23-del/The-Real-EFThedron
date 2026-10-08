"""Direct head enclosures for an exact positive combination of complete witnesses.

Only the listed spins are reintegrated. Other spins retain the convex sum
of complete base costs. Equal source/threshold gap intervals are required.
"""
from pathlib import Path
import sys,bisect
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/boundary_length_20261007/head'))
import certify_head as driver
from finite_fast import *
from replay_left_family_fast import CellCache
from head_certification_spin0 import support_upper
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dy(d):return F(int(d[0]))*F(2)**int(d[1])
def upq(v):return dy(v.upper().man_exp())
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--spin',type=int,choices=[0,2],required=True)
    ap.add_argument('--bits',type=int,default=256);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    assert not a.out.exists();ctx.prec=a.bits;tic=time.time();ell=a.spin;deps={}
    def read(p):
        p=ROOT/p;deps[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    cp=[Path('output/boundary_length_20261007/complete_angle85_banded_certificate.json'),Path('output/competitive_20261007/complete_g2_certificate.json')]
    comp=[read(p)for p in cp]
    for c in comp:
        for p,h in c['sha256_manifest'].items():assert sha(ROOT/p)==h;deps[p]=h
    wp=[ROOT/'output/full_direction_20261007/right/angle85_220.npz',ROOT/'output/competitive_20261007/robust/g2_robust_stage1.npz']
    bases=[np.load(p)for p in wp]
    target=ROOT/'output/weak_regions_20261007/upper/angle86_refined.npz';tar=np.load(target)
    direction=[F(float(v))for v in tar['direction']]
    bd=[[F(float(v))for v in w['direction']]for w in bases]
    weights=[direction[0]/bd[0][0]];weights.append(direction[1]-weights[0]*bd[0][1]);assert all(x>0 for x in weights)
    assert bd[1]==[F(0),F(1)] and all(sum(weights[i]*bd[i][j]for i in [0,1])==direction[j]for j in [0,1])
    assert np.array_equal(bases[0]['names'],bases[1]['names'])
    multipliers=[sum(weights[i]*F(float(bases[i]['multipliers'][j]))for i in [0,1])for j in range(len(bases[0]['names']))]
    paths=[Path(f'output/boundary_length_20261007/head/angle85/spin{ell}/head_l{ell}.json'),Path(f'output/all_direction_certification_20261007/head/g2_robust_stage1/{"low0_refined"if ell==0 else "low"}/head_l{ell}.json')]
    records=[read(p)for p in paths]
    moments=[ROOT/'output/boundary_length_20261007/head/angle85/moments128',ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1/moments128']
    for i in [0,1]:
        assert sha(wp[i])==comp[i]['witness_sha256']==records[i]['witness_sha256']
        assert sha(moments[i]/'manifest.json')==records[i]['source_manifest_sha256']
        deps[str((moments[i]/'manifest.json').relative_to(ROOT))]=sha(moments[i]/'manifest.json')
    models=[FiniteCentered(wp[i],ell,moments[i],N=128)for i in [0,1]];caches=[CellCache(m)for m in models]
    assert records[0]['gap_beta_intervals']==records[1]['gap_beta_intervals']
    gaps=[tuple(map(F,x))for x in records[0]['gap_beta_intervals']]
    # The recorded Arb ball strings enclose the original gap bounds; parse
    # each with outward rounding and take its upper endpoint again.
    gap_components=[dy(r['gap_upper_dyadic'])if 'gap_upper_dyadic'in r else upq(arb(r['gap_upper']))for r in records]
    gap=sum(w*g for w,g in zip(weights,gap_components));total=gap
    oldcells=[[ (F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in r['cells_data']]for r in records]
    starts=[[c[0]for c in cells]for cells in oldcells]
    knots=sorted(set([x for cells in oldcells for lo,hi,up in cells for x in [lo,hi]]+[x for g in gaps for x in g]))
    cells=[(lo,hi)for lo,hi in zip(knots,knots[1:])if not any(lo>=ga and hi<=gb for ga,gb in gaps)]
    driver.exact_partition(cells,gaps);data=[];W=list(map(exact,weights));fac=models[0].factor
    for index,(lo,hi)in enumerate(cells):
        width=exact(hi-lo);fallback=F(0)
        for i in [0,1]:
            j=bisect.bisect_right(starts[i],lo)-1;ol,oh,ou=oldcells[i][j];assert ol<=lo<hi<=oh
            fallback+=weights[i]*ou;caches[i].clear()
        fields=[c.fields(lo,hi)for c in caches];U=sum((W[i]*fields[i][0]for i in [0,1]),arb(0));V=sum((W[i]*fields[i][1]for i in [0,1]),arb(0))
        options=[(fallback,'convex_base_cells')];ri=fac*width*support_upper(U,V)
        if ri.is_finite():options.append((upq(ri),'mixture_riemann'))
        jets=[c.series(lo,hi)for c in caches]
        uj=sum((W[i]*jets[i][0]for i in [0,1]),arb_series([]));vj=sum((W[i]*jets[i][1]for i in [0,1]),arb_series([]));err=sum((W[i]*jets[i][2]for i in [0,1]),arb(0))
        es=arb(0)
        for x in [lo,hi]:
            ends=[c.fields(x,x)for c in caches];eu=sum((W[i]*ends[i][0]for i in [0,1]),arb(0));ev=sum((W[i]*ends[i][1]for i in [0,1]),arb(0));es+=support_upper(eu,ev)
        rem=width**3/8*(2*uj[2].abs_upper()+4*vj[2].abs_upper())+4*width*err
        tr=fac*(width*es/2+rem)
        if tr.is_finite():options.append((upq(tr),'mixture_convex_endpoints'))
        if hi<F(1,2)or(ell and F(1,2)<=hi<F(99,100)):
            nn=[c.norm(lo,hi)for c in caches];nu=sum((W[i]*nn[i][0]for i in [0,1]),arb(0));nv=sum((W[i]*nn[i][1]for i in [0,1]),arb(0))
            if nv.upper()<0:
                sm=fac*width*nu.abs_upper()**2*nn[0][2].abs_upper()*nn[0][3].abs_upper()**2/(2*(-nv.upper()))
                if sm.is_finite():options.append((upq(sm),'mixture_negative_V'))
        value,method=min(options);assert value>=0;total+=value;data.append(dict(a=str(lo),b=str(hi),upper_exact=str(value),method=method))
        if(index+1)%1000==0:print('MIXTURE',ell,a.bits,index+1,len(cells),time.time()-tic,flush=True)
    deps.update(driver.code_hashes());deps[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    deps['research/all_direction_certification_20261007/head/replay_left_family_fast.py']=sha(ROOT/'research/all_direction_certification_20261007/head/replay_left_family_fast.py')
    r=dict(status='PASS_EXACT_POSITIVE_MIXTURE_HEAD_COMPONENT',spin=ell,precision_bits=a.bits,objective_exact=list(map(str,direction)),weights_exact=list(map(str,weights)),base_certificate_paths=list(map(str,cp)),base_witness_sha256=[sha(p)for p in wp],names=bases[0]['names'].tolist(),multipliers_exact=list(map(str,multipliers)),gap_intervals=records[0]['gap_beta_intervals'],gap_components_exact=list(map(str,gap_components)),gap_upper_exact=str(gap),upper_exact=str(total),upper_float=float(total),cells=data,source_record_paths=list(map(str,paths)),dependencies=deps,elapsed_seconds=time.time()-tic,scope='Only this spin, allenergies; other spins require positivehomogeneity/subadditivity assembly from completebasecertificates.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2));print('COMPLETE',ell,a.bits,float(total),time.time()-tic,flush=True)
if __name__=='__main__':main()
