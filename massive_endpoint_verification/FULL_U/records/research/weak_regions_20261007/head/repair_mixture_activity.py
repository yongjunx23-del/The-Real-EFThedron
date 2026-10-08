"""Replace normalized candidates outside the Centered real-source validity window.
Only nine selected spin2 cells need repair; all other saved mixture cells remain
unchanged. Actual interval fields and convex-endpoint formula respect source
activity. Raw records are kept intact as diagnostics, not full certificates.
"""
from mixture_head import *
import copy

def repair(bits):
    ctx.prec=bits;tic=time.time();folder=ROOT/'output/weak_regions_20261007/mixture86';path=folder/f'head_l2_{bits}.json';r=json.loads(path.read_text());assert r['spin']==2 and r['precision_bits']==bits
    for p,h in r['dependencies'].items():assert sha(ROOT/p)==h,p
    wp=[ROOT/'output/full_direction_20261007/right/angle85_220.npz',ROOT/'output/competitive_20261007/robust/g2_robust_stage1.npz'];moments=[ROOT/'output/boundary_length_20261007/head/angle85/moments128',ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1/moments128'];models=[FiniteCentered(wp[i],2,moments[i],N=128)for i in[0,1]];caches=[CellCache(m)for m in models];weights=list(map(F,r['weights_exact']));W=list(map(exact,weights));fac=models[0].factor
    sources=[json.loads((ROOT/p).read_text())for p in r['source_record_paths']];oldcells=[[(F(c['a']),F(c['b']),dy(c['upper_dyadic']))for c in d['cells_data']]for d in sources];starts=[[c[0]for c in rows]for rows in oldcells]
    result=copy.deepcopy(r);changes=[]
    for index,c in enumerate(r['cells']):
        lo,hi=F(c['a']),F(c['b'])
        if c['method']!='mixture_negative_V' or hi*hi<=F(79,119):continue
        width=exact(hi-lo);fallback=F(0)
        for i in[0,1]:
            j=bisect.bisect_right(starts[i],lo)-1;ol,oh,ou=oldcells[i][j];assert ol<=lo<hi<=oh;fallback+=weights[i]*ou;caches[i].clear()
        fields=[cc.fields(lo,hi)for cc in caches];U=sum((W[i]*fields[i][0]for i in[0,1]),arb(0));V=sum((W[i]*fields[i][1]for i in[0,1]),arb(0));options=[(fallback,'convex_base_cells')];ri=fac*width*support_upper(U,V)
        if ri.is_finite():options.append((upq(ri),'mixture_riemann'))
        jets=[cc.series(lo,hi)for cc in caches];uj=sum((W[i]*jets[i][0]for i in[0,1]),arb_series([]));vj=sum((W[i]*jets[i][1]for i in[0,1]),arb_series([]));err=sum((W[i]*jets[i][2]for i in[0,1]),arb(0));es=arb(0)
        for x in[lo,hi]:
            ends=[cc.fields(x,x)for cc in caches];eu=sum((W[i]*ends[i][0]for i in[0,1]),arb(0));ev=sum((W[i]*ends[i][1]for i in[0,1]),arb(0));es+=support_upper(eu,ev)
        rem=width**3/8*(2*uj[2].abs_upper()+4*vj[2].abs_upper())+4*width*err;tr=fac*(width*es/2+rem)
        if tr.is_finite():options.append((upq(tr),'mixture_convex_endpoints'))
        value,method=min(options);assert value>=0;new=dict(a=c['a'],b=c['b'],upper_exact=str(value),method=method);result['cells'][index]=new;changes.append(dict(index=index,old=copy.deepcopy(c),replacement=new))
    assert len(changes)==9
    assert all(c['method']!='mixture_negative_V'or F(c['b'])**2<=F(79,119)for c in result['cells'])
    total=F(result['gap_upper_exact'])+sum((F(c['upper_exact'])for c in result['cells']),F(0));result.update(upper_exact=str(total),upper_float=float(total),source_activity_repair=True,raw_record=str(path.relative_to(ROOT)),raw_record_sha256=sha(path),changed_cells=changes,repair_seconds=time.time()-tic,raw_record_not_a_complete_certificate=True)
    result['dependencies'].update({str(path.relative_to(ROOT)):sha(path),str(Path(__file__).relative_to(ROOT)):sha(__file__)})
    out=folder/f'head_l2_{bits}_source_activity_repaired.json';assert not out.exists();out.write_text(json.dumps(result,indent=2));print('REPAIRED',bits,'cells',len(changes),'upper',float(total),'delta',float(total-F(r['upper_exact'])),'seconds',time.time()-tic,flush=True)
if __name__=='__main__':
    for bits in[256,384]:repair(bits)
