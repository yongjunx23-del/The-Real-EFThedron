"""Direct low-spin interval bounds for the exact83+b6 mixture.

Unequal endpoint/source neighborhoods are covered by their union; a whole
existing nonnegative cell bound is charged whenever that cell meets a new gap.
"""
from pathlib import Path
import sys,time,json,hashlib,bisect,argparse
from upper95_common import *
sys.path[:0]=[str(ROOT/'research/boundary_length_20261007/head'),str(ROOT/'research/full_direction_20261007')]
import certify_head as driver
from finite_fast import *
from head_family_extended import DirectionalZero,DirectionalLow
from head_certification_spin0 import SpinZero,support_upper
from head_certification_lowspin import LowSpin

class Cache:
    def __init__(self,m):self.m=m;self.points={};self.clear()
    def clear(self):self.fields_saved={};self.series_saved={};self.norm_saved={}
    def fields(self,a,b):
        cache=self.points if a==b else self.fields_saved
        if (a,b)not in cache:cache[a,b]=self.m.fields(a,b)
        return cache[a,b]
    def series(self,a,b):
        if (a,b)not in self.series_saved:self.series_saved[a,b]=self.m.fields_series(a,b)
        return self.series_saved[a,b]
    def norm(self,a,b):
        if (a,b)not in self.norm_saved:self.norm_saved[a,b]=self.m.normalized_interval(a,b)
        return self.norm_saved[a,b]

def upq(v):return dy(v.upper().man_exp())
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--spin',type=int,choices=[0,2],required=True)
    ap.add_argument('--bits',type=int,default=256);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    assert not args.out.exists();ell=args.spin;ctx.prec=args.bits;tic=time.time();deps={}
    def read(p):
        p=Path(p);p=p if p.is_absolute()else ROOT/p;deps[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    certificates=[read(p)for p in CERTS]
    for c in certificates:
        assert c['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'
        for p,h in c['sha256_manifest'].items():assert sha(ROOT/p)==h;deps[p]=h
    prep=read(OUT/'preparation.json');assert prep['status']=='PASS_EXACT_LEGACY_SOURCE_CACHE_REBUILD'and prep['all_nested_exact_fields_equal']
    for p,h in prep['sha256_manifest'].items():assert sha(ROOT/p)==h;deps[p]=h
    target=np.load(OUT/'objective.npz');direction=[F(float(v))for v in target['direction']]
    weights=exact_weights(direction,[list(map(F,c['objective_exact']))for c in certificates])
    raw=[read(p)for p in source_paths(ell)]
    for j in[0,1]:assert raw[j]['witness_sha256']==certificates[j]['witness_sha256']==sha(ROOT/WITNESSES[j])
    oldcells,oldgaps,oldbounds=canonical_records(raw,ell)
    cells,gaps,gap_components,gap_ledger=union_gap_data(oldcells,oldgaps,oldbounds)
    moments=ROOT/'output/peak_certification_20261007/upper83/head/moments128'
    assert sha(moments/'manifest.json')==raw[0]['source_manifest_sha256'];deps[str((moments/'manifest.json').relative_to(ROOT))]=sha(moments/'manifest.json')
    first=FiniteCentered(ROOT/WITNESSES[0],ell,moments,N=128)
    legacy_cache=ROOT/'output/certification_head/spin0_source_N64.pkl';before=sha(legacy_cache)
    assert before==prep['source_cache_sha256']
    legacy=SpinZero(ROOT/WITNESSES[1])if ell==0 else LowSpin(ROOT/WITNESSES[1],ell)
    second=(DirectionalZero if ell==0 else DirectionalLow)(legacy,F(6),F(107,5))
    assert sha(legacy_cache)==before
    models=[first,second];caches=[Cache(m)for m in models];W=list(map(exact,weights));fac=4*(2*ell+1)/arb.pi()
    gap=sum((w*g for w,g in zip(weights,gap_components)),F(0));total=gap;data=[]
    starts=[[a for a,b,v in cc]for cc in oldcells]
    for index,(lo,hi)in enumerate(cells):
        width=exact(hi-lo);fallback=F(0)
        for j in[0,1]:
            k=bisect.bisect_right(starts[j],lo)-1;a,b,v=oldcells[j][k];assert a<=lo<hi<=b
            fallback+=weights[j]*v;caches[j].clear()
        fields=[c.fields(lo,hi)for c in caches]
        U=sum((W[j]*fields[j][0]for j in[0,1]),arb(0));V=sum((W[j]*fields[j][1]for j in[0,1]),arb(0))
        options=[(fallback,'convex_base_cells')];ri=fac*width*support_upper(U,V)
        if ri.is_finite():options.append((upq(ri),'mixture_riemann'))
        jets=[c.series(lo,hi)for c in caches]
        uj=sum((W[j]*jets[j][0]for j in[0,1]),arb_series([]));vj=sum((W[j]*jets[j][1]for j in[0,1]),arb_series([]));err=sum((W[j]*jets[j][2]for j in[0,1]),arb(0))
        es=arb(0)
        for x in[lo,hi]:
            ends=[c.fields(x,x)for c in caches];eu=sum((W[j]*ends[j][0]for j in[0,1]),arb(0));ev=sum((W[j]*ends[j][1]for j in[0,1]),arb(0));es+=support_upper(eu,ev)
        rem=width**3/8*(2*uj[2].abs_upper()+4*vj[2].abs_upper())+4*width*err;tr=fac*(width*es/2+rem)
        if tr.is_finite():options.append((upq(tr),'mixture_convex_endpoints'))
        # Both normalized models include the Roy imaginary source here only
        # while mu<=119/10. Spin0 uses only direct fields and series above.
        if ell==2 and hi*hi<=F(79,119):
            nn=[c.norm(lo,hi)for c in caches];nu=sum((W[j]*nn[j][0]for j in[0,1]),arb(0));nv=sum((W[j]*nn[j][1]for j in[0,1]),arb(0))
            if nv.upper()<0:
                sm=fac*width*nu.abs_upper()**2*nn[0][2].abs_upper()*nn[0][3].abs_upper()**2/(2*(-nv.upper()))
                if sm.is_finite():options.append((upq(sm),'mixture_negative_V'))
        value,method=min(options);assert value>=0;total+=value
        data.append(dict(a=str(lo),b=str(hi),upper_exact=str(value),method=method))
        if(index+1)%1000==0:print('MIXTURE95',ell,args.bits,index+1,len(cells),time.time()-tic,flush=True)
    deps.update(driver.code_hashes())
    for p in[Path(__file__),Path(__file__).with_name('upper95_common.py'),ROOT/'research/full_direction_20261007/head_family_extended.py',legacy_cache]:deps[str(p.relative_to(ROOT))]=sha(p)
    definition=[]
    for j,wp in enumerate(WITNESSES):
        st=np.load(ROOT/wp);definition.append(dict(base_witness=wp,base_witness_sha256=sha(ROOT/wp),names=st['names'].tolist(),base_multipliers_exact=[str(F(float(v)))for v in st['multipliers']],weight_exact=str(weights[j]),additional_X4_exact='0'if j==0 else'107/5'))
    result=dict(status='PASS_EXACT_POSITIVE_MIXTURE_HEAD_COMPONENT',spin=ell,precision_bits=args.bits,objective_exact=list(map(str,direction)),weights_exact=list(map(str,weights)),base_certificate_paths=CERTS,base_witness_sha256=[sha(ROOT/p)for p in WITNESSES],coefficient_definition=definition,gap_intervals=[[str(a),str(b)]for a,b in gaps],gap_components_exact=list(map(str,gap_components)),gap_cover_ledger=gap_ledger,gap_upper_exact=str(gap),upper_exact=str(total),upper_float=float(total),cells=data,source_record_paths=source_paths(ell),dependencies=deps,elapsed_seconds=time.time()-tic,scope='Only this spin, every energy. Remaining head and infinite domains use the exact positive combination of complete83 andb6 certificates.')
    args.out.write_text(json.dumps(result,indent=2));print('COMPLETE',ell,args.bits,float(total),time.time()-tic,flush=True)
if __name__=='__main__':main()
