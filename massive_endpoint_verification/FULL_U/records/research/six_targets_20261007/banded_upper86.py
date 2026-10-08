"""Exact merge of saved real-norm cells and same-witness local UV margins.

No new kernel or integration approximation is introduced. The minimum
certified margin across all 204 finite spins bounds each aggregate A cell.
"""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left
import argparse, hashlib, json, time
ROOT=Path(__file__).resolve().parents[2]
ARCHIVED_ROOT='C:/Users/laure/Documents/ChatGPT/Unitarity codex/'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def workspace_path(name):
    """Resolve archived workspace paths; never read a different old checkout."""
    native=Path(name)
    if native.is_absolute():
        try:
            relative=native.resolve().relative_to(ROOT)
            return ROOT/relative
        except ValueError:
            pass
    s=str(name).replace('\\','/')
    if (len(s)>1 and s[1]==':') or s.startswith('/'):
        assert s.startswith(ARCHIVED_ROOT),'Unrecognized archived absolute workspace path: '+s
        s=s[len(ARCHIVED_ROOT):]
    p=(ROOT/s).resolve();p.relative_to(ROOT);return p

def dyadic_ceiling(x,bits=40):
    assert x>0
    e=x.numerator.bit_length()-x.denominator.bit_length()
    if F(2)**e>x:e-=1
    step=F(2)**(e-bits+1);q=x/step;n=-(-q.numerator//q.denominator)
    out=n*step;assert x<=out<x+step
    return out

class Coercivity:
    def __init__(self,angle):
        self.manifest={};self.visited=set();self.angle=angle
        base=ROOT/'output/six_targets_20261007/upper86'
        u=self.read(base/'tail'/'uv_complete.json')
        assert u['status']=='PASS_COMPLETE_DIRECTIONAL_UV_COERCIVITY'
        self.SHA=u['witness_sha256'];self.objective=u['objective_exact']
        assert len(self.objective)==2 and all(F(x)==F(x) for x in self.objective)
        assert u['finite_first_even_spin']==42 and u['finite_last_even_spin']==448 and u['infinite_first_even_spin']==450
        assert u['energy_domain']==['49','infinity']
        c=self.read(base/'compact'/'complete_compact_coercivity.json')
        self.common(c,'PASS_COMPLETE_COMPACT_COERCIVITY')
        assert c['minimum_even_spin']==42 and c['mu_domain_exact']==['4','49'] and c['lower_energy_endpoint_open'] is True
        self.compact=F(c['kappa_exact']);assert self.compact>0
        f=self.read(workspace_path(u['finite_component']));self.common(f,'PASS_COMPLETE_DIRECTIONAL_FINITE_UV')
        assert f['first_even_spin']==42 and f['last_even_spin']==448 and f['mu_domain']==['49','infinity']
        assert f['spins']==list(range(42,450,2))
        large=[]
        for file in f['large_tau_components']:
            d=self.read(workspace_path(file));self.common(d,'PASS_ANGLE'+angle+'_UNBOUNDED_TAU')
            assert d['lmin']<=108 and d['mu_domain']==['49','infinity'] and d['tau_domain']==['256','infinity']
            k=F(d['kappa_exact']);assert k>0 and F(d['V_over_p_upper_exact'])<-k;large.append((d['precision_bits'],k))
        assert sorted(x[0] for x in large)==[320,400]
        self.large=min(x[1] for x in large);assert self.large>=F(f['large_tau_kappa_exact'])>0
        self.cells={};self.references=[]
        for pair in f['paired_files']:
            pp=workspace_path(pair['producer']);rp=workspace_path(pair['replay'])
            assert sha(pp)==pair['producer_sha256'] and sha(rp)==pair['replay_sha256']
            first=self.read(pp);replay=self.read(rp)
            for d in [first,replay]:
                self.common(d,'PASS_ANGLE'+angle+'_FINITE_UV');assert not d['failed']
                assert all(F(x)<=0 for x in d['crossing_tail_coefficients'])
            assert replay['precision_bits']>first['precision_bits']
            assert replay['reference_certificate_sha256']==sha(pp)
            assert first['spins']==replay['spins']
            assert [r['spin'] for r in first['records']]==first['spins']
            assert [r['spin'] for r in replay['records']]==replay['spins']
            claimed=max(F(first['kappa_exact']),F(replay['kappa_exact']))
            for r,s in zip(first['records'],replay['records']):
                ell=r['spin'];assert ell==s['spin'] and ell not in self.cells
                assert r['status']==s['status']=='PASS' and len(r['cells'])==len(s['cells'])
                cells=[];cursor=F(0)
                for a,b in zip(r['cells'],s['cells']):
                    lo,hi=F(a['a']),F(a['b']);assert lo==F(b['a'])==cursor<hi==F(b['b'])
                    upper=max(F(a['upper_exact']),F(b['upper_exact']));assert upper<-claimed
                    cells.append((lo,hi,-upper));cursor=hi
                assert cursor==min(F(256),F(ell*(ell+1),45))==F(r['tau_end_exact'])==F(s['tau_end_exact'])
                self.cells[ell]=cells
            self.references.append(dict(producer=str(pp.relative_to(ROOT)),replay=str(rp.relative_to(ROOT))))
        assert sorted(self.cells)==list(range(42,450,2))
        self.ends={ell:[c[1] for c in cells] for ell,cells in self.cells.items()}

    def common(self,d,status):
        assert d['status']==status and d['witness_sha256']==self.SHA and d['objective_exact']==self.objective

    def read(self,path):
        path=workspace_path(path);key=str(path.relative_to(ROOT));h=sha(path)
        if key in self.manifest:assert self.manifest[key]==h
        self.manifest[key]=h;d=json.loads(path.read_text())
        if key not in self.visited:
            self.visited.add(key)
            for name,expected in d.get('sha256_manifest',{}).items():
                p=workspace_path(name);assert sha(p)==expected,name
                k=str(p.relative_to(ROOT))
                if k in self.manifest:assert self.manifest[k]==expected
                self.manifest[k]=expected
                if p.suffix=='.json':self.read(p)
        return d

    def coverage(self,ell,lo,hi):
        assert 0<=lo<hi<=1;out=[];T=ell*(ell+1);interval=None;large=False
        compact=lo*lo<F(45,49)
        if compact:out.append(self.compact)
        if hi*hi>F(45,49):
            ta=F(T)*(1-hi*hi)/(4*hi*hi)
            tb=min(F(T,45),F(T)*(1-lo*lo)/(4*lo*lo)) if lo else F(T,45)
            if tb>256:
                assert ell>=108;out.append(self.large);large=True
            tb=min(tb,F(256))
            if ta<=tb:
                cs=self.cells[ell];i=bisect_left(self.ends[ell],ta);start=i
                assert i<len(cs) and cs[i][0]<=ta<=cs[i][1]
                while i<len(cs) and cs[i][0]<=tb:out.append(cs[i][2]);i+=1
                assert i>start and cs[i-1][1]>=tb
                interval=[ell,start,i]
        assert out
        return min(out),compact,large,interval

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--angle',choices=['upper86'],required=True);ap.add_argument('--out');args=ap.parse_args();tic=time.time()
    co=Coercivity(args.angle);base=ROOT/'output/six_targets_20261007/upper86/real'
    source=base/'finite_split_replay.json';d=co.read(source)
    assert d['status']=='PASS_REAL_NORM_SPLIT_AND_PRECISION_REPLAY' and d['witness_sha256']==co.SHA
    assert d['mode']=='finite' and d['first_even_spin']==42 and d['last_even_spin']==448
    assert d['replay_precision_bits']>d['precision_bits']
    cursor=F(0);total=F(0);records=[];sums={'compact':F(0),'remainder':F(0)};costs={k:F(0) for k in sums}
    for c in d['cells']:
        lo,hi=F(c['lo']),F(c['hi']);assert cursor==lo<hi<=1
        first,replay,retained=map(F,[c['first_upper_exact'],c['replay_upper_exact'],c['retained_upper_exact']])
        assert 0<=first and 0<=replay and retained==max(first,replay)
        cover=[(ell,co.coverage(ell,lo,hi)) for ell in range(42,450,2)]
        limiting,entry=min(cover,key=lambda x:x[1][0]);kap=entry[0];assert kap>0
        weight=dyadic_ceiling(1/(2*kap));assert weight>=1/(2*kap)
        cost=retained*weight;total+=cost;sums[c['region']]+=retained;costs[c['region']]+=cost;cursor=hi
        ranges=[v[3] for ell,v in cover if v[3] is not None];large=[ell for ell,v in cover if v[2]]
        records.append(dict(lo=str(lo),hi=str(hi),region=c['region'],A_upper_exact=str(retained),kappa_exact=str(kap),limiting_spin=limiting,reciprocal_weight_upper_exact=str(weight),cost_upper_exact=str(cost),compact_piece_included=entry[1],large_tau_spins=large,finite_tau_cell_index_ranges=ranges,spins_checked=204))
    assert cursor==1 and all(sums[k]==F(d['A_upper_exact'][k]) for k in sums)
    co.manifest[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    out=dict(status='PASS_BANDED_COERCIVITY_REAL_COST',witness_sha256=co.SHA,objective_exact=co.objective,first_even_spin=42,last_even_spin=448,energy_domain=['4','infinity'],cost_upper_exact=str(total),cost_upper_float=float(total),cost_by_real_region_exact={k:str(v) for k,v in costs.items()},cells=records,compact_kappa_exact=str(co.compact),large_tau_kappa_exact=str(co.large),finite_UV_pairs=co.references,sha256_manifest=co.manifest,seconds=time.time()-tic,reciprocal_rounding='1/(2*kappa) is rounded upward to40 significant binary bits with exact rational ceiling inequalities; all products and sums are exact dyadics.',cover_ledger='Each finite_tau_cell_index_ranges entry [ell,first,stop] refers to every touching closed paired UV cell at indices first through stop-1. Both closed neighbours are retained at exact seams. Compact and large-tau pieces are included as explicitly recorded. Kappa is the exact minimum across all204 spins and all overlapping pieces.',scope='Complete finite42..448 allenergy support upper from existing same-witness real-norm cells; no new integration. Requires separate finitehead and infinite-spin support upper.')
    target=Path(args.out) if args.out else base/'finite_banded_cost.json';target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k not in ['cells','sha256_manifest','cost_upper_exact']},indent=2))
if __name__=='__main__':main()
