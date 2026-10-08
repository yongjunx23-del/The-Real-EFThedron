"""Rigorous finite real norm weighted by certified per-spin UV coercivity.

Every input proof must already cover its asserted continuum domain. This
module independently checks the exact partition, paired replay and hashes,
then uses the weakest kappa touching each real-norm integration rectangle.
"""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left
import argparse,json,hashlib,time,heapq
from real220_finite import Finite,val
from real220 import ROOT,arb,ctx,exact,square,tight,endpoint

SHA='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class Coercivity:
    def __init__(self,descriptor):
        self.manifest={};d=self.read(descriptor);assert d['witness_sha256']==SHA
        self.compact=F(d['compact_kappa_exact']);self.large=F(d['large_tau_kappa_exact'])
        assert self.compact>0 and self.large>0
        assert d['objective_exact']==['-1','0']
        assert d['tau_definition']=='ell*(ell+1)/(mu-4)'
        compact=self.read(ROOT/d['compact_proof'])
        assert compact['status']=='PASS_COMPLETE_COMPACT_COERCIVITY'
        assert compact['witness_sha256']==SHA and compact['objective_exact']==['-1','0']
        assert compact['minimum_even_spin']==42 and compact['mu_domain_exact']==['4','49']
        assert compact['lower_energy_endpoint_open'] and self.compact<=F(compact['kappa_exact'])
        large=self.read(ROOT/d['large_tau_proof'])
        assert large['status']=='PASS_LEFT_LARGE_TAU' and large['witness_sha256']==SHA
        assert large['objective_exact']==['-1','0'] and large['lmin']<=108
        assert large['mu_domain']==['49','infinity'] and large['tau_domain']==['256','infinity']
        assert self.large<=F(large['kappa_exact'])<=-F(large['directed_upper_exact'])
        self.read(ROOT/'output/all_direction_certification_20261007/tail/left/refined/finite_source_N12_K80.json')
        self.cells={}
        for p,q in d['finite_UV_pairs']:
            first=self.read(ROOT/p);replay=self.read(ROOT/q)
            assert first['witness_sha256']==replay['witness_sha256']==SHA
            assert first['status']==replay['status']=='PASS_LEFT_FINITE_UV'
            assert first['objective_exact']==replay['objective_exact']==['-1','0']
            assert not first['failed'] and not replay['failed']
            assert replay['precision_bits']>first['precision_bits']
            assert replay['reference_certificate_sha256']==sha(ROOT/p)
            assert replay['spins']==first['spins']
            assert [r['spin']for r in first['records']]==first['spins']
            assert [r['spin']for r in replay['records']]==replay['spins']
            for r,s in zip(first['records'],replay['records']):
                ell=r['spin'];assert ell==s['spin'] and ell not in self.cells
                assert r['status']==s['status']=='PASS'
                assert len(r['cells'])==len(s['cells'])
                cells=[];cursor=F(0)
                for a,b in zip(r['cells'],s['cells']):
                    lo,hi=F(a['a']),F(a['b']);assert lo==F(b['a'])==cursor<hi==F(b['b'])
                    upper=max(F(a['upper_exact']),F(b['upper_exact']));assert upper<0
                    cells.append((lo,hi,-upper));cursor=hi
                assert cursor==min(F(256),F(ell*(ell+1),45))==F(r['tau_end_exact'])==F(s['tau_end_exact'])
                self.cells[ell]=cells
        assert sorted(self.cells)==list(range(42,449,2))
        self.ends={ell:[c[1]for c in cells]for ell,cells in self.cells.items()}

    def read(self,path):
        path=Path(path);self.manifest[str(path.relative_to(ROOT))]=sha(path);d=json.loads(path.read_text())
        for p,h in d.get('sha256_manifest',{}).items():
            file=Path(p);file=file if file.is_absolute()else ROOT/file
            assert sha(file)==h,p;self.manifest[str(file.relative_to(ROOT))]=h
        return d

    def kappa(self,ell,lo,hi):
        assert 0<=lo<hi<=1;out=[];T=ell*(ell+1)
        # Exact comparison avoids rounding the beta threshold sqrt(45/49).
        if lo*lo<F(45,49):out.append(self.compact)
        if hi*hi>F(45,49):
            ta=F(T)*(1-hi*hi)/(4*hi*hi)
            tb=min(F(T,45),F(T)*(1-lo*lo)/(4*lo*lo))if lo else F(T,45)
            if tb>256:
                assert ell>=108;out.append(self.large)
            tb=min(tb,F(256))
            if ta<=tb:
                cs=self.cells[ell];i=bisect_left(self.ends[ell],ta)
                assert i<len(cs) and cs[i][0]<=ta<=cs[i][1]
                while i<len(cs)and cs[i][0]<=tb:
                    out.append(cs[i][2]);i+=1
        assert out
        return min(out)

class WeightedFinite(Finite):
    def density_weighted(self,lo,hi,weights):
        assert len(weights)==len(self.spins) and all(w>0 for w in weights)
        mid=(lo+hi)/2;h=(hi-lo)/2;beta=exact(lo).union(exact(hi));z=beta*beta;e=1-z;bm=exact(mid);dx=arb(0,exact(h))
        fields=[(val(p,bm),val(dp,bm),val(ddp,beta))for p,dp,ddp in self.zfields]
        if any(not v.is_finite()for r in fields for v in r):return arb.pos_inf()
        sa=4/(1-lo*lo);sb=4/(1-hi*hi)if hi<1 else None;strips=[]
        for (A,B),terms in self.groups.items():
            if sa>=B or(sb is not None and sb<=A):continue
            assert hi<1
            sc=(4/e).intersection(exact(A).union(exact(B)));en=(2*(sc-exact(A))/exact(B-A)-1).intersection(arb(0,1));taper=(1-square(en)).abs_upper()
            l2=sum((square(tight(self.energy[A,B,j],en)).upper()/(2*j+1)for j in {j for j,_,_ in terms}),arb(0))
            dd=sc-4;T=exact(F(4,3));fac=2*dd/(arb.pi()*(T*(dd-T)).sqrt())
            pref=arb.pi()/exact(B-A)*taper*l2.sqrt()*beta*(1+z)**5/(16*e**4)
            strips.append((pref.abs_upper(),fac.upper()))
        powers=[(e/3)**k*z**(self.K-k)for k in range(self.K+1)];numerator=z**(self.K-1);ss=arb(0)
        for (ell,B),weight in zip(self.spins,weights):
            f0=sum((B[k]*fields[k][0]for k in range(4)),arb(0));f1=sum((B[k]*fields[k][1]for k in range(4)),arb(0));f2=sum((B[k]*fields[k][2]for k in range(4)),arb(0))
            value=(f0+f1*dx+f2*dx*dx/2).abs_upper()
            for pref,fac in strips:value+=pref*min(arb(1),(fac/exact(F(2*ell+1,2))).sqrt().upper(),key=float)
            den=sum((B[k]*powers[k]for k in range(self.K+1)),arb(0))
            ss+=exact(weight)*(2*ell+1)*value*value*numerator/den
        return arb(65536)/(27*arb.pi())*e**3*(2+z)**3/(1+z)**10*ss

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--descriptor',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--cells',type=int,default=4096);ap.add_argument('--bits',type=int,default=256);ap.add_argument('--replay-bits',type=int,default=384)
    args=ap.parse_args();tic=time.time();co=Coercivity(args.descriptor.resolve());path=ROOT/'output/full_direction_20261007/left/angle180_final.npz';assert sha(path)==SHA
    ctx.prec=args.bits;model=WeightedFinite(path,42,448,20)
    for lo,hi in [(F(1,100),F(11,1000)),(F(8,10),F(81,100)),(F(999,1000),F(9991,10000))]:
        a=model.density(lo,hi);b=model.density_weighted(lo,hi,[F(1)]*len(model.spins));assert a.overlaps(b)
    heap=[];total=F(0);ident=0
    def integrate(lo,hi):
        weights=[1/(2*co.kappa(ell,lo,hi))for ell,_ in model.spins]
        v=model.density_weighted(lo,hi,weights)
        return max(F(0),endpoint(v))*(hi-lo)if v.is_finite()else None
    def insert(lo,hi):
        nonlocal ident,total
        mass=integrate(lo,hi)
        if mass is None:
            assert hi-lo>F(1,2**30);m=(lo+hi)/2;insert(lo,m);insert(m,hi);return
        total+=mass;heapq.heappush(heap,(-float(mass),ident,lo,hi,mass));ident+=1
    for i in range(64):insert(F(i,64),F(i+1,64))
    while len(heap)<args.cells:
        _,_,lo,hi,mass=heapq.heappop(heap);total-=mass;m=(lo+hi)/2;insert(lo,m);insert(m,hi)
        if len(heap)%512==0:print('cells',len(heap),'cost',float(total),'seconds',time.time()-tic,flush=True)
    ctx.prec=args.replay_bits;assert args.replay_bits>args.bits;model=WeightedFinite(path,42,448,20)
    records=[];retained=F(0);cursor=F(0)
    for i,(_,_,lo,hi,mass)in enumerate(sorted(heap,key=lambda r:r[2])):
        assert cursor==lo<hi<=1;v=integrate(lo,hi);assert v is not None
        keep=max(mass,v);retained+=keep;cursor=hi
        records.append(dict(lo=str(lo),hi=str(hi),first_upper_exact=str(mass),replay_upper_exact=str(v),retained_upper_exact=str(keep)))
        if (i+1)%1024==0:print('replayed',i+1,'seconds',time.time()-tic,flush=True)
    assert cursor==1
    for p in [Path(__file__),Path(__file__).with_name('real220.py'),Path(__file__).with_name('real220_finite.py'),path]:co.manifest[str(p.relative_to(ROOT))]=sha(p)
    out=dict(status='PASS_LOCAL_COERCIVITY_WEIGHTED_REAL_COST',witness_sha256=SHA,objective_exact=['-1','0'],first_even_spin=42,last_even_spin=448,
        energy_domain=['4','infinity'],cost_upper_exact=str(retained),cost_upper_float=float(retained),precision_bits=args.bits,replay_precision_bits=args.replay_bits,
        compact_kappa_exact=str(co.compact),large_tau_kappa_exact=str(co.large),cells=records,sha256_manifest=co.manifest,seconds=time.time()-tic,
        scope='Finite spins42..448, all energies. Complete directional certificate additionally requires allspin>=450 and finitehead0..40. Local weights use the weakest directed V bound touching each closed beta cell; no interpolation.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['cells','sha256_manifest','cost_upper_exact']},indent=2))
if __name__=='__main__':main()
