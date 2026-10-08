"""All-spin support bound allowing positive V on bounded joint-tau bands.

For eta>0, disk unitarity implies
  h(U,V) <= U^2/(2 eta p) + 2 p (V/p+eta)_+.
On mu>=M, tau=l(l+1)/(mu-4), and even l>=L,
  sum_l integral_[a,b] p dnu
  <= (b-a)b I0(2sqrt((4/3)b))/(2 pi sqrt(1-4/M) L^2).
The latter uses Q>=l(l+1)/tau, the positive Legendre Taylor majorant,
and sum_all_integer_l>=L (2l+1)/(l^2(l+1)^2)=1/L^2.
All omitted tau>=256 and compact-energy regions must obey V/p<=-eta.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, pickle, sys, time
ROOT=Path(__file__).resolve().parents[2]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'research/tail_certification'),str(ROOT/'tmp/certification_deps')]
from certify_uv import ar,arb,ctx,errors,leg_coeffs,poly
from certify_compact import run as compact_run
OUT=ROOT/'output/full_direction_20261007/certification'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def up(q):
    assert q.is_finite();m,e=q.upper().man_exp();return F(int(m))*F(2)**int(e)
def uv_bounds(d,cases,L,M,step):
    coeff=leg_coeffs(L);d2=[coeff[k]*k*(k-1) for k in range(2,len(coeff))]
    ea,en=errors(d,M)
    A=-F(32,9)*d['alpha']-d['signed2']+ea
    B=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+en
    pp=poly(coeff,ar(F(4,3)*256));nmax=B+2*sum((v for v,t in d['correction']),F(0))
    big=ar(A)+ar(max(nmax,F(0)))/pp
    large={b:up(big+ar(b)) for b,delta in cases}
    cuts=[F(0)]
    while cuts[-1]<8:cuts.append(min(F(8),cuts[-1]+step))
    cuts.extend(F(i) for i in range(9,257));rho=ar(1-F(4,M)).sqrt();records=[]
    for lo,hi in zip(cuts,cuts[1:]):
        mcell=max(F(M),F(4)+F(L*(L+1))/hi)
        eac,enc=errors(d,mcell)
        Ac=-F(32,9)*d['alpha']-d['signed2']+eac
        Bc=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+enc
        al=ar(F(4,3)*lo);ah=ar(F(4,3)*hi)
        pl=poly(coeff,al);pu=(2*ah.sqrt()).bessel_i(0)
        cor=sum((min(2*v,hi*t) for v,t in d['correction']),F(0));num=Bc+cor
        ratio=ar(num)/(pl if num>=0 else pu)
        if lo==0 or lo>=8:xlo=arb(0)
        else:
            p2=poly(d2,al);pa=(2*ah.sqrt()).bessel_i(1)/ah.sqrt()
            neg=ar(F(3,2)*hi/(mcell-F(4,3)))*pa/pl
            raw=ar(lo*lo/F(2))*p2/pu-neg;xlo=raw if raw>0 else arb(0)
        bounds={str(b):str(up(ar(Ac)+ratio+ar(d['lx']+delta)*xlo+ar(b))) for b,delta in cases}
        mass=up(ar((hi-lo)*hi)*pu/(2*arb.pi()*rho*L**2))
        records.append(dict(tau_a=str(lo),tau_b=str(hi),p_mass_upper_exact=str(mass),V_over_p_upper_by_slope=bounds))
    return large,records

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--slopes',default='1/2,1,2,4');ap.add_argument('--ratio',default='18/5')
    ap.add_argument('--precision',type=int,default=280);ap.add_argument('--step',default='1/250')
    ap.add_argument('--output',default='positive_uv_tail.json');a=ap.parse_args()
    tic=time.time();ctx.prec=a.precision;L=42;M=49
    cases=[(F(b),F(a.ratio)*F(b)) for b in a.slopes.split(',')]
    frozen_path=ROOT/'output/certification_tail/fixed_tail_certificate.json'
    frozen=json.loads(frozen_path.read_text())
    assert frozen['status']=='All-spin/all-energy omitted-spin tail certified under the stated amplitude assumptions.'
    assert frozen['replayed_components_this_run'] and F(frozen['V_over_p_upper_rational'])==-F(31,50)
    hashes=dict(frozen['sha256_manifest']);hashes[str(frozen_path.relative_to(ROOT))]=sha(frozen_path)
    for p,h in hashes.items():assert sha(ROOT/p)==h,p
    source=pickle.loads((ROOT/'output/certification_tail/source_cache_v2_depth10.pkl').read_bytes())
    compact=pickle.loads((ROOT/'output/certification_tail/compact_source_cache.pkl').read_bytes())
    assert source['sha']==compact['witness_sha256']==frozen['witness_sha256']
    d=source['data'];assert d['correction_integrated_t'] and all(d['lx']+delta<0 for b,delta in cases)
    big,cells=uv_bounds(d,cases,L,M,F(a.step));A=F(frozen['A40_upper_exact']);results=[]
    eta_list=list(map(F,['1/1000','1/500','1/200','1/100','1/50','1/20','1/10','1/5','2/5','1/2','1','2']))
    for b,delta in cases:
        dd=dict(d);dd['lx']+=delta
        try:comp=compact_run(L,dd,compact,'.02')
        except RuntimeError as error:
            r=dict(b=str(b),delta_X4=str(delta),status='FAILED_COMPACT_ENVELOPE',reason=str(error))
            results.append(r);print(json.dumps(r),flush=True);continue
        compact_upper=max(up(arb(comp['low_energy_upper'])+ar(b)),up(arb(comp['middle_energy_upper'])+ar(b)))
        candidates=[]
        for eta in eta_list:
            if compact_upper+eta>0 or big[b]+eta>0:continue
            positive=F(0);active=[]
            for i,c in enumerate(cells):
                vp=F(c['V_over_p_upper_by_slope'][str(b)])
                if vp+eta>0:
                    positive+=2*(vp+eta)*F(c['p_mass_upper_exact']);active.append(i)
            total=A/(2*eta)+positive
            candidates.append(dict(eta_exact=str(eta),tail_support_upper_exact=str(total),tail_support_upper_float=float(total),
                real_term_exact=str(A/(2*eta)),positive_V_term_exact=str(positive),positive_cells=active))
        if candidates:
            best=min(candidates,key=lambda r:F(r['tail_support_upper_exact']))
            r=dict(b=str(b),delta_X4=str(delta),status='PASS',compact_V_over_p_upper_exact=str(compact_upper),
                large_tau_V_over_p_upper_exact=str(big[b]),best=best,eta_trials=candidates,compact=comp)
        else:r=dict(b=str(b),delta_X4=str(delta),status='FAILED_UNBOUNDED_REGION_SIGN',compact=comp,large_tau_V_over_p_upper_exact=str(big[b]))
        results.append(r);print(json.dumps(dict(b=str(b),status=r['status'],tail_support_upper_float=r.get('best',{}).get('tail_support_upper_float'))),flush=True)
    for p in [Path(__file__),ROOT/'research/directional_certification/certify_directional_tail.py']:
        hashes[str(p.relative_to(ROOT))]=sha(p)
    result=dict(status='PASS' if all(r['status']=='PASS' for r in results) else 'PARTIAL',witness_sha256=frozen['witness_sha256'],
        domain='All even ell>=42 and all mu>4',precision_bits=ctx.prec,results=results,uv_cells=cells,
        all_spin_real_norm_upper_exact=str(A),classification='Validated tail support; finite heads remain separate.',
        proof=__doc__,seconds=time.time()-tic,sha256_manifest=hashes)
    (OUT/a.output).write_text(json.dumps(result,indent=2));print(json.dumps(dict(seconds=result['seconds'])))
if __name__=='__main__':main()
