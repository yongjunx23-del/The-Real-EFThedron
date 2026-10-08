"""Reuse the fixed-source analytic all-spin tail envelope for (1,b), X4+d.

All old files are read-only. Exact source data are shared because only the
objective g2 coefficient and the homogeneous X4 multiplier change.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,pickle,hashlib,argparse,time
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'research/tail_certification'),str(ROOT/'tmp/certification_deps')]
from certify_uv import ar,arb,ctx,errors,leg_coeffs,poly
from certify_compact import run as compact_run

CASES=[(F(s),F(0) if F(s)<=0 else F(32,9)*F(s)) for s in ['-1','-1/2','-1/4','-1/10','0','1/10','1/4','1/2']]
def dyad(x):
    assert x.is_finite();m,e=x.upper().man_exp();return F(int(m))*F(2)**int(e)
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def key(b):return str(b).replace('-','m').replace('/','_')

def uv_family(d,L=42,M=49,step=F(1,250)):
    """Same inequalities as frozen certify_uv.certify, with affine b,d added."""
    assert L>=42 and L%2==0 and M>48 and d['correction_integrated_t']
    assert all(d['lx']+delta<0 for b,delta in CASES)
    coeff=leg_coeffs(L);d2=[coeff[k]*k*(k-1) for k in range(2,len(coeff))]
    ea,en=errors(d,M)
    A=-F(32,9)*d['alpha']-d['signed2']+ea
    B=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+en
    pp=poly(coeff,ar(F(4,3)*256))
    Nmax=B+2*sum((v for v,t in d['correction']),F(0))
    big=ar(A)+ar(max(Nmax,F(0)))/pp
    worst={b:dyad(big+ar(b)) for b,delta in CASES};where={b:'tau>=256' for b,delta in CASES}
    cuts=[F(0)]
    while cuts[-1]<8:cuts.append(min(F(8),cuts[-1]+step))
    cuts.extend(F(i) for i in range(9,257))
    saved=[]
    for lo,hi in zip(cuts,cuts[1:]):
        mcell=max(F(M),F(4)+F(L*(L+1))/hi)
        eac,enc=errors(d,mcell)
        Ac=-F(32,9)*d['alpha']-d['signed2']+eac
        Bc=-F(280,9)*d['alpha']+32*sum((v*(-1)**j for j,v in d['harm']),F(0))-d['signed1']+enc
        al=ar(F(4,3)*lo);ah=ar(F(4,3)*hi)
        pl=poly(coeff,al);pu=(2*ah.sqrt()).bessel_i(0)
        cor=sum((min(2*v,hi*t) for v,t in d['correction']),F(0))
        num=Bc+cor;ratio=ar(num)/(pl if num>=0 else pu)
        if lo==0 or lo>=8:xlo=arb(0)
        else:
            p2=poly(d2,al);pa=(2*ah.sqrt()).bessel_i(1)/ah.sqrt()
            neg=ar(F(3,2)*hi/(mcell-F(4,3)))*pa/pl
            raw=ar(lo*lo/F(2))*p2/pu-neg;xlo=raw if raw>0 else arb(0)
        for b,delta in CASES:
            value=ar(Ac)+ratio+ar(d['lx']+delta)*xlo+ar(b)
            bound=dyad(value)
            if bound>worst[b]:worst[b]=bound;where[b]=[str(lo),str(hi)]
        saved.append({'tau_a':str(lo),'tau_b':str(hi),'cross_lower':str(xlo),'source_upper':str((ar(Ac)+ratio).upper())})
    return {b:dict(upper_exact=str(worst[b]),worst_cell=where[b],large_tau_upper_exact=str(dyad(big+ar(b)))) for b,delta in CASES},saved

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,default=42);ap.add_argument('--output',type=Path,default=ROOT/'output/directional_certification/tail_family.json');args=ap.parse_args()
    ctx.prec=200;tic=time.time()
    frozen=json.loads((ROOT/'output/certification_tail/fixed_tail_certificate.json').read_text())
    for p,h in frozen['sha256_manifest'].items():assert sha(p)==h
    source=pickle.loads((ROOT/'output/certification_tail/source_cache_v2_depth10.pkl').read_bytes())
    cd=pickle.loads((ROOT/'output/certification_tail/compact_source_cache.pkl').read_bytes())
    assert source['sha']==cd['witness_sha256']==frozen['witness_sha256'];d=source['data']
    uv,cells=uv_family(d,L=args.L)
    results=[]
    for b,delta in CASES:
        if b<=0 and args.L==42:
            k=F(31,50)-b
            r=dict(b=str(b),delta_X4=str(delta),status='PASS',tail_kappa_exact=str(k),method='Frozen global coercivity shifted by negative g2 direction.')
        else:
            dd=dict(d);dd['lx']+=delta
            compact=compact_run(args.L,dd,cd,'.02')
            bounds=[F(uv[b]['upper_exact']),dyad(arb(compact['low_energy_upper'])+ar(b)),dyad(arb(compact['middle_energy_upper'])+ar(b))]
            upper=max(bounds);k=-upper
            r=dict(b=str(b),delta_X4=str(delta),status='PASS' if k>0 else 'FAILED_ENVELOPE',tail_kappa_exact=str(k),compact=compact,uv=uv[b],method='Replayed compact and unbounded UV analytic envelopes with modified negative crossing coefficient; exact b added.')
        if r['status']=='PASS' and args.L==42:
            A=F(frozen['A40_upper_exact']);cost=A/(2*F(r['tail_kappa_exact']))
            r['A_tail_upper_exact']=str(A);r['tail_support_upper_exact']=str(cost);r['tail_support_upper_float']=float(cost)
        results.append(r);print(json.dumps({k:v for k,v in r.items() if k not in ['compact','uv']}),flush=True)
    files=['research/directional_certification/certify_directional_tail.py','output/certification_tail/fixed_tail_certificate.json']
    hashes={p:sha(p) for p in files};hashes.update(frozen['sha256_manifest'])
    out=dict(witness_sha256=frozen['witness_sha256'],L=args.L,domain=f'all even ell>={args.L}, all mu>4',precision_bits=ctx.prec,seconds=time.time()-tic,results=results,uv_cells=cells,sha256_manifest=hashes)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
