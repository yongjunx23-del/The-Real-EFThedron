"""Combine fixed finite-real and finite-V interval certificates locally.

No kernel values or quadrature nodes are recomputed. On each beta-cell the
certified real mass is divided by twice the weakest certified coercivity
over every overlapping mu/w cell. All partition comparisons are rational.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tmp/certification_deps'))
from flint import arb,ctx,fmpq

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dyad(p):return F(int(p[0]))*F(2)**int(p[1])
def ar(f):return arb(fmpq(f.numerator,f.denominator))
def verify_partition(cells,lo,hi):
    assert cells[0][0]==lo and cells[-1][1]==hi
    assert all(a<b for a,b in cells)
    assert all(b==c for (a,b),(c,d) in zip(cells,cells[1:]))

def main():
    ctx.prec=280
    ap='output/certification_tail/finite_A8_40_tighter.json'
    vp='output/certification_head/V8_40_kappa14_exact.json'
    rp='output/certification_head/V8_40_kappa14_independent_replay.json'
    aa=json.loads((ROOT/ap).read_text());vv=json.loads((ROOT/vp).read_text())
    replay=json.loads((ROOT/rp).read_text())
    assert replay['status']=='PASS' and replay['certificate_sha256']==sha(vp)
    assert replay['every_spin_finite_and_below_threshold']
    assert aa['witness_sha256']==vv['witness_sha256']
    assert aa['spins']==vv['spins']==list(range(8,42,2))
    assert replay['spins']==vv['spins'] and replay['witness_sha256']==vv['witness_sha256']
    for path,digest in replay['sha256_manifest'].items():assert sha(path)==digest
    cells=sorted(aa['partition'],key=lambda c:F(c['a']))
    verify_partition([(F(c['a']),F(c['b'])) for c in cells],F(0),F(1))
    assert sum((dyad(c['mass_upper']) for c in cells),F(0))<=dyad(aa['A_upper_dyadic'])<F(47,200)
    segments=[]
    for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
        part=sorted((F(c['a']),F(c['b'])) for c in vv['leaves'] if c['variable']==kind)
        verify_partition(part,lo,hi)
    for i,c in enumerate(vv['leaves']):
        k=-dyad(c['upper_dyadic']);assert k>14
        a,b=F(c['a']),F(c['b'])
        if c['variable']=='mu':lo,hi=a,b
        else:lo,hi=4+1/b,(4+1/a if a else None)
        segments.append((lo,hi,k,i))
    segments.sort()
    assert segments[0][0]==4 and segments[-1][1] is None
    assert all(b==c for (a,b,k,i),(c,d,m,j) in zip(segments,segments[1:]))
    j=0;total=arb(0);records=[]
    for c in cells:
        a,b=F(c['a']),F(c['b']);lo=4/(1-a*a);hi=4/(1-b*b) if b<1 else None
        while segments[j][1] is not None and segments[j][1]<=lo:j+=1
        stop=j
        while stop+1<len(segments) and (hi is None or segments[stop+1][0]<hi):stop+=1
        k,idx=min((segments[i][2],segments[i][3]) for i in range(j,stop+1))
        mass=dyad(c['mass_upper']);assert mass>=0
        term=ar(mass)/(2*ar(k));assert term.is_finite()
        total+=term.upper()
        records.append({'a':str(a),'b':str(b),'worst_v_leaf_index':idx,'coercivity_exact':str(k),'support_upper_dyadic':list(map(str,term.upper().man_exp()))})
    bound=total.upper();assert bound<ar(F(47,5600))
    inputs=[ap,vp,rp,'research/tail_certification/combine_finite_local.py']
    out=dict(status='PASS',claim='Complete disk support of finite even spins8..40 bounded using local certified coercivity.',witness_sha256=aa['witness_sha256'],spins=aa['spins'],precision_bits=ctx.prec,support_upper=str(bound),support_upper_dyadic=list(map(str,bound.man_exp())),cells=len(records),proof='For each real-norm cell, every overlapping spectral V cell is checked. V<=-k*p and h<=U^2/(2*k*p) hold throughout that cell. Sum all directed cell bounds.',sha256={p:sha(p) for p in inputs},partition=records)
    op=ROOT/'output/certification_head/block8_40_local_coercivity.json';op.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='partition'},indent=2))

if __name__=='__main__':main()
