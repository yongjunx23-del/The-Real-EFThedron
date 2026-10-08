"""Exact-overlap merge of per-spin real norms and replayed b6 coercivity."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib
ROOT=Path(__file__).resolve().parents[4];sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'tmp/certification_deps'))
from flint import arb,ctx,fmpq
OUT=ROOT/'output/competitive_20261007/upper_arc/b6'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(pair):return F(int(pair[0]))*F(2)**int(pair[1])
def ar(q):return arb(fmpq(q.numerator,q.denominator))
def upper_pair(q):
    assert q.is_finite();return list(map(str,q.upper().man_exp()))
def partition(parts,lo,hi):
    assert parts[0][0]==lo and parts[-1][1]==hi and all(a<b for a,b in parts)
    assert all(b==c for(a,b),(c,d)in zip(parts,parts[1:]))
def main():
    ctx.prec=320;tic=time.time();hashes={}
    def read(p):
        q=json.loads(p.read_text());hashes[str(p.relative_to(ROOT))]=sha(p)
        for f,h in q['sha256_manifest'].items():assert sha(ROOT/f)==h,f;hashes[f]=h
        return q
    aa=read(ROOT/'output/directional_certification/finite_real_vectors.json');vv=read(OUT/'finite_v_replayed.json')
    assert aa['status']=='PASS' and vv['status']=='PASS_FRESH_ALL_CELL_REPLAY_AND_OUTWARD_MAXIMA'
    assert aa['witness_sha256']==vv['witness_sha256'] and aa['spins']==vv['spins']==list(range(8,42,2))
    assert vv['direction_exact']==['1','6'] and vv['crossing_increment_exact']=='107/5'
    spins=aa['spins'];acells=sorted(aa['partition'],key=lambda c:F(c['a']))
    partition([(F(c['a']),F(c['b']))for c in acells],F(0),F(1))
    for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
        partition(sorted((F(c['a']),F(c['b']))for c in vv['partition']if c['variable']==kind),lo,hi)
    segments=[]
    for c in vv['partition']:
        a,b=F(c['a']),F(c['b']);assert c['variable']in('mu','w')
        lo,hi=(a,b)if c['variable']=='mu'else(4+1/b,4+1/a if a else None)
        assert set(map(int,c['combined_upper_exact_by_spin']))==set(spins)
        k={l:-F(c['combined_upper_exact_by_spin'][str(l)])for l in spins};assert all(v>0 for v in k.values())
        segments.append((lo,hi,k))
    segments.sort(key=lambda c:c[0]);assert segments[0][0]==4 and segments[-1][1]is None
    assert all(a[1]==b[0]for a,b in zip(segments,segments[1:]))
    totals={l:arb(0)for l in spins};j=0
    for i,c in enumerate(acells):
        a,b=F(c['a']),F(c['b']);lo=4/(1-a*a);hi=4/(1-b*b)if b<1 else None
        while segments[j][1]is not None and segments[j][1]<=lo:j+=1
        stop=j
        while stop+1<len(segments)and(hi is None or segments[stop+1][0]<hi):stop+=1
        assert set(map(int,c['mass_upper_by_spin']))==set(spins)
        for l in spins:
            mass=dyad(c['mass_upper_by_spin'][str(l)]);assert mass>=0
            k=min(segments[z][2][l]for z in range(j,stop+1));v=ar(mass)/(2*ar(k));assert v.is_finite();totals[l]+=v.upper()
    byspin={str(l):upper_pair(totals[l])for l in spins}
    total=sum((ar(dyad(v))for v in byspin.values()),arb(0)).upper();up=upper_pair(total)
    hashes[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    out=dict(status='PASS_FINITE_SUPPORT_8_40',witness_sha256=aa['witness_sha256'],direction_exact=['1','6'],crossing_increment_exact='107/5',
      domain='All mu>4 and every even spin8..40',spins=spins,precision_bits=ctx.prec,real_cells=len(acells),spectral_cells=len(segments),
      finite_support_upper_exact=str(dyad(up)),finite_support_upper_float=float(total),support_upper_by_spin_dyadic=byspin,
      seconds=time.time()-tic,sha256_manifest=hashes,
      proof='U is unchanged from old150. Each spectral leaf uses the maximum of fresh320/fresh380 direct V/p upper bounds. For each real beta-cell/spin, minimize positive coercivity over all spectral leaves with positive-length overlap; omitted closed-boundary intersections are measure zero. Then H<=U^2/(2*kappa*p), summed with directed Arb320 arithmetic. The final w=0 leaf includes mu infinity.')
    (OUT/'finite_block.json').write_text(json.dumps(out,indent=2));print('PASS',float(total),time.time()-tic,flush=True)
if __name__=='__main__':main()
