"""Exact partition merge, using each spin's own real mass and coercivity."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[3]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'tmp/certification_deps'))
from flint import arb,ctx,fmpq
CASES=[(F(4),F(72,5))]

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dyad(pair):return F(int(pair[0]))*F(2)**int(pair[1])
def ar(q):return arb(fmpq(q.numerator,q.denominator))
def upper_pair(x):
    assert x.is_finite();return list(map(str,x.upper().man_exp()))
def partition(parts,lo,hi):
    assert parts[0][0]==lo and parts[-1][1]==hi
    assert all(a<b for a,b in parts)
    assert all(b==c for (a,b),(c,d) in zip(parts,parts[1:]))

def main():
    ctx.prec=280;tic=time.time()
    replaypath=ROOT/'output/competitive_20261007/robust/b4/finite_v_replay.json'
    replay=json.loads(replaypath.read_text())
    assert replay['status']=='PASS_FRESH_SIGN_REPLAY_AND_EXACT_REUSE'
    assert F(replay['maximum_combined_upper_overshoot_exact'])<F(1,10000)
    assert replay['certificate_sha256']==sha('output/competitive_20261007/robust/b4/finite_v_b4.json')
    files=['output/directional_certification/finite_real_vectors.json','output/competitive_20261007/robust/b4/finite_v_b4.json']
    aa,vv=[json.loads((ROOT/p).read_text()) for p in files]
    assert aa['status']==vv['status']=='PASS' and aa['witness_sha256']==vv['witness_sha256']
    assert aa['spins']==vv['spins']==list(range(8,42,2));spins=aa['spins']
    hashes={p:sha(p) for p in files}
    hashes[str(replaypath.relative_to(ROOT))]=sha(str(replaypath.relative_to(ROOT)))
    for p,h in replay['sha256_manifest'].items():assert sha(p)==h;hashes[p]=h
    for data in (aa,vv):
        for p,h in data['sha256_manifest'].items():assert sha(p)==h;hashes[p]=h
    acells=sorted(aa['partition'],key=lambda c:F(c['a']))
    partition([(F(c['a']),F(c['b'])) for c in acells],F(0),F(1))
    for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
        ranges=sorted((F(c['a']),F(c['b'])) for c in vv['partition'] if c['variable']==kind);partition(ranges,lo,hi)
    segments=[];failed={}
    for c in vv['partition']:
        assert set(map(int,c['upper_by_spin']))==set(spins)
        a,b=F(c['a']),F(c['b']);assert c['variable'] in ('mu','w')
        lo,hi=(a,b) if c['variable']=='mu' else (4+1/b,4+1/a if a else None)
        kappas={}
        for slope,delta in CASES:
            kappas[slope]={}
            for l in spins:
                bounds=c['upper_by_spin'][str(l)]
                k=-(dyad(bounds['V_upper'])+slope+delta*dyad(bounds['X4_upper'])+F(1,10000))
                if k<=0 and slope not in failed:
                    failed[slope]=dict(variable=c['variable'],a=str(a),b=str(b),spin=l,coercivity_exact=str(k))
                kappas[slope][l]=k
        segments.append((lo,hi,kappas))
    segments.sort(key=lambda s:s[0]);assert segments[0][0]==4 and segments[-1][1] is None
    assert all(a[1]==b[0] for a,b in zip(segments,segments[1:]))
    active=[(b,d) for b,d in CASES if b not in failed]
    totals={b:{l:arb(0) for l in spins} for b,d in active};j=0
    for i,c in enumerate(acells):
        a,b=F(c['a']),F(c['b']);lo=4/(1-a*a);hi=4/(1-b*b) if b<1 else None
        while segments[j][1] is not None and segments[j][1]<=lo:j+=1
        stop=j
        while stop+1<len(segments) and (hi is None or segments[stop+1][0]<hi):stop+=1
        assert set(map(int,c['mass_upper_by_spin']))==set(spins)
        masses={l:dyad(c['mass_upper_by_spin'][str(l)]) for l in spins}
        assert all(m>=0 for m in masses.values())
        for slope,delta in active:
            for l in spins:
                k=min(segments[z][2][slope][l] for z in range(j,stop+1))
                value=ar(masses[l])/(2*ar(k));assert value.is_finite()
                totals[slope][l]+=value.upper()
        if (i+1)%4000==0:print(json.dumps({'cells':i+1,'seconds':time.time()-tic}),flush=True)
    results=[]
    for slope,delta in CASES:
        if slope in failed:
            results.append(dict(b=str(slope),delta_X4=str(delta),status='FAILED_COERCIVITY',first_failure=failed[slope]));continue
        byspin={str(l):upper_pair(totals[slope][l]) for l in spins}
        total=sum((ar(dyad(v)) for v in byspin.values()),arb(0)).upper();up=upper_pair(total)
        results.append(dict(b=str(slope),delta_X4=str(delta),status='PASS',finite_support_upper_exact=str(dyad(up)),finite_support_upper_dyadic=up,finite_support_upper_float=float(total),support_upper_by_spin_dyadic=byspin))
    here='research/competitive_20261007/robust/combine_b4_finite_replayed.py';hashes[here]=sha(here)
    out=dict(status='PASS',witness_sha256=aa['witness_sha256'],domain='all mu>4 and finite even spins8..40',spins=spins,precision_bits=ctx.prec,real_cells=len(acells),spectral_cells=len(segments),results=results,seconds=time.time()-tic,independent_replay_uniform_V_slack_exact='1/10000',proof='U is unchanged. Added1/10000 to each combined V upper to dominate the independently replayed endpoints. Every source leaf bounds Vnew/p by V_upper+b+d*X4_upper, with d>=0. For each spin and each real beta-cell use the minimum coercivity over every exact-overlapping spectral leaf; sum real_mass/(2*coercivity) with directed rounding.',sha256_manifest=hashes)
    path=ROOT/'output/competitive_20261007/robust/b4/finite_block_replayed.json';path.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k!='sha256_manifest'},indent=2))

if __name__=='__main__':main()
