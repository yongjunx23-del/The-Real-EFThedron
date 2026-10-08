"""Independent partition/merge audit and sampled higher-precision vector replay.

This checks all stored inequalities and the exact partition merge. Fresh
kernel evaluations are sampled, not a claim of a second full vector replay.
"""
from pathlib import Path
from fractions import Fraction as F
from bisect import bisect_left, bisect_right
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent))
import certify_finite_vectors as m

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def q(pair):return F(int(pair[0]))*F(2)**int(pair[1])
def upper(v):
    assert v.is_finite()
    return q(v.upper().man_exp())
def check_partition(parts,lo,hi):
    assert parts[0][0]==lo and parts[-1][1]==hi
    assert all(a<b for a,b in parts)
    assert all(b==c for (a,b),(c,d) in zip(parts,parts[1:]))
def select(n,count):return sorted(set([0,n-1]+[(n-1)*i//(count-1) for i in range(count)]))

def main():
    tic=time.time();m.ctx.prec=340
    paths=['output/directional_certification/finite_real_vectors.json',
           'output/directional_certification/finite_v_vectors.json',
           'output/directional_certification/finite_block_family.json']
    aa,vv,combined=[json.loads((ROOT/p).read_text()) for p in paths]
    hashes={p:sha(p) for p in paths}
    for data in (aa,vv,combined):
        assert data['status']=='PASS'
        for p,h in data['sha256_manifest'].items():assert sha(p)==h;hashes[p]=h
    spins=list(range(8,42,2));assert aa['spins']==vv['spins']==combined['spins']==spins
    assert aa['witness_sha256']==vv['witness_sha256']==combined['witness_sha256']==m.EXPECTED
    acells=sorted(aa['partition'],key=lambda c:F(c['a']))
    check_partition([(F(c['a']),F(c['b'])) for c in acells],F(0),F(1))
    vparts={kind:sorted([c for c in vv['partition'] if c['variable']==kind],key=lambda c:F(c['a'])) for kind in ['mu','w']}
    assert sum(map(len,vparts.values()))==len(vv['partition'])
    check_partition([(F(c['a']),F(c['b'])) for c in vparts['mu']],F(4),F(101))
    check_partition([(F(c['a']),F(c['b'])) for c in vparts['w']],F(0),F(1,97))
    segments=[]
    for c in vv['partition']:
        assert set(c['upper_by_spin'])==set(map(str,spins))
        a,b=F(c['a']),F(c['b'])
        lo,hi=(a,b) if c['variable']=='mu' else (4+1/b,4+1/a if a else None)
        segments.append((lo,hi,c))
    segments.sort(key=lambda c:c[0])
    assert segments[0][0]==4 and segments[-1][1] is None
    assert all(a[1]==b[0] for a,b in zip(segments,segments[1:]))
    starts=[c[0] for c in segments];ends=[c[1] for c in segments[:-1]]
    cases=[(F(r['b']),F(r['delta_X4']),r) for r in combined['results']]
    assert [(b,d) for b,d,r in cases]==[(F(s),max(F(0),4*F(s))) for s in ['-1','-1/2','-1/4','-1/10','0','1/10','1/4','1/2']]
    kappas={};active=[]
    for b,d,r in cases:
        ks=[{l:-(q(c['upper_by_spin'][str(l)]['V_upper'])+b+d*q(c['upper_by_spin'][str(l)]['X4_upper'])) for l in spins} for lo,hi,c in segments]
        valid=all(v>0 for row in ks for v in row.values())
        assert valid==(r['status']=='PASS')
        if valid:kappas[b]=ks;active.append((b,d,r))
    totals={b:{l:m.arb(0) for l in spins} for b,d,r in active}
    for c in acells:
        assert set(c['mass_upper_by_spin'])==set(map(str,spins))
        a,b=F(c['a']),F(c['b']);lo=4/(1-a*a);hi=4/(1-b*b) if b<1 else None
        # Independent binary searches implement interior overlap. Endpoints
        # are still covered because both source partitions are closed.
        first=bisect_right(ends,lo)
        last=len(segments) if hi is None else bisect_left(starts,hi)
        assert first<last
        for slope,d,r in active:
            for l in spins:
                mass=q(c['mass_upper_by_spin'][str(l)]);assert mass>=0
                k=min(row[l] for row in kappas[slope][first:last]);assert k>0
                totals[slope][l]+= (m.ar(mass)/(2*m.ar(k))).upper()
    comparisons=[]
    for b,d,r in active:
        fresh=upper(sum((v.upper() for v in totals[b].values()),m.arb(0)))
        stored=F(r['finite_support_upper_exact'])
        assert abs(fresh-stored)<F(1,10**12)
        cap=(max(fresh,stored)*10**12).__ceil__()*F(1,10**12)+F(1,10**12)
        comparisons.append(dict(b=str(b),delta_X4=str(d),fresh_upper_exact=str(fresh),stored_upper_exact=str(stored),common_decimal_cap=str(cap),common_decimal_cap_float=float(cap)))
    print(json.dumps({'event':'merge_pass','cases':len(comparisons),'seconds':time.time()-tic}),flush=True)

    bank=m.FiniteV();vchecks=[];max_v_excess=F(0);max_x_excess=F(0)
    for kind,cells in vparts.items():
        for i in select(len(cells),65):
            c=cells[i];a,b=F(c['a']),F(c['b']);x=m.arb(m.ar((a+b)/2),m.ar((b-a)/2))
            vals=bank.fields(x) if kind=='mu' else bank.uv_fields(x)
            for l,v in vals.items():
                if kind=='mu':
                    delta=x-4;rs=m.tight(bank.roy.polys[l][1],delta)
                    cross=m.tight(bank.cross[l],delta)/((x-m.ar(m.C))*rs)
                else:
                    rs=m.tight(m.rev(bank.roy.polys[l][1]),x)
                    cross=x*x*m.tight(m.rev(bank.cross[l],l-1),x)/((1+2*m.ar(m.C)*x)*rs)
                ve=upper(v)-q(c['upper_by_spin'][str(l)]['V_upper'])
                xe=upper(cross)-q(c['upper_by_spin'][str(l)]['X4_upper'])
                max_v_excess=max(max_v_excess,ve);max_x_excess=max(max_x_excess,xe)
                assert ve<F(1,10**5) and xe<F(1,10**12),(kind,i,l,float(ve),float(xe))
            vchecks.append(dict(variable=kind,index=i,a=str(a),b=str(b)))
    real=m.VectorReal(ROOT/'output/angles_resummed/angle_000p000.npz');achecks=[];max_a_excess=F(0)
    for i in select(len(acells),97):
        c=acells[i];a,b=F(c['a']),F(c['b']);dens=real.densities(a,b)
        for l,v in dens.items():
            fresh=upper(m.ar(b-a)*v);stored=q(c['mass_upper_by_spin'][str(l)])
            max_a_excess=max(max_a_excess,fresh-stored)
            assert fresh-stored<F(1,10**12),(i,l,float(fresh-stored))
        achecks.append(dict(index=i,a=str(a),b=str(b)))
    here=str(Path(__file__).relative_to(ROOT));hashes[here]=sha(here)
    for p,h in hashes.items():assert sha(p)==h
    out=dict(status='PASS',claim='Exact all-cell partition/merge review, all saved affine inequalities checked; sampled fresh kernel evaluations at340bits, not a second complete vector replay.',
        witness_sha256=m.EXPECTED,precision_bits=m.ctx.prec,real_cells=len(acells),spectral_cells=len(segments),spins=spins,
        comparisons=comparisons,sampled_v_cells=vchecks,sampled_real_cells=achecks,
        sampled_v_max_upper_excess_exact=str(max_v_excess),sampled_x_max_upper_excess_exact=str(max_x_excess),sampled_real_max_upper_excess_exact=str(max_a_excess),
        sha256_manifest=hashes,seconds=time.time()-tic)
    p=ROOT/'output/directional_certification/finite_vectors_independent_audit.json';p.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k not in ['sha256_manifest','sampled_v_cells','sampled_real_cells']},indent=2))

if __name__=='__main__':main()
