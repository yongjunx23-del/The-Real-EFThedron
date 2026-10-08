"""Replay the original left head and enclose ten exact objective shifts.

The witness is defined by base NPZ SHA plus objective(-1,-b); lambda and
source tables are unchanged. Original outputs remain independent of the
additional family records. No source/kernel quadrature is rebuilt.
"""
from finite_gap import *
from head_certification_spin0 import support_upper

LEFT_SHA='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
# Frozen ten-slope variant: only this list and descriptive count differ.
SLOPES=[F(1,4),F(1,2),F(1),F(2),F(4),F(8),F(16),F(32),F(64),F(128)]

def pscaled(m,beta):
    """Exact p/e^2, p=P_l(1+2C/(mu-4))/(mu-C)^3; e=1-beta^2."""
    z=beta*beta;e=1-z;dbar=4*z;Qbar=4-exact(C)*e
    R=sum((exact(c)*dbar**(m.ell-k)*(exact(C)*e)**k for k,c in enumerate(m.coeff)),z*0)
    return e*R/(dbar**m.ell*Qbar**3)

class CellCache:
    """Cache pure field enclosures while the original trapezoid runs.

    The source activity/PV mode is selected inside the original model.
    We never mutate its direction and never reuse an endpoint *support*
    value for a shifted direction. Endpoint field balls are reusable.
    """
    def __init__(self,m):
        self.m=m;self.points={};self.clear()
        self.old_fields=m.fields;self.old_series=m.fields_series;self.old_norm=m.normalized_interval
        m.fields=self.fields;m.fields_series=self.series;m.normalized_interval=self.norm
    def clear(self):self.intervals={};self.jets={};self.norms={}
    def fields(self,a,b,omit_sources=()):
        key=(a,b,tuple(omit_sources));cache=self.points if a==b and not omit_sources else self.intervals
        if key not in cache:cache[key]=self.old_fields(a,b,omit_sources)
        return cache[key]
    def series(self,a,b):
        if (a,b) not in self.jets:self.jets[a,b]=self.old_series(a,b)
        return self.jets[a,b]
    def norm(self,a,b):
        if (a,b) not in self.norms:self.norms[a,b]=self.old_norm(a,b)
        return self.norms[a,b]

def shifted_bounds(m,cache,a,b,original_upper):
    """Each returned bound encloses the complete cell at its exact shift."""
    width=exact(b-a);ball=arb(exact((a+b)/2),exact((b-a)/2))
    P=pscaled(m,ball);u,v,_=cache.fields(a,b)
    ends=[]
    for endpoint in (a,b):
        eu,ev,_=cache.fields(endpoint,endpoint)
        ends.append((eu,ev,pscaled(m,exact(endpoint))))
    try:
        uj,vj,err=cache.series(a,b)
        pj=pscaled(m,arb_series([ball,1],prec=3))
        jet_ok=all(x.is_finite() for x in [uj[2],vj[2],pj[2],err])
    except (ValueError,ZeroDivisionError):jet_ok=False
    # These are exactly the original applicability conditions. In particular,
    # ell0 must NOT use the old normalized real-source branch above beta=.5.
    norm=None
    if b<F(1,2) or (m.ell and F(1,2)<=b<F(99,100)):
        try:norm=cache.norm(a,b)
        except (ValueError,ZeroDivisionError):pass
    out=[]
    for slope in SLOPES:
        shift=exact(slope)
        # Since p>0 for mu>4, H(U,V-b*p)<=H(U,V). The original
        # certificate is therefore always a valid fallback for every b>=0.
        options=[original_upper]
        if out:options.append(out[-1])  # H also decreases between sorted slopes.
        riemann=m.factor*width*support_upper(u,v-shift*P)
        if riemann.is_finite():options.append(riemann)
        if jet_ok:
            vshift=vj-shift*pj
            rem=width**3/8*(2*uj[2].abs_upper()+4*vshift[2].abs_upper())+4*width*err
            endsupport=sum((support_upper(eu,ev-shift*ep) for eu,ev,ep in ends),arb(0))
            trap=m.factor*(width/2*endsupport+rem)
            if trap.is_finite():options.append(trap)
        if norm is not None:
            nu,nv,pinv,e=norm;vn=nv-shift
            if vn.upper()<0:
                small=m.factor*width*nu.abs_upper()**2*pinv.abs_upper()*e.abs_upper()**2/(2*(-vn.upper()))
                if small.is_finite():options.append(small)
        out.append(min(options,key=float).upper())
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--moments',type=Path,required=True)
    ap.add_argument('--record',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--family-out',type=Path,required=True);ap.add_argument('--ceiling',required=True)
    ap.add_argument('--bits',type=int,default=256);a=ap.parse_args();ctx.prec=a.bits;tic=time.time()
    path=ROOT/'output/full_direction_20261007/left/angle180_final.npz';r=json.loads(a.record.read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest()==r['witness_sha256']==LEFT_SHA
    for name,digest in r['dependencies'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    assert hashlib.sha256((a.moments/'manifest.json').read_bytes()).hexdigest()==r['source_manifest_sha256']
    klass=FiniteGap if any('finite_gap.py' in x for x in r['dependencies']) else FiniteQuartic if any('finite_quartic.py' in x for x in r['dependencies']) else FiniteCentered
    m=klass(path,r['spin'],a.moments,N=r.get('N',128));assert m.base.direction==[F(-1),F(0)]
    gaps=[tuple(map(F,x)) for x in r['gap_beta_intervals']];assert gaps==m.gaps(30)
    cells=[(F(c['a']),F(c['b'])) for c in r['cells_data']]
    partition=sorted(gaps+cells);assert partition[0][0]==0 and partition[-1][1]==1
    assert all(lo<hi for lo,hi in partition) and all(x[1]==y[0] for x,y in zip(partition,partition[1:]))
    gapmass=sum((m.gap_bound(lo,hi) for lo,hi in gaps),arb(0))
    total=gapmass;family=[gapmass for _ in SLOPES];records=[];cache=CellCache(m)
    for k,(lo,hi) in enumerate(cells):
        cache.clear();up,_=m.trapezoid(lo,hi);assert up.is_finite()
        vals=shifted_bounds(m,cache,lo,hi,up);total+=up
        for i,value in enumerate(vals):family[i]+=value
        records.append(dict(a=str(lo),b=str(hi),original_upper_dyadic=list(map(str,up.upper().man_exp())),
            shifted_upper_dyadic=[list(map(str,v.man_exp())) for v in vals]))
        if (k+1)%500==0:print('REPLAYED_FAMILY',m.ell,k+1,'of',len(cells),'seconds',time.time()-tic,flush=True)
    passed=bool(total<arb(a.ceiling));assert passed,(str(total),a.ceiling)
    deps={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('finite220.py'),Path(__file__).with_name('forward220.py')]}
    original=dict(witness_sha256=m.sha,source_record=str(a.record),source_record_sha256=hashlib.sha256(a.record.read_bytes()).hexdigest(),
        model_class=klass.__name__,spin=m.ell,cells_replayed=len(cells),gaps_replayed=len(gaps),partition_covers_beta_zero_to_one_exactly=True,
        recomputed_upper=str(total.upper()),upper_dyadic=list(map(str,total.upper().man_exp())),ceiling=a.ceiling,passed=True,
        elapsed_seconds=time.time()-tic,scope='Original left finite single-spin all-energy component; family results are separate.',dependency_sha256=deps)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(original,indent=2))
    familyout=dict(status='PASS_FIXED_PARTITION_FAMILY_HEAD_COMPONENT',base_witness_sha256=m.sha,spin=m.ell,
        witness_definition='Same exact base multipliers and source rows, objective(-1,-b); V_b=V_base-b*p. No different-witness source hash is asserted.',
        slopes=[str(b) for b in SLOPES],objectives=[['-1',str(-b)] for b in SLOPES],
        upper_dyadic=[list(map(str,v.upper().man_exp())) for v in family],upper_float=[float(v.upper()) for v in family],
        gap_bound_reused_by_monotonicity=str(gapmass.upper()),gap_bound_dyadic=list(map(str,gapmass.upper().man_exp())),
        gap_intervals=r['gap_beta_intervals'],cells=records,original_replay_sha256=hashlib.sha256(a.out.read_bytes()).hexdigest(),
        source_record_sha256=hashlib.sha256(a.record.read_bytes()).hexdigest(),source_manifest_sha256=r['source_manifest_sha256'],
        dependencies={**r['dependencies'],**deps},precision_bits=a.bits,seconds=time.time()-tic,
        scope='One finite spin at all energies for ten exact objective shifts. Complete directions require every finite head and a correctly shifted same-base all-spin tail.')
    a.family_out.parent.mkdir(parents=True,exist_ok=True);a.family_out.write_text(json.dumps(familyout,indent=2))
    print(json.dumps({k:v for k,v in familyout.items() if k not in ['cells','dependencies','upper_dyadic']},indent=2),flush=True)
if __name__=='__main__':main()
