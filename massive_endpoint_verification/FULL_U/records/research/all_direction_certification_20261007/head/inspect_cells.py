from finite_taylor import *
ctx.prec=256;m=FiniteTaylor(ROOT/'output/full_direction_20261007/right/g2_guarded192_refined.npz',2,ROOT/'output/all_direction_certification_20261007/head/g2/moments128',N=128)
gaps=m.gaps();points=sorted(set([F(i,256) for i in range(257)]+[F(math.ceil(math.sqrt(1-4/float(m.uv_start))*2**20),2**20)]+[x for g in gaps for x in g]));out=[]
for a,b in zip(points,points[1:]):
    if any(a>=ga and b<=gb for ga,gb in gaps):continue
    u,l=m.trapezoid(a,b)
    if not u.is_finite():continue
    out.append(dict(a=str(a),b=str(b),upper=float(u),midpoint=float(l)))
print(json.dumps(sorted(out,key=lambda q:q['upper'],reverse=True)[:15],indent=2),flush=True)
