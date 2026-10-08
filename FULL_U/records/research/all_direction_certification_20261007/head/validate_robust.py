from finite_centered import *
from angle_audit import fields
ctx.prec=256;tic=time.time();path=ROOT/'output/competitive_20261007/robust/g2_robust_stage1.npz';st=np.load(path)
mom=ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1/moments128'
energies=np.array([4.0001,4.02,4.5,5.,7.,10.,11.8,11.99,12.1,18.,23.9,24.1,36.,47.9,48.1,80.,127.,128.,160.,1000.,1e6]);rows=[]
for ell in (0,2,4,8,14,40):
    m=FiniteCentered(path,ell,mom,N=128)
    uu,vv,pp=fields(st['names'],st['multipliers'],st['direction'],energies,[ell],kernel_order=128,angular_order=128)
    rec=[]
    for i,s in enumerate(energies):
        beta=F(float(math.sqrt(1-4/s)));u,v,_=m.fields(beta,beta);e=1-float(beta)**2
        un=float(u)*e**2;vn=float(v)*e**2/pp[0,i]
        rec.append(dict(mu=float(s),U_relative=abs(un-uu[0,i])/(1+abs(uu[0,i])),Vnorm_relative=abs(vn-vv[0,i])/(1+abs(vv[0,i]))))
    rows.append(dict(spin=ell,points=rec));print(ell,max(q['U_relative'] for q in rec),max(q['Vnorm_relative'] for q in rec),flush=True)
out=ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1/independent_field_comparison.json'
out.write_text(json.dumps(dict(witness_sha256=m.sha,rows=rows,seconds=time.time()-tic,scope='Independent numerical comparison diagnostic, not itself a rigorous integral certificate'),indent=2))
