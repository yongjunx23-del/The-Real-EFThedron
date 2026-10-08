from pathlib import Path
import sys,math,json,time
import numpy as np
from forward220 import Forward220,F,exact,ROOT
from flint import ctx
from angle_audit import normalized_inventory
ctx.prec=224
path=ROOT/'output/full_direction_20261007/right/g2_guarded192_refined.npz'
st=np.load(path);keep=[i for i,n in enumerate(st['names']) if str(n).startswith(('Dt','X5','X6','X7','X8'))]
names=st['names'][keep];lam=st['multipliers'][keep]
energies=np.array([4.0001,4.02,4.5,5.,10.,20.,100.,1000.,1e6])
tic=time.time();records=[]
for ell in (0,2,4,14,40):
    m=Forward220(path,ell);u0,v0,p,u,v=normalized_inventory(names,energies,[ell],96,96)
    uu=np.einsum('i,ijk->jk',lam,u)[0];vv=np.einsum('i,ijk->jk',lam,v)[0]
    for i,mu in enumerate(energies):
        beta=F(float(math.sqrt(1-4/mu)));b=exact(beta);e=1-b*b
        ub,vs=m.derivative(b);vs+=m.higher(b)-exact(m.direction[1])*e*(sum(exact(c)*(4*b*b)**(ell-k)*(exact(F(4,3))*e)**k for k,c in enumerate(m.coeff)))/((4*b*b)**ell*(4-exact(F(4,3))*e)**3)
        unum=float(ub/b*e*e);vnum=float(vs*e*e)/p[0,i]
        ur=abs(unum-uu[i])/(1+abs(uu[i]));vr=abs(vnum-vv[i])/(1+abs(vv[i]))
        records.append(dict(ell=ell,mu=float(mu),U_rel=ur,Vn_rel=vr))
        assert ur<1e-7 and vr<1e-7,(records[-1],unum,uu[i],vnum,vv[i])
    print('VALID',ell,max(r['U_rel'] for r in records if r['ell']==ell),max(r['Vn_rel'] for r in records if r['ell']==ell),flush=True)
out=ROOT/'output/all_direction_certification_20261007/head/forward_validation.json'
out.write_text(json.dumps(dict(records=records,seconds=time.time()-tic,scope='Independent float comparison is a code diagnostic, not a proof'),indent=2))

