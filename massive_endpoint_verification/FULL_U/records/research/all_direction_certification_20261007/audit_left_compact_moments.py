"""Independent exact beta-integral T_j checks and shared-table overlap audit."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import sys,json,re,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'tmp/certification_deps'),str(ROOT/'research')]
from flint import fmpq_poly


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def interval(q):return [F(int(m))*F(2)**int(e) for m,e in q]


def moment_power_legendre(j,n):
    return F(0) if j<n else F(factorial(j)**2,factorial(j-n)*factorial(j+n+1))


def main():
    folder=ROOT/'output/competitive_20261007/left/compact/moments'
    old=ROOT/'output/all_direction_certification_20261007/head/left/moments128'
    if not old.exists():
        choices=list((ROOT/'output/all_direction_certification_20261007/head/g2_robust_stage1').rglob('roy_moments.json'))
        assert len(choices)==1;old=choices[0].parent
    witness=ROOT/'output/full_direction_20261007/left/angle180_final.npz'
    st=np.load(witness);rows=[(str(n),F(float(v))) for n,v in zip(st['names'],st['multipliers'])]
    X=fmpq_poly([0,1]);P=[fmpq_poly([1]),X]
    for n in range(1,6):P.append(((2*n+1)*X*P[-1]-n*P[-2])/(n+1))
    records=[];hashes={str(witness.relative_to(ROOT)):sha(witness),str(Path(__file__).relative_to(ROOT)):sha(__file__)}
    for label in ['roy','strip_12_24','strip_24_48']:
        path=folder/f'{label}_moments.json';q=json.loads(path.read_text());assert q['metadata']['witness_sha256']==sha(witness)
        hashes[str(path.relative_to(ROOT))]=sha(path)
        for j,ball in enumerate(q['T']):
            exact=F(0)
            for name,lam in rows:
                if label=='roy':
                    m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',name)
                    if not m:continue
                    B,J,n=F(m[1]),int(m[2]),int(m[3])
                    angular=sum((F(str(co))*F(factorial(b)*factorial(j),factorial(b+j+1)) for b,co in enumerate(P[J].coeffs())),F(0))
                    exact+=2*lam*(-(B-4)/2)**j*moment_power_legendre(j,n)*angular
                else:
                    m=re.fullmatch(r'Taper1_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',name)
                    if not m or label!=f'strip_{m[1]}_{m[2]}':continue
                    J,n=int(m[4]),int(m[5]);energy=F(2,3) if n==0 else -F(2,15) if n==2 else F(0)
                    exact+=lam*(-F(4,3))**j*energy*(-1)**J*moment_power_legendre(j,J)
            lo,hi=interval(ball);assert lo<=exact<=hi,(label,j)
        before=old/f'{label}_moments.json';qq=json.loads(before.read_text());assert qq['metadata']['witness_sha256']==sha(witness)
        overlaps=0
        for key in ['T','M_t','M_u','M_s']:
            left=q[key];right=qq[key]
            pairs=list(zip(left,right)) if key in ['T','M_t'] else [(x,y) for a,b in zip(left,right) for x,y in zip(a,b)]
            for x,y in pairs:
                a,b=interval(x);c,d=interval(y);assert max(a,c)<=min(b,d),(label,key,overlaps)
                overlaps+=1
        hashes[str(before.relative_to(ROOT))]=sha(before)
        records.append(dict(label=label,exact_beta_integral_T_checks=len(q['T']),shared_moment_overlap_checks=overlaps,
            new_order=q['metadata']['order'],earlier_order=qq['metadata']['order']))
    out=dict(status='PASS_INDEPENDENT_EXACT_BETA_MOMENTS_AND_SHARED_TABLES',witness_sha256=sha(witness),records=records,sha256_manifest=hashes,
        scope='All T_j are checked against exact factorial/beta integrals. Shared moments are compared to the earlier table constructed using a different validated Gauss rule. The new table still relies on its asserted polynomial degree and validated Gauss-node/weight enclosures for remaining entries.')
    (folder/'independent_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='sha256_manifest'},indent=2))


if __name__=='__main__':main()

