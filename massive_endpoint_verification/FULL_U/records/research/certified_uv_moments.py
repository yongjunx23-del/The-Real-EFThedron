"""Exact ball-polynomial contractions for the finite-spin UV source series.

Source denominator coefficients of degree zero and one vanish identically
by s+t+u=4 and c=4/3; the returned series starts with degree two.
"""
from certified_roy_moments import ROOT,ar,load,encode,decode
from fractions import Fraction as F
from flint import arb,ctx
from pathlib import Path
import math,json,hashlib

def build(source,jmax=6,N=32):
    meta=source['metadata'];roy='B' in meta and 'A' not in meta
    kind='roy' if roy else 'strip';A=F(4) if roy else F(meta['A']);B=F(meta['B'])
    s_center=F(4) if roy else F(meta['s_center']);u_center=F(4) if roy else F(meta['u_center_minus'])
    tmax=(B-4)/2 if roy else F(4,3);norm=F(meta['absolute_W_mass_bound']);c=F(4,3)
    T=source['T'];ms=source['M_s'];mu=source['M_u']
    assert len(T)>jmax+N+1 and len(mu[0])>=N+2 and len(ms[0])>=N+2
    def shift(row,k,constant,sign=1):
        return sum((ar(math.comb(k,r)*constant**(k-r)*sign**r)*row[r] for r in range(k+1)),arb(0))
    def crossed(j,k):return sum((ar(math.comb(k,r)*(2*c)**(k-r)*(-1)**r)*T[j+r] for r in range(k+1)),arb(0))
    aa=[]
    for j in range(jmax+1):
        aa.append([shift(ms[j],k+2,s_center)+shift(mu[j],k+2,-u_center)-ar(c**(k+2))*T[j]-crossed(j,k+2) for k in range(N)])
    bb=[T[k+2]+crossed(0,k+2)-2*ar(c**(k+2))*T[0] for k in range(N)]
    mm={'kind':kind,'A':str(A),'B':str(B),'Tmax':str(tmax),'source_W_norm_bound':str(norm),
        'jmax':jmax,'N':N,'witness_sha256':meta['witness_sha256'],'precision_bits':ctx.prec,
        'A_definition':'A[j][k]=integral W*t^j*(s^(k+2)+u^(k+2)-c^(k+2)-(2c-t)^(k+2))',
        'B_definition':'B[k]=integral W*(t^(k+2)+(2c-t)^(k+2)-2c^(k+2))',
        'UV_remainder':'At x>=0 and B*x<1, |error A_j(x)| <=4*norm*Tmax^j*B^2*(B*x)^N/(1-B*x); |error B(x)| <=4*norm*B^2*(B*x)^N/(1-B*x).',
        'assumptions':'All |s|,|u|,|t|,c,|2c-t|<=B on the source; checked for these windows.',
        'scaled_V':'-x/16*(sum_j C_lj*(x/(1-4*x))^j*A_j(x)+P_l(Z)*B(x)), x=1/mu; V scaled by (1-beta^2)^-2'}
    assert max(tmax,2*c+tmax,c)<B and A>=4
    return {'metadata':mm,'A':aa,'B':bb}

def save(d,path):
    out={'metadata':d['metadata'],'A':[[encode(x) for x in row] for row in d['A']],'B':[encode(x) for x in d['B']]}
    Path(path).write_text(json.dumps(out,separators=(',',':')))

def load_uv(path):
    d=json.loads(Path(path).read_text())
    return {'metadata':d['metadata'],'A':[[decode(x) for x in row] for row in d['A']],'B':[decode(x) for x in d['B']]}

if __name__=='__main__':
    ctx.prec=260;base=ROOT/'output/certification_head/N64'
    for name in ('roy','strip_12_24','strip_24_48'):
        src=base/(name+'_moments.json');d=build(load(src));d['metadata']['source_table_sha256']=hashlib.sha256(src.read_bytes()).hexdigest()
        out=ROOT/'output/certification_head'/(name+'_uv_moments.json');save(d,out)
        print(name,'A00',d['A'][0][0],'B0',d['B'][0])
