"""Source moment and norm data for compact-energy tail inequalities."""
import sys,re,pickle,math,hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import sympy as sy
from certify_uv import bernstein,integral,split_box,ROOT

def builder_hash():
 return hashlib.sha256(Path(__file__).read_bytes()+(Path(__file__).parent/'certify_uv.py').read_bytes()+(ROOT/'research/certify_roy_dangerous_sign.py').read_bytes()).hexdigest()

def compute(path,depth=10,N=40):
 st=np.load(path);x,y=sy.symbols('x y');roy=0;strips={};B=F(119,10);Ds={}
 for name,val in zip(st['names'],st['multipliers']):
  v=F(float(val));m=re.fullmatch(r'Roy_B([^_]+)_J(\d+)_n(\d+)',str(name))
  if m:roy+=2*sy.Rational(v)*sy.legendre(int(m[3]),2*x-1)*sy.legendre(int(m[2]),y);continue
  m=re.fullmatch(r'Taper(\d+)_Strip_A([^_]+)_B([^_]+)_T([^_]+)_j(\d+)_n(\d+)',str(name))
  if m:
   power,A,E,j,n=int(m[1]),F(m[2]),F(m[3]),int(m[5]),int(m[6]);assert power==1
   key=(A,E);strips[key]=strips.get(key,0)+sy.Rational(v)*(1-(2*x-1)**2)*sy.legendre(n,2*x-1)*sy.legendre(j,2*y-1)
   Ds[key]=Ds.get(key,F(0))+abs(v)*(4+n*(n+1))
 out={'roy_pos':F(0),'roy_neg':F(0),'roy_abs':F(0),'strip_abs':F(0),'strip_pos':F(0),'strip_neg':F(0),'strip_D':Ds,'moments':[F(0)]*N}
 for group,W in [('roy',roy)]+[('strip',v) for v in strips.values()]:
  bs=bernstein(W,x,y)
  for box,a,b,d,e in split_box(bs,depth):
   scale=(b-a)*(e-d)/(len(box)*len(box[0]))
   out[group+'_pos']+=scale*sum((sum((max(q,F(0)) for q in row),F(0)) for row in box),F(0))
   out[group+'_neg']+=scale*sum((sum((max(-q,F(0)) for q in row),F(0)) for row in box),F(0))
  out[group+'_abs']=out[group+'_pos']+out[group+'_neg']
  p=sy.Poly(sy.expand(W),x,y);terms=[(i,j,F(v)) for (i,j),v in p.terms()]
  if group=='roy':
   hs=[]
   for r in range(2*N+1):
    hs.append(sum((v/F(i+r+1)*F(math.factorial(j)*math.factorial(r),math.factorial(j+r+1)) for i,j,v in terms),F(0)))
   for k in range(N):
    n=2*k+2
    out['moments'][k]+=2*sum((F(math.comb(n,r))*F(4,3)**(n-r)*((B-4)/2)**r*hs[r] for r in range(n+1)),F(0))
  else:
   hs=[sum((v/F(i+1)/F(j+r+1) for i,j,v in terms),F(0)) for r in range(2*N+1)]
   for k in range(N):
    n=2*k+2
    out['moments'][k]+=2*F(4,3)**n*sum((F(math.comb(n,r))*(-1)**r*2**(n-r)*hs[r] for r in range(n+1)),F(0))
 out['witness_sha256']=hashlib.sha256(Path(path).read_bytes()).hexdigest()
 out['builder_sha256']=builder_hash()
 return out
if __name__=='__main__':
 p=ROOT/'output/angles_resummed/angle_000p000.npz';d=compute(p);f=ROOT/'output/certification_tail/compact_source_cache.pkl';f.write_bytes(pickle.dumps(d))
 print({k:float(v) for k,v in d.items() if isinstance(v,F)});print('D', {str(k):float(v) for k,v in d['strip_D'].items()});print('moments',[float(v) for v in d['moments'][:4]])
