"""Independent fixed-leaf finite-V replay with a hash-bound exact-guard record."""
import sys,json,hashlib,time,argparse
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'research/tail_certification')]
from certified_finite_v import FiniteV,arb,ctx,ar,EXPECTED

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(p):return F(int(p[0]))*F(2)**int(p[1])

def main():
 ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path)
 ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 deps=[Path(__file__),ROOT/'research/tail_certification/certified_finite_v.py',
  ROOT/'research/tail_certification/certify_finite_v_domain.py',
  ROOT/'research/tail_certification/certify_finite_real.py',
  ROOT/'research/tail_certification/certify_uv.py',ROOT/'research/certify_real_tail.py',
  ROOT/'research/certified_roy_moments.py',ROOT/'output/angles_resummed/angle_000p000.npz',
  ROOT/'output/certification_head/roy_moments.json',
  ROOT/'output/certification_head/strip_12_24_moments.json',
  ROOT/'output/certification_head/strip_24_48_moments.json',
  ROOT/'output/certification_tail/compact_source_cache.pkl']
 hashes={str(p.relative_to(ROOT)):sha(p) for p in deps};certsha=sha(args.certificate)
 d=json.loads(args.certificate.read_text());assert d['witness_sha256']==EXPECTED
 assert d['spins']==list(range(8,42,2));kappa=F(d['kappa']);assert kappa>0
 for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
  ranges=sorted((F(c['a']),F(c['b'])) for c in d['leaves'] if c['variable']==kind)
  assert ranges[0][0]==lo and ranges[-1][1]==hi
  assert all(a<b for a,b in ranges)
  assert all(x[1]==y[0] for x,y in zip(ranges,ranges[1:]))
 assert all(c['variable'] in ('mu','w') for c in d['leaves'])
 ctx.prec=280;bank=FiniteV();tic=time.time();worst=None;count=0
 for i,c in enumerate(d['leaves']):
  a,b=F(c['a']),F(c['b']);x=arb(ar((a+b)/2),ar((b-a)/2))
  vals=bank.fields(x) if c['variable']=='mu' else bank.uv_fields(x)
  assert set(vals)==set(d['spins'])
  for l,v in vals.items():
   assert v.is_finite() and v.upper()<-ar(kappa),(i,c,l,v)
   bound=dyad(v.upper().man_exp());count+=1
   if worst is None or bound>worst:worst=bound
  if (i+1)%1000==0:print(json.dumps({'leaves':i+1,'seconds':time.time()-tic}),flush=True)
 assert sha(args.certificate)==certsha
 assert all(sha(ROOT/p)==h for p,h in hashes.items())
 out={'status':'PASS','certificate_path':str(args.certificate.resolve().relative_to(ROOT)),
  'certificate_sha256':certsha,'witness_sha256':EXPECTED,'kappa':str(kappa),
  'acceptance_policy':'Every individual spin interval finite and its directed Arb upper strictly below -kappa',
  'precision_bits':ctx.prec,'spins':d['spins'],'leaves':len(d['leaves']),
  'spin_cell_inequalities':count,'exact_gapless_domain_partitions':True,
  'domain':'all mu>4; mu cells [4,101], reciprocal cells w=1/(mu-4) in [0,1/97]',
  'worst_upper_exact':str(worst),'seconds':time.time()-tic,'sha256_manifest':hashes}
 args.out.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='sha256_manifest'},indent=2))

if __name__=='__main__':main()
