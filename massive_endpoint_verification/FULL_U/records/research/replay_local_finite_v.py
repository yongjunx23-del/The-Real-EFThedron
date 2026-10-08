"""Independently replay each stored finite-V local bound with explicit slack."""
import argparse,hashlib,json,sys,time
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'research/tail_certification')]
from certified_finite_v import FiniteV,arb,ctx,ar,EXPECTED

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dyad(q):return F(int(q[0]))*F(2)**int(q[1])

def main():
 ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path)
 ap.add_argument('--slack',default='0.000001');ap.add_argument('--out',type=Path,required=True)
 args=ap.parse_args();path=args.certificate.resolve();certsha=sha(path)
 d=json.loads(path.read_text());assert d['witness_sha256']==EXPECTED
 assert d['spins']==list(range(8,42,2));slack=F(args.slack);assert slack>=0
 deps=[Path(__file__),ROOT/'research/tail_certification/certified_finite_v.py',
  ROOT/'research/tail_certification/certify_finite_v_domain.py',
  ROOT/'research/tail_certification/certify_finite_real.py',
  ROOT/'research/tail_certification/certify_uv.py',ROOT/'research/certify_real_tail.py',
  ROOT/'research/certified_roy_moments.py',ROOT/'output/angles_resummed/angle_000p000.npz',
  ROOT/'output/certification_head/roy_moments.json',
  ROOT/'output/certification_head/strip_12_24_moments.json',
  ROOT/'output/certification_head/strip_24_48_moments.json',
  ROOT/'output/certification_tail/compact_source_cache.pkl']
 hashes={str(p.relative_to(ROOT)):sha(p) for p in deps}
 for kind,lo,hi in [('mu',F(4),F(101)),('w',F(0),F(1,97))]:
  ranges=sorted((F(c['a']),F(c['b'])) for c in d['leaves'] if c['variable']==kind)
  assert ranges[0][0]==lo and ranges[-1][1]==hi
  assert all(a<b for a,b in ranges)
  assert all(x[1]==y[0] for x,y in zip(ranges,ranges[1:]))
 assert all(c['variable'] in ('mu','w') for c in d['leaves'])
 ctx.prec=280;bank=FiniteV();tic=time.time();worst_excess=None;count=0;min_coercivity=None
 for i,c in enumerate(d['leaves']):
  a,b=F(c['a']),F(c['b']);x=arb(ar((a+b)/2),ar((b-a)/2))
  vals=bank.fields(x) if c['variable']=='mu' else bank.uv_fields(x)
  assert set(vals)==set(d['spins']);saved=dyad(c['upper_dyadic']);limit=saved+slack
  assert limit<0
  if min_coercivity is None or -limit<min_coercivity:min_coercivity=-limit
  for l,v in vals.items():
   assert v.is_finite(),(i,c,l,v)
   fresh=dyad(v.upper().man_exp());assert fresh<=limit,(i,l,fresh,limit)
   excess=fresh-saved;count+=1
   if worst_excess is None or excess>worst_excess:worst_excess=excess
  if (i+1)%1000==0:print(json.dumps({'leaves':i+1,'seconds':time.time()-tic}),flush=True)
 assert sha(path)==certsha and all(sha(ROOT/p)==h for p,h in hashes.items())
 out={'status':'PASS','certificate_path':str(path.relative_to(ROOT)),
  'certificate_sha256':certsha,'witness_sha256':EXPECTED,'precision_bits':ctx.prec,
  'kernel_sha256':sha(ROOT/'research/tail_certification/certified_finite_v.py'),
  'producer_sha256':sha(ROOT/'research/tail_certification/certify_finite_v_domain.py'),
  'spins':d['spins'],'cells':len(d['leaves']),'spin_cell_inequalities':count,
  'uniform_slack_exact':str(slack),'every_spin_finite':True,
  'every_spin_upper_le_saved_leaf_upper_plus_slack':True,
  'maximum_fresh_upper_minus_saved_upper_exact':str(worst_excess),
  'minimum_relaxed_coercivity_exact':str(min_coercivity),
  'exact_gapless_domain_partitions':True,
  'claim':'For each stored spectral leaf and every listed spin, V/p is no larger than saved upper_dyadic plus uniform_slack_exact.',
  'seconds':time.time()-tic,'sha256_manifest':hashes}
 args.out.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='sha256_manifest'},indent=2))

if __name__=='__main__':main()
