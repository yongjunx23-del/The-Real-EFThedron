"""Scoped old150 b=4,6 extension. This does not certify either full objective."""
from pathlib import Path
from fractions import Fraction as F
import sys,pickle,json,hashlib,time,argparse
ROOT=Path(__file__).resolve().parents[3]
sys.dont_write_bytecode=True
sys.path[:0]=[str(Path(__file__).parent),str(ROOT/'research/tail_certification')]
from positive_uv_jensen import uv_bounds,up,ar,arb,ctx
from compact_jensen import run
OUT=ROOT/'output/all_direction_certification_20261007/tail/old150_extension'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=280);args=ap.parse_args()
 ctx.prec=args.precision;tic=time.time();OUT.mkdir(parents=True,exist_ok=True)
 fixed_path=ROOT/'output/certification_tail/fixed_tail_certificate.json';fixed=json.loads(fixed_path.read_text())
 assert fixed['replayed_components_this_run'] and F(fixed['V_over_p_upper_rational'])==-F(31,50)
 hashes=dict(fixed['sha256_manifest']);hashes[str(fixed_path.relative_to(ROOT))]=sha(fixed_path)
 for p,h in hashes.items():assert sha(ROOT/p)==h,p
 src=ROOT/'output/certification_tail/source_cache_v2_depth10.pkl';cp=ROOT/'output/certification_tail/compact_source_cache.pkl'
 saved=pickle.loads(src.read_bytes());d=saved['data'];cd=pickle.loads(cp.read_bytes())
 assert saved['sha']==cd['witness_sha256']==fixed['witness_sha256']
 cases=[(F(4),F(72,5)),(F(6),F(108,5))]
 assert all(d['lx']+dd<0 for b,dd in cases)
 big,cells=uv_bounds(d,cases,42,49,F(1,250));results=[]
 A=F(fixed['A40_upper_exact'])
 for b,dd in cases:
  trial=dict(d);trial['lx']+=dd;L=46 if b==4 else 108
  compact=run(L,trial,cd,'.02')
  cu=max(up(arb(compact['low_energy_upper'])+ar(b)),up(arb(compact['middle_energy_upper'])+ar(b)))
  assert cu<0
  choices=[]
  for eta in map(F,['1/1000','1/500','1/200','1/100','1/50','1/20','1/10','1/5','2/5','1/2','1','2']):
   assert big[b]+eta<0
   pos=sum((2*max(F(0),F(c['V_over_p_upper_by_slope'][str(b)])+eta)*F(c['p_mass_upper_exact']) for c in cells),F(0))
   val=A/(2*eta)+pos
   choices.append(dict(eta_exact=str(eta),positive_V_term_exact=str(pos),real_term_exact=str(A/(2*eta)),UV_support_upper_exact=str(val)))
  best=min(choices,key=lambda r:F(r['UV_support_upper_exact']))
  results.append(dict(b=str(b),delta_X4=str(dd),status='SCOPED_COMPONENTS_PASS',
   compact_spin_min=L,compact_mu_domain=['4','49'],compact_V_upper_exact=str(cu),compact=compact,
   UV_spin_min=42,UV_mu_domain=['49','infinity'],large_tau_upper_exact=str(big[b]),UV_support_best=best,
   uncovered_compact_spins=list(range(42,L,2)),UV_trials=choices))
 for p in [Path(__file__),Path(__file__).parent/'positive_uv_jensen.py',Path(__file__).parent/'compact_jensen.py']:
  hashes[str(p.relative_to(ROOT))]=sha(p)
 out=dict(status='SCOPED_COMPONENTS_PASS',witness_sha256=fixed['witness_sha256'],precision_bits=ctx.prec,
  results=results,uv_cells=cells,real_norm_A_all_spins42_up_exact=str(A),seconds=time.time()-tic,
  claim='UV support bound for all even ell>=42, mu>=49; compact V<0 for ell>=46 (b4) or108 (b6). This is not a complete directional coupling certificate.',
  missing='Finite compact bridge spins, finite heads and full assembly/replay remain required.',sha256_manifest=hashes)
 path=OUT/f'partial_components_{ctx.prec}.json';path.write_text(json.dumps(out,indent=2))
 print(json.dumps(dict(status=out['status'],seconds=out['seconds'],results=[dict(b=r['b'],compact_spin_min=r['compact_spin_min'],compact_upper=float(F(r['compact_V_upper_exact'])),UV_support=float(F(r['UV_support_best']['UV_support_upper_exact']))) for r in results]),indent=2))
if __name__=='__main__':main()
