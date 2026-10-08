"""Rebuild and assemble the FIXED upper-witness all-spin tail certificate.

This does not certify the finite-spin head. Defaults to a full replay. Use
--assemble-only only to check and assemble the already generated components.
"""
import sys,json,hashlib,argparse,subprocess
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tmp/certification_deps'))
from flint import arb,ctx
ctx.prec=160
HERE=Path(__file__).resolve().parent
OUT=ROOT/'output/certification_tail'
WIT=ROOT/'output/angles_resummed/angle_000p000.npz'
EXPECTED='bbbc5990b459e2b3adf276f53fd6105d86fe2274b7b1b2313195c55bb376511c'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def command(args):
 result=subprocess.run([sys.executable,*map(str,args)],cwd=ROOT,text=True,capture_output=True)
 if result.returncode:
  print(result.stdout);print(result.stderr);raise RuntimeError('Component replay failed: '+str(args))
 print('PASS',Path(args[0]).name,flush=True)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--assemble-only',action='store_true');ap.add_argument('--rebuild',action='store_true',help='Full replay (also default)');a=ap.parse_args()
 assert sha(WIT)==EXPECTED,'This verifier is restricted to the fixed resummed witness.'
 if not a.assemble_only:
  command([HERE/'certify_expanded_roy_sign.py',WIT,'--output',OUT/'expanded_roy_sign.json'])
  command([HERE/'certify_uv.py','--witness',WIT,'--depth','10','--L','42','--M','49','--dynamic','--step','.004','--output',OUT/'uv_L42_M49_v2.json'])
  command([HERE/'compact_source_data.py'])
  command([HERE/'certify_compact.py','--L','42','--output',OUT/'compact_L42.json'])
  command([ROOT/'research/certify_real_tail.py',WIT,'--L','40','--target','2e-5','--max-cells','8000','--combined-strips','--out',OUT/'real_A40.json'])
 sign=json.loads((OUT/'expanded_roy_sign.json').read_text());uv=json.loads((OUT/'uv_L42_M49_v2.json').read_text());co=json.loads((OUT/'compact_L42.json').read_text());aa=json.loads((OUT/'real_A40.json').read_text())
 assert all(x['witness_sha256']==EXPECTED for x in [sign,uv,co]) and aa['sha256']==EXPECTED
 assert uv['L']==co['L']==42 and uv['M']==49 and uv['dynamic_mu_constraint']
 assert aa['L']==40 and aa['combined_strip_integration'] and aa['target_met']
 assert sign['W_certified_lower_bound']>0
 kappa=Fraction(31,50)
 for value in [uv['small_tau_max_upper'],uv['large_tau_upper'],co['low_energy_upper'],co['middle_energy_upper']]:
  assert arb(value)<-arb(31)/50, 'Insufficient coercivity margin'
 mantissa,exponent=aa['combined_A_upper_dyadic'];A=Fraction(int(mantissa))*Fraction(2)**int(exponent)
 assert A<Fraction(1,50000)
 tail=A/(2*kappa)
 files=[WIT,HERE/'certify_uv.py',HERE/'certify_compact.py',HERE/'compact_source_data.py',HERE/'certify_expanded_roy_sign.py',HERE/'verify_fixed_tail.py',HERE/'PROOF.txt',ROOT/'research/certify_roy_dangerous_sign.py',ROOT/'research/certify_real_tail.py',OUT/'source_cache_v2_depth10.pkl',OUT/'compact_source_cache.pkl',OUT/'expanded_roy_sign.json',OUT/'uv_L42_M49_v2.json',OUT/'compact_L42.json',OUT/'real_A40.json']
 manifest={str(p.relative_to(ROOT)):sha(p) for p in files}
 result={'witness_sha256':EXPECTED,'physical_domain':'all mu>4 and all even l>=42',
  'V_over_p_upper_rational':'-31/50','A40_upper_exact':str(A),
  'A40_upper_outward_rounded':'0.000020000000000000000000',
  'disk_support_tail_upper_exact':str(tail),'disk_support_tail_upper_outward_rounded':'0.00001612904',
  'formula':'sqrt(U^2+V^2)+V <= U^2/(2*kappa*p), integrate and sum all even l>=42',
  'external_g2_cap_used':False,'replayed_components_this_run':not a.assemble_only,
  'status':'All-spin/all-energy omitted-spin tail certified under the stated amplitude assumptions.',
  'limitation':'The finite-spin head l=0,2,...,40 still requires an independent certified integration; this tail result alone is not a complete exclusion bound.',
  'sha256_manifest':manifest}
 (OUT/'fixed_tail_certificate.json').write_text(json.dumps(result,indent=2))
 print(json.dumps({k:v for k,v in result.items() if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
