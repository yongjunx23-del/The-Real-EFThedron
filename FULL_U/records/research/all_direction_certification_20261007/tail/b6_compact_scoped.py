"""Paired analytic compact intervals for b6, using the reviewed Jensen proof."""
from b6_uniform_uv import *
from b6_compact_scoped_kernel import run_scoped

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--precision',type=int,default=320);a=ap.parse_args();ctx.prec=a.precision
 src=ROOT/'output/certification_tail/source_cache_v2_depth10.pkl';cp=ROOT/'output/certification_tail/compact_source_cache.pkl'
 d=pickle.loads(src.read_bytes())['data'];cd=pickle.loads(cp.read_bytes());d['lx']+=F(107,5)
 result=run_scoped(42,d,cd,'.02');low=upper(arb(result['low_energy_upper'])+6);middle=upper(arb(result['middle_energy_upper'])+6)
 assert low<-1 and middle<-1
 prep=json.loads((FOLDER/'prepared_components.json').read_text());manifest=dict(prep['sha256_manifest'])
 for p in[Path(__file__),Path(run_scoped.__code__.co_filename),src,cp]:manifest[str(p.relative_to(ROOT))]=sha(p)
 for p,h in manifest.items():assert sha(ROOT/p)==h,p
 out=dict(status='PASS_B6_COMPACT_SCOPED',base_witness_sha256=sha(W),objective_exact=['1','6'],delta_X4_exact='107/5',minimum_even_spin=42,mu_domains=[['4','5'],['8','49']],first_endpoint_open=True,kappa_exact='1',low_upper_exact=str(low),middle_upper_exact=str(middle),details=result,precision_bits=a.precision,sha256_manifest=manifest)
 (FOLDER/f'compact_scoped_{a.precision}.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in['details','sha256_manifest']},indent=2))

if __name__=='__main__':main()
