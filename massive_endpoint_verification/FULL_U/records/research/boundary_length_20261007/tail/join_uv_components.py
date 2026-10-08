"""Join two already-frozen UV component packets without rewriting either."""
from assemble_uv import ROOT,CONFIG,sha,F
from pathlib import Path
import json,argparse
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--angle',required=True,choices=list(CONFIG));a=ap.parse_args();cfg=CONFIG[a.angle];folder=ROOT/f'output/boundary_length_20261007/tail/angle{a.angle}';manifest={};seen=set()
 def ingest(p,expected=None):
  p=Path(p).resolve();key=str(p.relative_to(ROOT));h=sha(p)
  if expected is not None:assert h==expected,key
  if key in manifest:assert manifest[key]==h,key
  manifest[key]=h
  if p.suffix!='.json':return None
  d=json.loads(p.read_text())
  if key not in seen:
   seen.add(key)
   for q,hh in d.get('sha256_manifest',{}).items():ingest(ROOT/q,hh)
  return d
 source=ingest(folder/'source_exact_rebuild.json')
 assert source['status']=='PASS_EXACT_RATIONAL_SOURCE_REBUILD'and source['witness_sha256']==cfg['sha']and source['objective_exact']==cfg['objective']and source['moment_count']==972 and source['source_remainder_coefficients']==9 and source['all_entries_equal']is True
 fp=folder/'finite_uv_complete.json';ip=folder/'infinite_uv_complete.json';f=ingest(fp);i=ingest(ip)
 for d in[f,i]:assert d['witness_sha256']==cfg['sha']and d['objective_exact']==cfg['objective']and d['mu_domain']==['49','infinity']
 assert f['status']=='PASS_COMPLETE_DIRECTIONAL_FINITE_UV'and f['spins']==list(range(42,450,2))and f['first_even_spin']==42 and f['last_even_spin']==448
 assert i['status']=='PASS_COMPLETE_DIRECTIONAL_INFINITE_UV'and i['lmin']==450
 assert F(f['minimum_kappa_exact'])>0 and F(i['kappa_exact'])==F(cfg['kappa'])
 ingest(Path(__file__))
 d=dict(status='PASS_COMPLETE_DIRECTIONAL_UV_COERCIVITY',witness_sha256=cfg['sha'],objective_exact=cfg['objective'],finite_first_even_spin=42,finite_last_even_spin=448,infinite_first_even_spin=450,finite_kappa_exact=f['minimum_kappa_exact'],infinite_kappa_exact=i['kappa_exact'],energy_domain=['49','infinity'],finite_component=str(fp.relative_to(ROOT)),infinite_component=str(ip.relative_to(ROOT)),finite_local_kappa_rule=f['local_bound_rule'],sha256_manifest=manifest,scope='Every even spin>=42, all mu>=49. Separate compact coercivity and same-witness real-norm bounds are required for a complete physical support bound.')
 p=folder/'uv_complete.json';assert not p.exists(),'Never overwrite a frozen combined packet'
 p.write_text(json.dumps(d,indent=2));print(json.dumps({k:v for k,v in d.items()if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
