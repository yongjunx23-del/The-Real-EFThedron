"""Extract the stronger exact large-tau margin from two frozen enclosures."""
from left_large_tau import ROOT,OUT,WITNESS
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 manifest={};ups=[];files=[]
 for bits in[400,480]:
  p=OUT/f'large_tau_L108_{bits}.json';d=json.loads(p.read_text());files.append(str(p.relative_to(ROOT)))
  assert d['status']=='PASS_UNBOUNDED_TAU'and d['witness_sha256']==sha(WITNESS)and d['objective_exact']==['-1','0']
  assert d['lmin']==108 and d['mu_domain']==['49','infinity']and d['tau_domain']==['256','infinity']and F(d['conservative_objective_offset_exact'])==1
  u=F(d['V_over_p_upper_exact']);assert u<-3000;ups.append(u)
  for f,h in d['sha256_manifest'].items():assert sha(ROOT/f)==h;manifest[f]=h
  manifest[str(p.relative_to(ROOT))]=sha(p)
 manifest[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
 out=dict(status='PASS_LEFT_LARGE_TAU',witness_sha256=sha(WITNESS),objective_exact=['-1','0'],lmin=108,mu_domain=['49','infinity'],tau_domain=['256','infinity'],kappa_exact='3000',directed_upper_exact=str(max(ups)),paired_precision_bits=[400,480],source_certificates=files,sha256_manifest=manifest)
 (OUT/'large_tau_complete.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='sha256_manifest'},indent=2))
if __name__=='__main__':main()
