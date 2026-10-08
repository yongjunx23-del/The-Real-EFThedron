"""Higher-precision UV-envelope replay, and compact-envelope replay.

Only the declared positive-direction cases passing the all-spin proof are
accepted. A 1e-10 allowance avoids assuming nesting of Arb enclosures.
"""
import json,pickle,sys,hashlib
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent))
import certify_directional_tail as m
m.CASES=[(F(1,10),F(2,5)),(F(1,4),F(1))]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    inp=ROOT/'output/directional_certification/tail_family_4b.json'
    saved=json.loads(inp.read_text())
    for p,h in saved['sha256_manifest'].items():assert sha(ROOT/p)==h
    data=pickle.loads((ROOT/'output/certification_tail/source_cache_v2_depth10.pkl').read_bytes())['data']
    compact=pickle.loads((ROOT/'output/certification_tail/compact_source_cache.pkl').read_bytes())
    m.ctx.prec=320
    uv,_=m.uv_family(data)
    results=[];slack=F(1,10**10)
    for b,d in m.CASES:
        old=next(r for r in saved['results'] if F(r['b'])==b)
        assert old['status']=='PASS' and F(old['delta_X4'])==d
        changed=dict(data);changed['lx']+=d
        cc=m.compact_run(42,changed,compact,'.02')
        values=[F(uv[b]['upper_exact']),m.dyad(m.arb(cc['low_energy_upper'])+m.ar(b)),m.dyad(m.arb(cc['middle_energy_upper'])+m.ar(b))]
        bound=max(values);old_bound=-F(old['tail_kappa_exact'])
        assert bound<=old_bound+slack<0
        results.append(dict(b=str(b),delta_X4=str(d),fresh_global_upper_exact=str(bound),
            accepted_global_upper_exact=str(old_bound+slack),kappa_exact=str(-old_bound-slack)))
    out=dict(status='PASS',claim='Independent invocation of all positive passing tail envelopes; UV replay at320bits and compact replay at frozen kernel precision.',
        results=results,uniform_slack_exact=str(slack),input_sha256=sha(inp),source_sha256=sha(__file__))
    target=ROOT/'output/directional_certification/tail_replay.json';target.write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
