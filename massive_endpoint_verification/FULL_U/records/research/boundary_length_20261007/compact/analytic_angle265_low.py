"""Low-energy allspin42 subdomain of the exact all-powers compact envelope."""
import analytic_angle265 as model
from pathlib import Path
import hashlib,json
model.PAIRS=[('4','9/2')]
model.DEST=model.DEST/'low'
model.DEST.mkdir(parents=True,exist_ok=True)
if __name__=='__main__':
    model.main()
    p=model.DEST/f'analytic_piecewise_{model.ctx.prec}.json'
    j=json.loads(p.read_text())
    assert j['records'][0]['minimum_even_spin']==42
    q=Path(__file__).resolve();j['sha256_manifest'][str(q.relative_to(model.ROOT))]=hashlib.sha256(q.read_bytes()).hexdigest()
    p.write_text(json.dumps(j,indent=2))
