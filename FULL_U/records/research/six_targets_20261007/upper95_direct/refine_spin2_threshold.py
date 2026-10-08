"""Replace only the overly conservative spin2 threshold cap, then full replay.
Interior cell formulas are identical for Centered and ConstantCentered.
Old interior enclosures remain valid; the newly exposed interval is integrated.
"""
from prepare_head import *
import copy
sys.path.insert(0,str(OLD))
from certify_head_constant import driver
def qdy(v):return Fraction(int(v[0]))*Fraction(2)**int(v[1])
def main():
    driver.ctx.prec=256;tic=time.time();source=BASE/'spin2/head_l2.json';raw=json.loads(source.read_text())
    assert raw['final']and raw['model']=='centered'and raw['spin']==2 and raw['witness_sha256']==EXPECTED
    for p,h in raw['dependencies'].items():assert sha(ROOT/p)==h,p
    assert sha(BASE/'moments128/manifest.json')==raw['source_manifest_sha256']
    model=driver.MODELS['constant_centered'](WITNESS,2,BASE/'moments128',N=128)
    gaps=model.gaps(30);oldgaps=[tuple(map(Fraction,x))for x in raw['gap_beta_intervals']]
    assert gaps[1:]==oldgaps[1:]and gaps[0][0]==oldgaps[0][0]==0 and gaps[0][1]<oldgaps[0][1]
    lo,hi=gaps[0][1],oldgaps[0][1];n=64;new=[]
    for i in range(n):
        a=lo+(hi-lo)*i/n;b=lo+(hi-lo)*(i+1)/n
        up,_=model.trapezoid(a,b);assert up.is_finite()and up>=0
        new.append(dict(a=str(a),b=str(b),upper_dyadic=driver.dyadic(up)))
    cells=sorted(copy.deepcopy(raw['cells_data'])+new,key=lambda c:Fraction(c['a']))
    driver.exact_partition([(Fraction(c['a']),Fraction(c['b']))for c in cells],gaps)
    gapmass=sum((model.gap_bound(a,b)for a,b in gaps),driver.arb(0))
    total=gapmass+sum((driver.exact(qdy(c['upper_dyadic']))for c in cells),driver.arb(0))
    result=copy.deepcopy(raw);result.update(model='constant_centered',upper=str(total.upper()),upper_dyadic=driver.dyadic(total),target_met=bool(total<driver.arb(raw['target'])),cells=len(cells),gap_upper=str(gapmass.upper()),gap_upper_dyadic=driver.dyadic(gapmass),gap_beta_intervals=[[str(a),str(b)]for a,b in gaps],threshold=model.cap_status,cells_data=cells,threshold_refinement_source=str(source.relative_to(ROOT)),threshold_refinement_source_sha256=sha(source),new_threshold_cells=n,old_interior_cells_retained=len(raw['cells_data']),elapsed_seconds=time.time()-tic,midpoint_estimate='Not recomputed: exact upper accounting only')
    result['dependencies'].update(driver.code_hashes());result['dependencies'].update({str(p.relative_to(ROOT)):sha(p)for p in[source,Path(__file__)]})
    out=BASE/'spin2_refined/head_l2.json';assert not out.exists();out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2))
    print('PASS_REFINED_CAP',str(total.upper()),'gap',str(gapmass.upper()),'seconds',time.time()-tic,flush=True)
if __name__=='__main__':main()
