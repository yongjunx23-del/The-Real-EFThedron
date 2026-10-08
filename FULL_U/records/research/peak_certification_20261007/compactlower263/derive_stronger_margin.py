"""Derive kappa3 only after checking every paired compact enclosure."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,copy
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'output/peak_certification_20261007/lower263/compact';K=F(3)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_bound(path,manifest):
    path=path.resolve();d=json.loads(path.read_text());manifest[str(path.relative_to(ROOT))]=sha(path)
    for p,h in d.get('sha256_manifest',{}).items():assert sha(ROOT/p)==h,p
    return d

def main():
    source=OUT/'complete_compact_coercivity.json';base=json.loads(source.read_text());assert base['status']=='PASS_COMPLETE_COMPACT_COERCIVITY'and F(base['kappa_exact'])==1
    manifest=dict(base['sha256_manifest'])
    for p,h in manifest.items():assert sha(ROOT/p)==h,p
    manifest[str(source.relative_to(ROOT))]=sha(source);manifest[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    witness=base['witness_sha256'];ledger=[];all_uppers=[];energy_end=F(4)
    for c in base['coverage_rectangles']:
        lo,hi=map(F,c['mu_domain_exact']);assert lo==energy_end and lo<hi;energy_end=hi;L=c['analytic_minimum_even_spin'];rec={'mu_domain_exact':c['mu_domain_exact'],'analytic_minimum_even_spin':L,'analytic':[],'finite':[]}
        for rel in c['analytic_pair']:
            d=load_bound(ROOT/rel,manifest);assert d['witness_sha256']==witness and d['objective_exact']==base['objective_exact'];assert d['status']=='PASS_PIECEWISE_ALLSPIN_COMPACT_ENVELOPES'
            rows=[r for r in d['records']if (r['mu_domain_exact']==c['analytic_record'] if 'analytic_record'in c else len(d['records'])==1)];assert len(rows)==1;r=rows[0];a,b=map(F,r['mu_domain_exact']);assert a<=lo<hi<=b and r['minimum_even_spin']==L
            up=F(r['V_over_p_upper_exact']);assert up<-K
            if r['method']=='low':assert F(r['spin_derivative_upper_exact'])<0
            all_uppers.append(up);rec['analytic'].append({'source':rel,'precision_bits':d['precision_bits'],'V_over_p_upper_exact':str(up)})
        fin=c.get('finite_component',{'spins':c.get('finite_spins',[])})
        assert fin['spins']==list(range(42,L,2))
        if L>42:
            for which in ['producer','replay']:
                rel=fin[which];d=load_bound(ROOT/rel,manifest);assert d['witness_sha256']==witness and d['spins']==fin['spins'] and list(map(F,d['mu_domain']))==[lo,hi];assert d['status']=='PASS_FINITE_COMPACT_BRIDGE'
                worst=None;n=0
                for ell,cells in d['cells_by_spin'].items():
                    assert int(ell) in fin['spins'];last=lo
                    for cell in sorted(cells,key=lambda q:F(q['a'])):
                        a,b,up=map(F,[cell['a'],cell['b'],cell['V_over_p_upper_exact']]);assert a==last and a<b and up<-K;last=b;n+=1;worst=up if worst is None else max(up,worst)
                    assert last==hi
                assert set(map(int,d['cells_by_spin']))==set(fin['spins']);all_uppers.append(worst);rec['finite'].append({'source':rel,'precision_bits':d['precision_bits'],'checked_spin_cells':n,'V_over_p_upper_exact':str(worst)})
        ledger.append(rec)
    assert energy_end==49;worst=max(all_uppers);assert worst<-K
    margin={'status':'PASS_STRONGER_COMPACT_MARGIN_FROM_ALL_SAVED_ENCLOSURES','source_complete_component':str(source.relative_to(ROOT)),'source_sha256':sha(source),'witness_sha256':witness,'objective_exact':base['objective_exact'],'kappa_exact':str(K),'worst_V_over_p_upper_exact':str(worst),'minimum_verified_margin_exact':str(-worst),'coverage':ledger,'derivation':'Every original and replay analytic rectangle and every original and replay finite spin cell is strictly below -3. The exact adjacent energy cover and every even-spin cover are independently rechecked here. Shared source mathematics is unchanged.'}
    lp=OUT/'stronger_margin_ledger.json';assert not lp.exists();lp.write_text(json.dumps(margin,indent=2));manifest[str(lp.relative_to(ROOT))]=sha(lp)
    result=copy.deepcopy(base);result.update(kappa_exact=str(K),source_complete_component=str(source.relative_to(ROOT)),source_complete_component_sha256=sha(source),source_kappa_exact='1',stronger_margin_ledger=str(lp.relative_to(ROOT)),minimum_verified_margin_exact=str(-worst),scope='Only compact coercivity V/p<-3 for every even spin>=42 and4<mu<=49. No real norm, UV, head integral or complete directional support.',sha256_manifest=manifest)
    dest=OUT/'complete_compact_coercivity_stronger.json';assert not dest.exists();dest.write_text(json.dumps(result,indent=2));print('PASS kappa',K,'minimum margin',float(-worst),'files',len(manifest),'SHA256',sha(dest))
if __name__=='__main__':main()
