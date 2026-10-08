"""Remove inapplicable normalized-source estimates from shifted spin2 records.

Raw records remain unchanged. Each affected contribution is replaced by its
authenticated original base-cell upper, valid under the nonpositive shift.
"""
from pathlib import Path
from fractions import Fraction as F
import copy,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/weak_regions_20261007'
LIMIT=F(79,119)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dy(x):return F(int(x[0]))*F(2)**int(x[1])

def corrected_record(raw,base,raw_path):
    assert raw['status']=='PASS_EXACT_SHIFT_FINITE_HEAD_COMPONENT' and raw['spin']==base['spin']==2
    assert raw['source_record_sha256']==sha(ROOT/raw['source_record'])
    assert base['status']=='PASS_FINITE_SPIN_ENCLOSURE' and base['final']
    assert base['witness_sha256']==raw['base_witness_sha256']
    assert all(F(x)<0 for x in raw['g2_shift_exact'])
    assert len(raw['cells'])==len(base['cells_data'])
    result=copy.deepcopy(raw);changes=[]
    for index,(cell,original) in enumerate(zip(result['cells'],base['cells_data'])):
        a,b=F(cell['a']),F(cell['b']);assert(a,b)==(F(original['a']),F(original['b']))
        original_upper=dy(original['upper_dyadic']);assert original_upper>=0
        for j,method in enumerate(cell['methods']):
            if method=='shifted_negative_V' and b*b>LIMIT:
                previous=F(cell['upper_exact'][j]);assert 0<=previous<=original_upper
                changes.append(dict(cell_index=index,objective_index=j,a=str(a),b=str(b),old_method=method,old_upper_exact=str(previous),new_method='base_monotonicity',new_upper_exact=str(original_upper),increase_exact=str(original_upper-previous)))
                cell['methods'][j]='base_monotonicity';cell['upper_exact'][j]=str(original_upper)
    count=len(raw['objective_exact']);gap=F(raw['base_gap_upper_exact'])
    values=[gap+sum((F(c['upper_exact'][j])for c in result['cells']),F(0))for j in range(count)]
    assert all(values[j]>=F(raw['upper_exact'][j])for j in range(count))
    result['upper_exact']=list(map(str,values));result['upper_float']=list(map(float,values))
    result['dependencies'][str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    result['correction']=dict(status='PASS_DETERMINISTIC_MONOTONICITY_REPAIR',raw_record=str(Path(raw_path).resolve().relative_to(ROOT)),raw_record_sha256=sha(raw_path),repair_script=str(Path(__file__).relative_to(ROOT)),repair_script_sha256=sha(__file__),source_activity_beta_squared_limit_exact=str(LIMIT),source_activity_energy_cutoff_exact='119/10',affected_objective_cell_count=len(changes),affected_distinct_cell_count=len(set(c['cell_index']for c in changes)),changes=changes,increase_exact=[str(values[j]-F(raw['upper_exact'][j]))for j in range(count)],scope='Only selected spin2 normalized-negative-V choices above the Roy source endpoint are replaced. Every other cell, gap, coefficient and objective is unchanged. Raw outputs are invalidated separately and retained for provenance.')
    assert not any(c['methods'][j]=='shifted_negative_V' and F(c['b'])**2>LIMIT for c in result['cells']for j in range(count))
    return result

def main():
    dest=OUT/'head_corrected';dest.mkdir(exist_ok=True)
    records=[]
    for bits in [256,384]:
        raw_path=OUT/'head'/f'head_l2_{bits}.json';raw=json.loads(raw_path.read_text());base=json.loads((ROOT/raw['source_record']).read_text())
        result=corrected_record(raw,base,raw_path);path=dest/f'head_l2_{bits}.json';assert not path.exists();path.write_text(json.dumps(result,indent=2))
        records.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),bits=bits,**result['correction']))
        print(bits,'affected cells',result['correction']['affected_distinct_cell_count'],'objective cells',result['correction']['affected_objective_cell_count'],'increase',list(map(float,map(F,result['correction']['increase_exact']))))
    invalid=[]
    for name in ['complete_shift269p0404.json','complete_shift269p5.json']:
        path=OUT/name;assert path.exists();invalid.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path),reason='Depends on raw shifted spin2 choices that applied the Centered normalized real-source branch above its activity cutoff.'))
    notice=dict(status='RAW_SHIFT_CERTIFICATES_INVALIDATED_USE_CORRECTED_REPLACEMENTS',invalidated_complete_certificates=invalid,corrected_head_records=records,unaffected=['complete_fresh269p5.json','all prior frozen baseline certificates'],explanation='The original Centered head integrator used safe applicability conditions. A new shifted-head wrapper broadened them incorrectly. The repair restores valid original-cell monotonicity bounds on every affected selected cell; it makes no assumption that the numerical error was harmless. Raw PASS labels are retained only as historical provenance and must not be selected as valid certificates.')
    (dest/'CORRECTION_NOTICE.json').write_text(json.dumps(notice,indent=2))
if __name__=='__main__':main()
