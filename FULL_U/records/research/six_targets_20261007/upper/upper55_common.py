"""Exact enlargement of a source gap using existing nonnegative cell bounds."""
from fractions import Fraction as F
import copy

def dy(pair):return F(int(pair[0]))*F(2)**int(pair[1])
def as_dyadic(v):
    assert v.denominator & (v.denominator-1) == 0
    return [str(v.numerator),str(-(v.denominator.bit_length()-1))]

def align_gaps(records,ell):
    records=copy.deepcopy(records)
    if ell==0:
        assert records[0]['gap_beta_intervals']==records[1]['gap_beta_intervals']
        return records,dict(method='identical_gap_partitions',moved_cells=[])
    assert ell==2
    a,b=records
    assert list(map(F,a['gap_beta_intervals'][0]))==[F(0),F(1,1024)]
    assert list(map(F,b['gap_beta_intervals'][0]))==[F(0),F(1,256)]
    assert a['gap_beta_intervals'][1:]==b['gap_beta_intervals'][1:]
    moved=[];kept=[];cursor=F(1,1024);added=F(0)
    for cell in a['cells_data']:
        lo,hi=F(cell['a']),F(cell['b'])
        if hi<=F(1,256):
            assert lo==cursor<hi;cursor=hi;v=dy(cell['upper_dyadic']);assert v>=0;added+=v;moved.append(cell)
        else:
            assert lo>=F(1,256);kept.append(cell)
    assert cursor==F(1,256)and len(moved)==64
    oldgap=dy(a['gap_upper_dyadic']);newgap=oldgap+added
    a['gap_beta_intervals']=copy.deepcopy(b['gap_beta_intervals'])
    a['gap_upper_dyadic']=as_dyadic(newgap);a['cells_data']=kept;a['cells']=len(kept)
    assert dy(a['gap_upper_dyadic'])+sum((dy(c['upper_dyadic'])for c in kept),F(0))<=dy(a['upper_dyadic'])
    return records,dict(method='extend_threshold_gap_with_original_cell_bounds',old_gap_exact=str(oldgap),added_cell_bounds_exact=str(added),new_gap_exact=str(newgap),moved_cells=moved,old_endpoint='1/1024',new_endpoint='1/256')
