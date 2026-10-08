"""Exact replacement of the finite head by its fully replayed cap refinement."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math
from assemble_upper95_direct import Inputs,ROOT,BASE,HEAD,SHA,sha,head_cost,ceiling_decimal
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,default=ROOT/'output/six_targets_20261007/complete_upper95_direct_certificate.json')
    ap.add_argument('--out',type=Path,default=ROOT/'output/six_targets_20261007/complete_upper95_direct_refined_certificate.json')
    args=ap.parse_args();assert not args.out.exists();inputs=Inputs()
    old=inputs.read(args.source);assert old['status']=='PASS_COMPLETE_CONTINUUM_DIRECTIONAL_BOUND'and old['witness_sha256']==SHA
    direction=old['objective_exact'];parts={k:F(v)for k,v in old['components_exact'].items()}
    assert all(v>=0 for v in parts.values())and sum(parts.values(),F(0))==F(old['support_upper_exact'])
    original=inputs.read(HEAD/'head_summary.json');selected=inputs.read(HEAD/'head_summary_selected.json')
    oldhead=head_cost(inputs,HEAD/'head_summary.json',SHA,direction)
    newhead=head_cost(inputs,HEAD/'head_summary_selected.json',SHA,direction)
    assert oldhead==parts['head_0_40']and newhead<oldhead
    assert [r['spin']for r in original['records']]==[r['spin']for r in selected['records']]==list(range(0,41,2))
    for x,y in zip(original['records'],selected['records']):
        if x['spin']!=2:assert x==y
    review=inputs.read(BASE/'cap_refinement_peer_review.json')
    assert review['status']=='PASS_THRESHOLD_REFINEMENT_IDENTITY_COVERAGE_AND_ACCOUNTING_REVIEW'and review['witness_sha256']==SHA
    assert review['new_cells']==64 and review['new_threshold_beta_upper_exact']=='1/1024'and review['old_threshold_beta_upper_exact']=='1/256'
    for path in [HEAD/'spin2/head_l2.json',HEAD/'spin2_refined/head_l2.json']:
        assert review['sha256_manifest'][str(path.relative_to(ROOT))]==sha(path)
    commands=inputs.read(HEAD/'spin2_refined/commands.json')
    assert commands['status']=='PASS_REFINED_SPIN2_AND_FULL_REPLAY'and commands['witness_sha256']==SHA
    assert len(commands['commands'])==2 and all(x['exit_code']==0 for x in commands['commands'])
    inputs.bind(ROOT/'research/six_targets_20261007/upper95_direct/run_threshold_refinement.py',commands['script_sha256'])
    for p in [Path(__file__),ROOT/'research/six_targets_20261007/assemble_upper95_direct.py',ROOT/'research/boundary_length_20261007/assemble_complete_v3.py',ROOT/'research/six_targets_20261007/upper95_direct/refinement_adaptation_manifest.json',ROOT/'research/six_targets_20261007/REVIEW_UPPER95_THRESHOLD_REFINEMENT.md']:
        inputs.bind(p)
    parts['head_0_40']=newhead;total=sum(parts.values(),F(0))
    assert total==F(old['support_upper_exact'])-oldhead+newhead
    assert total-newhead==F(old['tail_upper_exact'])
    dx,dy=map(F,direction);norm2=dx*dx+dy*dy;scale=2**256
    normlo=F(math.isqrt((norm2.numerator*scale*scale)//norm2.denominator),scale)
    assert normlo>0 and normlo**2<=norm2<(normlo+F(1,scale))**2
    assert str(norm2)==old['direction_norm_squared_exact']and str(normlo)==old['direction_norm_lower_exact']
    outward=ceiling_decimal(total,9);normalized=ceiling_decimal(F(outward)/normlo,9)
    result=dict(old);result.update(inequality=f'({direction[0]}) g0+({direction[1]}) g2 <= {outward}',support_upper_exact=str(total),support_upper_outward_9dp=outward,support_normalized_upper_outward_9dp=normalized,support_normalized_float=float(normalized),components_exact={k:str(v)for k,v in parts.items()},components_float={k:float(v)for k,v in parts.items()},head_threshold_refinement=dict(source_certificate=str(args.source.resolve().relative_to(ROOT)),source_certificate_sha256=sha(args.source),replaced_head_exact=str(oldhead),replacement_head_exact=str(newhead),improvement_exact=str(oldhead-newhead),improvement_float=float(oldhead-newhead),scope='All unchanged spins retain their complete enclosures; spin2 uses max of the refined256-bit producer and full384-bit replay. No source-window extrapolation.'),sha256_manifest=inputs.manifest)
    args.out.write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k]for k in ['status','support_normalized_float','components_float','head_threshold_refinement']},indent=2))
if __name__=='__main__':main()
