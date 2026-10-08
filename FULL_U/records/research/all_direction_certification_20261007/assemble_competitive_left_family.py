"""Attach the complete unchanged left tail to verified exact objective shifts.

This is a consequence of H(U,V-b*p)<=H(U,V), not an optimization claim.
Only complete finite-head packets enter the final list.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util, hashlib, json, math, sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research'))
from angle_benchmarks import support_report
SHA='e30dca80d94b35fd8c44b78bea02cd3ffe4b344ccdb1eb40b3da655a89f8cc4b'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def upward(x,digits):
    scale=10**digits;n=-((-x.numerator*scale)//x.denominator)
    return f'{n//scale}.{n%scale:0{digits}d}'
def main():
    # Re-execute the packet validator rather than trusting its previous status.
    module_path=ROOT/'research/all_direction_certification_20261007/head/assemble_left_family.py'
    spec=importlib.util.spec_from_file_location('family_validator',module_path)
    validator=importlib.util.module_from_spec(spec);spec.loader.exec_module(validator);validator.main()
    manifest={};seen=set()
    def bind(p,h=None):
        p=Path(p);p=(p if p.is_absolute()else ROOT/p).resolve();key=str(p.relative_to(ROOT));value=sha(p)
        if h is not None:assert value==h,key
        if key in manifest:assert manifest[key]==value,key
        manifest[key]=value;return p
    def read(p,h=None):
        p=bind(p,h);d=json.loads(p.read_text());key=str(p)
        if key not in seen:
            seen.add(key)
            for q,v in d.get('sha256_manifest',{}).items():
                q=bind(q,v)
                if q.suffix=='.json':read(q,v)
        return d
    base=read('output/competitive_20261007/left/complete_left_certificate.json')
    assert base['status']=='PASS_COMPLETE_CONTINUUM_LEFT_BOUND'and base['witness_sha256']==SHA
    assert base['objective_exact']==['-1','0']
    tail=F(base['tail_upper_exact']);assert tail>=0
    heads=read('output/all_direction_certification_20261007/head/left/family_packets.json')
    assert heads['base_witness_sha256']==SHA
    bind(module_path,heads['assembly_sha256'])
    for c in heads['source_checks']:
        assert c['exact_rebuild_match'];bind(c['path'],c['sha256']);bind(c['rebuilt_path'],c['sha256'])
    results=[];omitted=[]
    for p in heads['packets']:
        b=F(p['slope']);assert b>0 and p['exact_objective']==['-1',str(-b)]
        if p['status']!='COMPLETE_REPLAYED_FINITE_HEAD':
            omitted.append(dict(slope=str(b),missing_spins=p['missing_spins']));continue
        assert p['base_witness_sha256']==SHA and not p['missing_spins']
        assert p['included_spins']==list(range(0,41,2))
        bind(p['base_witness_path'],SHA);bind(p['source_manifest_path'],p['source_manifest_sha256'])
        for name,h in p['dependency_sha256'].items():bind(name,h)
        total=F(0)
        for c in p['selected_components']:
            r=read(c['component_path'],c['component_sha256']);k=c['selected_slope_index']
            assert r['spin']==c['spin']and F(r['slopes'][k])==b
            value=validator.dyadic(r['upper_dyadic'][k]);assert value==F(c['upper_exact']);total+=value
            read(c['source_record_path'],c['source_record_sha256']);read(c['original_replay_path'],c['original_replay_sha256'])
        assert total==F(p['head_upper_exact'])
        cost=total+tail;rounded=upward(cost,6);C=F(rounded);assert C>=cost
        # Rational lower bound for sqrt(1+b*b), hence upper unit-normal support.
        scale=10**30;n=math.isqrt((b.numerator**2+b.denominator**2)*scale**2)
        norm_lower=F(n,b.denominator*scale);assert 0<norm_lower and norm_lower**2<=1+b*b
        unit=upward(C/norm_lower,12)
        theta=(math.degrees(math.atan2(-float(b),-1))%360)
        dual=support_report(theta);primal=support_report(theta,'primal')
        results.append(dict(direction_exact=['-1',str(-b)],slope_exact=str(b),theta_degrees=theta,
            head_upper_exact=str(total),unchanged_tail_upper_exact=str(tail),support_upper_exact=str(cost),
            support_upper_outward_6dp=rounded,unit_normal_support_upper_outward_12dp=unit,
            graphical_comparison=dict(dual=dual,primal=primal,relative_to_dual_estimate_percent=100*(float(unit)/dual['support']-1))))
    assert results,'No complete family member is available.'
    bind(__file__);bind('research/all_direction_certification_20261007/REVIEW_LEFT_OBJECTIVE_FAMILY.md')
    out=dict(status='PASS_COMPLETE_CONTINUUM_LEFT_FAMILY',base_witness_sha256=SHA,
        witness_definition='Unchanged exact 220-row multipliers; objective(-1,-b). U unchanged; V shifts by -b*p.',
        justification='p>0 and H(U,V)=sqrt(U^2+V^2)+V is nondecreasing in V, so the complete base tail bounds every b>0 shifted tail.',
        coverage='Every even spin ell>=0, every energy mu>4 for each listed result.',
        caveat='Interval-enclosed fixed witness consequences, not independently optimized angles. Shifted finite heads use directed interval arithmetic; base partitions and gaps have precision replays.',
        results=results,omitted_incomplete=omitted,sha256_manifest=manifest)
    dest=ROOT/'output/competitive_20261007/left/complete_left_family_certificate.json';dest.write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items()if k not in ['sha256_manifest','results']},indent=2))
    print(json.dumps([{k:r[k]for k in ['slope_exact','support_upper_outward_6dp','unit_normal_support_upper_outward_12dp']}for r in results],indent=2))
if __name__=='__main__':main()
