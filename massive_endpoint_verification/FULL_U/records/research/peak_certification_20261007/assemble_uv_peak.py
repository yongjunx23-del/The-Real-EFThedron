"""Exact coverage and hash binding for a new directional UV coercivity proof.

This assembles only mu>=49 and even spins>=42. It asserts no support
integral, compact-energy theorem, or complete directional bound by itself.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,time
ROOT=Path(__file__).resolve().parents[2]
CONFIG={'upper83': {'sha': '68e39591f3f9a9d3b5cdfe3ade3924fecb227a329de4a9543729f0607e8f663a', 'objective': ['535119985632953/4503599627370496', '8943390007067311/9007199254740992'], 'panels': 16, 'groups': [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (9, 9), (10, 10), (11, 11), (12, 12), (13, 13), (14, 14), (15, 15)], 'kappa': '1/200', 'suffix': '1_200', 'uniform_stem': 'uniform'}, 'lower263': {'sha': '8931658c2a353a2bbf9b9cb074db490ba5d15f044c80837236d3fd9dd801c541', 'objective': ['-2109641626282615/18014398509481984', '-8945221711437375/9007199254740992'], 'panels': 1, 'groups': [(0, 0)], 'kappa': '1/50', 'suffix': '1_50', 'uniform_stem': 'uniform1'}}
GROUPS=[(42,90),(92,142),(144,192),(194,244),(246,294),(296,346),(348,396),(398,448)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--angle',choices=list(CONFIG),required=True);ap.add_argument('--scope',choices=['all','finite','infinite'],default='all');ap.add_argument('--verify-only',action='store_true');a=ap.parse_args();tic=time.time();cfg=CONFIG[a.angle];outdir=ROOT/'output/peak_certification_20261007'/a.angle/'tail';manifest={};seen=set()
 def ingest(path,expected=None):
  p=Path(path).resolve();key=str(p.relative_to(ROOT));h=sha(p)
  if expected is not None:assert h==expected,key
  if key in manifest:assert manifest[key]==h,key
  manifest[key]=h
  if p.suffix!='.json':return None
  d=json.loads(p.read_text())
  if key not in seen:
   seen.add(key)
   for q,hh in d.get('sha256_manifest',{}).items():ingest(ROOT/q,hh)
  return d
 def common(d,status):
  assert d['status']==status and d['witness_sha256']==cfg['sha']and d['objective_exact']==cfg['objective']
 def partition(cells,end,kappa):
  assert cells and F(cells[0]['a'])==0 and F(cells[-1]['b'])==end
  assert all(F(c['a'])<F(c['b'])and F(c['upper_exact'])<-kappa for c in cells)
  assert all(F(x['b'])==F(y['a'])for x,y in zip(cells,cells[1:]))
 def pair_records(p,r,key):
  assert len(p)==len(r)
  n=0;worst=None
  for x,y in zip(p,r):
   assert x[key]==y[key]and len(x['cells'])==len(y['cells'])
   for xx,yy in zip(x['cells'],y['cells']):
    assert F(xx['a'])==F(yy['a'])and F(xx['b'])==F(yy['b'])
    u=max(F(xx['upper_exact']),F(yy['upper_exact']));worst=u if worst is None else max(worst,u);n+=1
  return n,worst
 prep=ingest(outdir/'prepared_components.json');common(prep,f'PASS_ANGLE{a.angle}_PREPARED_COMPONENTS')
 source=ingest(outdir/'source_exact_rebuild.json');common(source,'PASS_EXACT_RATIONAL_SOURCE_REBUILD');assert source['all_entries_equal']is True and source['moment_count']==972 and source['source_remainder_coefficients']==9
 large=[]
 for bits in[320,400]:
  p=outdir/f'large_tau_{bits}.json';d=ingest(p);common(d,f'PASS_ANGLE{a.angle}_UNBOUNDED_TAU')
  assert d['precision_bits']==bits and d['lmin']==108 and d['mu_domain']==['49','infinity']and d['tau_domain']==['256','infinity']
  assert F(d['kappa_exact'])==1 and F(d['V_over_p_upper_exact'])<-1
  large.append(str(p.relative_to(ROOT)))
 assert F(106*107,45)<256
 finite=None;infinite=None
 if a.scope in['all','finite']:
  files=[];allspins=[];total=0;worst=None
  for lo,hi in GROUPS:
   pp=outdir/f'finite_{lo}_{hi}_400_k{cfg["suffix"]}.json';rp=outdir/f'finite_{lo}_{hi}_480_k{cfg["suffix"]}.json';p=ingest(pp);r=ingest(rp)
   for bits,d in[(400,p),(480,r)]:
    common(d,f'PASS_ANGLE{a.angle}_FINITE_UV');assert not d['failed']and d['precision_bits']==bits and F(d['kappa_exact'])==F(cfg['kappa'])
    assert d['spins']==list(range(lo,hi+1,2))and len(d['records'])==len(d['spins'])and all(F(x)<=0 for x in d['crossing_tail_coefficients'])
    for l,rec in zip(d['spins'],d['records']):
     end=min(F(256),F(l*(l+1),45));assert rec['spin']==l and rec['status']=='PASS'and F(rec['tau_end_exact'])==end
     partition(rec['cells'],end,F(cfg['kappa']))
   assert r['reference_certificate_sha256']==sha(pp)
   n,w=pair_records(p['records'],r['records'],'spin');total+=n;worst=w if worst is None else max(worst,w);allspins.extend(p['spins'])
   files.append(dict(spins=[lo,hi],producer=str(pp.relative_to(ROOT)),replay=str(rp.relative_to(ROOT)),producer_sha256=sha(pp),replay_sha256=sha(rp)))
  assert allspins==list(range(42,450,2))
  finite=dict(status='PASS_COMPLETE_DIRECTIONAL_FINITE_UV',witness_sha256=cfg['sha'],objective_exact=cfg['objective'],spins=allspins,first_even_spin=42,last_even_spin=448,mu_domain=['49','infinity'],minimum_kappa_exact=str(min(-worst,F(1))),weakest_paired_upper_exact=str(worst),cell_count=total,local_bound_rule='On each closed tau cell use minus max(producer upper_exact,replay upper_exact). For tau>=256 use kappa=1.',paired_files=files,large_tau_components=large,large_tau_kappa_exact='1')
 if a.scope in['all','infinite']:
  panels=[];files=[];total=0;worst=None
  for lo,hi in cfg['groups']:
   pp=outdir/f'{cfg["uniform_stem"]}_p{lo}_{hi}_400_k{cfg["suffix"]}.json';rp=outdir/f'{cfg["uniform_stem"]}_p{lo}_{hi}_480_k{cfg["suffix"]}.json';p=ingest(pp);r=ingest(rp)
   for bits,d in[(400,p),(480,r)]:
    common(d,f'PASS_ANGLE{a.angle}_UNIFORM_UV');assert not d['failed']and d['precision_bits']==bits and F(d['kappa_exact'])==F(cfg['kappa'])
    assert d['lmin']==450 and d['tau_domain']==['0','256']and d['eta_panels']==cfg['panels']and d['panel_indices']==list(range(lo,hi+1))
    assert len(d['records'])==len(d['panel_indices'])and all(F(x)<=0 for x in d['crossing_tail_coefficients'])
    for i,rec in zip(d['panel_indices'],d['records']):
     assert rec['panel']==i and rec['status']=='PASS'and F(rec['eta_lo_exact'])==F(i,cfg['panels']*450*451)and F(rec['eta_hi_exact'])==F(i+1,cfg['panels']*450*451)
     partition(rec['cells'],F(256),F(cfg['kappa']))
   assert r['reference_certificate_sha256']==sha(pp)
   n,w=pair_records(p['records'],r['records'],'panel');total+=n;worst=w if worst is None else max(worst,w);panels.extend(p['panel_indices'])
   files.append(dict(panels=[lo,hi],producer=str(pp.relative_to(ROOT)),replay=str(rp.relative_to(ROOT)),producer_sha256=sha(pp),replay_sha256=sha(rp)))
  assert panels==list(range(cfg['panels']))and F(4)+F(450*451,256)>49
  infinite=dict(status='PASS_COMPLETE_DIRECTIONAL_INFINITE_UV',witness_sha256=cfg['sha'],objective_exact=cfg['objective'],lmin=450,mu_domain=['49','infinity'],kappa_exact=cfg['kappa'],eta_panels=cfg['panels'],cell_count=total,weakest_paired_upper_exact=str(worst),paired_files=files,large_tau_components=large,large_tau_kappa_exact='1')
 for rel in [
  'research/all_direction_certification_20261007/tail/contracted_forward.py',
  'research/all_direction_certification_20261007/tail/finite_uv_bridge220.py',
  'research/all_direction_certification_20261007/tail/uniform_source_uv.py',
  'research/all_direction_certification_20261007/tail/finite_source_corrections.py',
  'research/all_direction_certification_20261007/tail/compact_highspin220.py',
  'research/lower_g0_upper_g2_20261006/tail_finite_uv_region.py',
  'research/lower_g0_upper_g2_20261006/tail_certify_joint_limit.py',
  'research/lower_g0_upper_g2_20261006/matrix_crossing_series.py',
  'research/lower_g0_upper_g2_20261006/matrix_source_norms.py',
  'research/lower_g0_upper_g2_20261006/matrix_roy_guards.py',
  'research/crossing_kernels.py',
  'research/tail_certification/certify_uv.py',
  'research/tail_certification/compact_source_data.py']:
  ingest(ROOT/rel)
 ingest(ROOT/'research/boundary_length_20261007/REVIEW_NEW_ENVELOPES.md')
 ingest(ROOT/'research/boundary_length_20261007/tail/UV_PROOF_NOTES.txt')
 ingest(Path(__file__).with_name('UV_TAIL_PROOF_AND_REVIEW.md'))
 ingest(outdir/'exact_rebuild_and_algebra_audit.json')
 ingest(outdir/'uniform_queue_results.json')
 ingest(outdir.parent/'tail_auxiliary_campaign.json')
 ingest(Path(__file__))
 for name,d in[('finite_uv_complete',finite),('infinite_uv_complete',infinite)]:
  if d:
   d.update(sha256_manifest=dict(manifest),seconds=time.time()-tic,scope='UV coercivity only; no compact-energy or support-integral claim.')
   if a.verify_only:
    saved=ingest(outdir/(name+'.json'));assert all(saved.get(k)==v for k,v in d.items()if k not in['seconds','sha256_manifest'])
   else:(outdir/(name+'.json')).write_text(json.dumps(d,indent=2))
   print(json.dumps({k:v for k,v in d.items()if k not in['sha256_manifest','paired_files','spins']},indent=2))
 if a.scope=='all':
  for name in['finite_uv_complete','infinite_uv_complete']:ingest(outdir/(name+'.json'))
  d=dict(status='PASS_COMPLETE_DIRECTIONAL_UV_COERCIVITY',witness_sha256=cfg['sha'],objective_exact=cfg['objective'],finite_first_even_spin=42,finite_last_even_spin=448,infinite_first_even_spin=450,finite_kappa_exact=finite['minimum_kappa_exact'],infinite_kappa_exact=infinite['kappa_exact'],energy_domain=['49','infinity'],finite_component=str((outdir/'finite_uv_complete.json').relative_to(ROOT)),infinite_component=str((outdir/'infinite_uv_complete.json').relative_to(ROOT)),finite_local_kappa_rule=finite['local_bound_rule'],sha256_manifest=manifest,scope='Every even spin>=42, all mu>=49. Combine separately with compact coercivity and same-witness real-norm bounds.')
  if a.verify_only:
   saved=ingest(outdir/'uv_complete.json');assert all(saved.get(k)==v for k,v in d.items()if k!='sha256_manifest');print('PASS_READ_ONLY_FROZEN_PACKET_VERIFICATION')
  else:(outdir/'uv_complete.json').write_text(json.dumps(d,indent=2))
if __name__=='__main__':main()
