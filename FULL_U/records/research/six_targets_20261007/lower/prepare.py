"""Define three exact rational mixtures and diagnose low-spin head gains."""
from pathlib import Path
from fractions import Fraction as F
import os,sys,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'output/six_targets_20261007/lower'
CASES={'lower15':('255','259'),'lower45':('259','261'),'lower60':('261','263')}
BASES={
 '255':('output/boundary_length_20261007/complete_angle255_banded_certificate.json','output/boundary_length_20261007/head/angle255'),
 '259':('output/boundary_length_20261007/complete_angle259_banded_certificate.json','output/boundary_length_20261007/head/angle259'),
 '261':('output/boundary_length_20261007/complete_angle261p3814_certificate.json','output/boundary_length_20261007/head/angle261p3814'),
 '263':('output/peak_certification_20261007/complete_lower263_certificate.json','output/peak_certification_20261007/lower263/head')}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dy(v):return F(int(v[0]))*F(2)**int(v[1])
def cost(U,V,W):
 h=np.hypot(U,V)+V;n=V<0;h[n]=U[n]*(U[n]/(np.hypot(U[n],V[n])-V[n]));return h*W
def per_spin(values,meta):
 offset=0;result={}
 for block in meta['blocks']:
  ns,ne=len(block['spins']),block['energy_nodes'];end=offset+ns*ne
  rows=values[offset:end].reshape(ns,ne)
  result.update({str(ell):float(rows[j].sum())for j,ell in enumerate(block['spins'])});offset=end
 assert offset==len(values);return result
def main():
 targets=json.loads((ROOT/'output/six_targets_20261007/target_map.json').read_text());targets={r['tag']:r for r in targets['targets']}
 prepared=[]
 for tag,labels in CASES.items():
  folder=OUT/tag;folder.mkdir(parents=True,exist_ok=True)
  assert not (folder/'definition.json').exists()
  cc=[json.loads((ROOT/BASES[x][0]).read_text())for x in labels]
  wp=[ROOT/c['witness']for c in cc];ww=[np.load(p)for p in wp]
  assert all(sha(p)==c['witness_sha256']for p,c in zip(wp,cc))
  assert np.array_equal(ww[0]['names'],ww[1]['names'])and len(ww[0]['names'])==220
  base_d=[list(map(F,c['objective_exact']))for c in cc]
  assert all([F(float(v))for v in w['direction']]==b for w,b in zip(ww,base_d))
  d=list(map(lambda x:F(float(x)),targets[tag]['physical_outward_unit_normal']))
  a,b=base_d[0];c,e=base_d[1];det=a*e-b*c
  weights=[(d[0]*e-d[1]*c)/det,(a*d[1]-b*d[0])/det]
  assert min(weights)>0 and [sum(weights[i]*base_d[i][j]for i in[0,1])for j in[0,1]]==d
  multipliers=[sum(weights[i]*F(float(ww[i]['multipliers'][j]))for i in[0,1])for j in range(220)]
  heads=[json.loads((ROOT/BASES[x][1]/'head_summary.json').read_text())for x in labels]
  rows=[{r['spin']:r for r in s['records']}for s in heads]
  oldcosts=[]
  for i,s in enumerate(heads):
   costs={}
   for ell,r in rows[i].items():
    p=json.loads((ROOT/r['path']).read_text());q=json.loads((ROOT/r['replay']).read_text());costs[ell]=max(dy(p['upper_dyadic']),dy(q['upper_dyadic']))
   assert sum(costs.values(),F(0))==F(cc[i]['components_exact']['head_0_40']);oldcosts.append(costs)
  definition=dict(status='EXACT_RATIONAL_WITNESS_DEFINITION_NOT_A_CERTIFICATE',witness_format='exact_rational_positive_combination',tag=tag,target_station=targets[tag],base_labels=list(labels),base_certificate_paths=[BASES[x][0]for x in labels],base_certificate_sha256=[sha(ROOT/BASES[x][0])for x in labels],base_witness_paths=[str(p.relative_to(ROOT))for p in wp],base_witness_sha256=[sha(p)for p in wp],base_head_directories=[BASES[x][1]for x in labels],weights_exact=list(map(str,weights)),objective_exact=list(map(str,d)),names=ww[0]['names'].tolist(),multipliers_exact=list(map(str,multipliers)),base_retained_costs=[{str(k):str(v)for k,v in costs.items()}for costs in oldcosts],old_convex_total_exact=str(sum(weights[i]*F(cc[i]['support_upper_exact'])for i in[0,1])),preparation_source_sha256=sha(Path(__file__)))
  (folder/'definition.json').write_text(json.dumps(definition,indent=2));prepared.append(definition)
 for order in[192,384]:
  cache=ROOT/f'output/full_direction_20261007/right/dense220_n{order}k{96 if order==192 else 128}_cache.npz';z=np.load(cache);meta=json.loads(cache.with_suffix('.json').read_text())
  for d in prepared:
   assert d['names']==z['names'].tolist()
   a,b=map(float,map(F,d['objective_exact']));lam=np.array([float(F(v))for v in d['multipliers_exact']]);U=-a*z['ub']+lam@z['U'];V=-a*z['vb']+b*z['p']+lam@z['V'];values=cost(U,V,z['W']);spins=per_spin(values,meta)
   weights=list(map(F,d['weights_exact']));base=[{int(k):F(v)for k,v in row.items()}for row in d['base_retained_costs']]
   old=F(d['old_convex_total_exact']);removed=sum(weights[i]*sum((base[i][ell]for ell in[0,2]),F(0))for i in[0,1]);estimate=float(old-removed)+spins['0']+spins['2']
   report=dict(status='NUMERICAL_MIXTURE_DIAGNOSTIC_NOT_A_BOUND',tag=d['tag'],n=order,head_total=float(values.sum()),per_spin=spins,old_convex=float(old),estimated_complete_with_only_spins0_2_reintegrated=estimate,old_polygon_support=d['target_station']['old_certified_support'],estimated_improvement_vs_old_polygon=d['target_station']['old_certified_support']-estimate,cache=str(cache.relative_to(ROOT)),cache_sha256=sha(cache))
   (OUT/d['tag']/f'numerical_n{order}.json').write_text(json.dumps(report,indent=2));print(d['tag'],order,'head',report['head_total'],'prospective',estimate,'gain',report['estimated_improvement_vs_old_polygon'],flush=True)
if __name__=='__main__':main()
