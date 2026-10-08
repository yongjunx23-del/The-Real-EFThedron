"""Continuous exact-rational geometry on the observed digitized dual branches.

No new physical certificate is claimed. Segment parameters, halfplanes,
line-projection domains and envelope breakpoints are exact rationals. Arb
is used only for square-root metric lengths and the associated ratios.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,bisect,csv,hashlib,json,math,sys,time
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tmp/certification_deps'))
from flint import arb,ctx

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ar(q):q=F(q);return arb(q.numerator)/arb(q.denominator)
def endpoint(v):m,e=v.man_exp();return F(int(m))*F(2)**int(e)
def enclosure(v):
 assert v.is_finite()
 return dict(lower_exact=str(endpoint(v.lower())),upper_exact=str(endpoint(v.upper())),midpoint_float=float(v.mid()))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def clip(interval,m,b):
 """Intersect a closed rational t interval with m*t<=b."""
 if interval is None:return None
 lo,hi=interval
 if m>0:hi=min(hi,b/m)
 elif m<0:lo=max(lo,b/m)
 elif b<0:return None
 return None if lo>hi else(lo,hi)
def size(interval):return F(0)if interval is None else interval[1]-interval[0]
def signed_key(u,L2):
 s=(u>0)-(u<0)
 return s,s*u*u*L2

def polygon_audit(cuts,saved):
 vertices=set()
 for i,(a,b,c)in enumerate(cuts):
  assert a or b
  for d,e,f in cuts[i+1:]:
   det=a*e-b*d
   if not det:continue
   v=((c*e-b*f)/det,(a*f-c*d)/det)
   if all(x*v[0]+y*v[1]<=z for x,y,z in cuts):vertices.add(v)
 assert vertices==set(tuple(map(F,v))for v in saved),'Saved polygon differs from exact halfplane intersection'
 assert len(vertices)>=3
 vv=list(vertices);assert any((q[0]-vv[0][0])*(r[1]-vv[0][1])!=(q[1]-vv[0][1])*(r[0]-vv[0][0])for q in vv[1:]for r in vv[1:])
 # A nontrivial planar polyhedral recession cone has a boundary ray parallel
 # to a constraint line, unless it is the whole plane (already excluded).
 for a,b,c in cuts:
  for v in[(b,-a),(-b,a)]:assert any(x*v[0]+y*v[1]>0 for x,y,z in cuts),'Unbounded recession direction'
 return dict(exact_pairwise_vertices=len(vertices),saved_vertices_equal=True,bounded_recession_cone=True,full_dimensional=True)

def segment(p,q,side,scale,cuts,index):
 delta=tuple(y-x for x,y in zip(p,q));d=tuple(x/s for x,s in zip(delta,scale));L2=dot(d,d);assert L2>0
 n=(-d[1],d[0])if side=='upper'else(d[1],-d[0]);ray=tuple(x*s for x,s in zip(n,scale))
 upper=[];lower=[];parallel=[];inside=(F(0),F(1))
 for k,(a,b,C)in enumerate(cuts):
  A=(a,b);v=C-dot(A,p);E=dot(A,delta);D=dot(A,ray)
  inside=clip(inside,E,v)
  if D==0:parallel.append((E,v,k))
  elif D>0:upper.append((v/D,-E/D,k))
  else:lower.append((v/D,-E/D,k))
 assert upper and lower,'Bounded full-dimensional polygon has both normal-line endpoints'
 domain=(F(0),F(1))
 for E,v,k in parallel:domain=clip(domain,E,v)
 for qa,ra,ka in lower:
  for qb,rb,kb in upper:domain=clip(domain,ra-rb,qb-qa)
 data=dict(segment=index,start_exact=list(map(str,p)),end_exact=list(map(str,q)),metric_length_squared_exact=str(L2),inside_reference_t=None if inside is None else list(map(str,inside)),normal_line_domain_t=None if domain is None else list(map(str,domain)),reference_excluded_parameter_fraction_exact=str(1-size(inside)),normal_line_missing_parameter_fraction_exact=str(1-size(domain)),metric_normal_unnormalized_exact=list(map(str,n)),physical_ray_unnormalized_exact=list(map(str,ray)),physical_normal_angle_degrees=math.degrees(math.atan2(float(n[1]/scale[1]),float(n[0]/scale[0])))%360)
 internal=dict(L2=L2,p=p,delta=delta,n=n,ray=ray,upper=upper,lower=lower,domain=domain,inside=inside,data=data)
 if domain is None:
  data.update(upper_envelope_pieces=[],maximum_signed_gap=None);return internal
 lo,hi=domain;candidates={lo,hi}
 for j,(q0,r0,k0)in enumerate(upper):
  for q1,r1,k1 in upper[j+1:]:
   if r0!=r1:
    t=(q1-q0)/(r0-r1)
    if lo<t<hi:candidates.add(t)
 def val(t):return min(q0+r0*t for q0,r0,k0 in upper)
 scored=[(val(t),t)for t in sorted(candidates)];peak=max(scored)
 knots={lo,hi}
 for u,t in scored:
  active=[x for x in upper if x[0]+x[1]*t==u]
  if len({x[1]for x in active})>1:knots.add(t)
 knots=sorted(knots);pieces=[]
 for a,b in zip(knots,knots[1:]):
  mid=(a+b)/2;q0,r0,k0=min(upper,key=lambda x:x[0]+x[1]*mid)
  assert all(q0+r0*t<=q1+r1*t for t in[a,b]for q1,r1,k1 in upper)
  pieces.append(dict(t_interval_exact=[str(a),str(b)],active_cut=k0,intercept_exact=str(q0),slope_exact=str(r0)))
 u,t=peak;assert max(q0+r0*t for q0,r0,k0 in lower)<=u
 at=tuple(x+t*y for x,y in zip(p,delta));hit=tuple(x+u*y for x,y in zip(at,ray))
 assert all(a*hit[0]+b*hit[1]<=c for a,b,c in cuts)
 data.update(exact_candidate_count=len(candidates),upper_envelope_pieces=pieces,maximum_signed_gap=dict(segment_t_exact=str(t),normal_parameter_exact=str(u),gap_squared_exact=str(u*u*L2),sign=(u>0)-(u<0),gap=enclosure(ar(u)*ar(L2).sqrt()),reference_point_exact=list(map(str,at)),polygon_intersection_exact=list(map(str,hit)),active_cuts=[k for qa,ra,k in upper if qa+ra*t==u],at_reference_segment_endpoint=t in[F(0),F(1)],at_normal_line_domain_endpoint=t in[lo,hi]))
 internal.update(peak_u=u,peak_t=t,peak_key=signed_key(u,L2));return internal

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',default='output/boundary_length_20261007/comparison/certified_polygon_comparison.json');ap.add_argument('--out',default='output/boundary_length_20261007/geometry_continuous_audit.json');ap.add_argument('--md',default='output/boundary_length_20261007/CONTINUOUS_GEOMETRY_AUDIT.md');ap.add_argument('--stations',type=int,default=401);args=ap.parse_args();ctx.prec=192;tic=time.time()
 source=ROOT/args.input;raw=source.read_bytes();source_hash=hashlib.sha256(raw).hexdigest();packet=json.loads(raw);assert packet['bounded_all_directions']and packet['status']=='COMPLETED_CUTS_AND_EXACT_POLYGON_CONSEQUENCES'
 manifest={str(source.relative_to(ROOT)):source_hash,str(Path(__file__).relative_to(ROOT)):sha(__file__)}
 for p,h in packet['source_sha256'].items():assert sha(ROOT/p)==h;manifest[p]=h
 cuts=[tuple(F(c[k])for k in['a_exact','b_exact','C_exact'])for c in packet['cuts']]
 poly=polygon_audit(cuts,packet['vertices_exact']);curves={}
 for side in['upper','lower']:
  path=ROOT/f'research/benchmarks_manuscript_dual_{side}.csv';manifest[str(path.relative_to(ROOT))]=sha(path)
  with path.open(newline='')as f:reader=csv.reader(f);next(reader);curves[side]=[tuple(F(x)for x in row)for row in reader]
  assert all(p[0]<q[0]for p,q in zip(curves[side],curves[side][1:]))
 points=sum(curves.values(),[]);origin=tuple(min(p[j]for p in points)for j in[0,1]);scale=tuple(max(p[j]for p in points)-origin[j]for j in[0,1]);assert all(s>0 for s in scale)
 branches={};allmodels=[];alllength=arb(0);allexcluded=arb(0);allmissing=arb(0);sample_maxima=[]
 for side,pts in curves.items():
  models=[segment(p,q,side,scale,cuts,j)for j,(p,q)in enumerate(zip(pts,pts[1:]))]
  lengths=[ar(m['L2']).sqrt()for m in models];total=sum(lengths,arb(0));excluded=sum((L*ar(1-size(m['inside']))for m,L in zip(models,lengths)),arb(0));missing=sum((L*ar(1-size(m['domain']))for m,L in zip(models,lengths)),arb(0))
  peak=max((m for m in models if m['domain']is not None),key=lambda m:m['peak_key']);cum=[arb(0)]
  for L in lengths:cum.append(cum[-1]+L)
  for j,m in enumerate(models):
   m['data']['normalized_metric_length']=enclosure(lengths[j]);m['data']['branch_arc_interval']=dict(start=enclosure(cum[j]/total),end=enclosure(cum[j+1]/total))
   if m['domain']is not None:m['data']['maximum_signed_gap']['branch_arc_fraction']=enclosure((cum[j]+lengths[j]*ar(m['peak_t']))/total)
   m['branch']=side
  # Recompute approximately equal-arclength stations, but evaluate each
  # selected rational t exactly and compare its algebraic gap to the peak.
  cf=[float(v.mid())for v in cum];lf=[float(v.mid())for v in lengths];samples=[];missed=0
  for k in range(args.stations):
   s=(k+.5)*cf[-1]/args.stations;j=min(len(models)-1,bisect.bisect_right(cf,s)-1);t=F((s-cf[j])/lf[j]);assert 0<=t<=1;m=models[j]
   if m['domain']is None or not(m['domain'][0]<=t<=m['domain'][1]):missed+=1;continue
   u=min(q+r*t for q,r,i in m['upper']);key=signed_key(u,m['L2']);assert key<=peak['peak_key'];samples.append((key,u,m['L2'],j,t))
  sample=max(samples)if samples else None
  gap=ar(max(F(0),peak['peak_u']))*ar(peak['L2']).sqrt();samplegap=ar(max(F(0),sample[1]))*ar(sample[2]).sqrt()if sample else None
  if sample:sample_maxima.append(sample)
  branches[side]=dict(observed_segments=len(models),observed_normalized_length=enclosure(total),excluded_reference_length=enclosure(excluded),excluded_reference_length_fraction=enclosure(excluded/total),normal_line_missing_length=enclosure(missing),normal_line_missing_length_fraction=enclosure(missing/total),continuous_max_positive_normal_gap=enclosure(gap),continuous_peak_segment=peak['data']['segment'],continuous_peak=peak['data']['maximum_signed_gap'],station_comparison=dict(stations=args.stations,missing_lines_at_stations=missed,maximum_positive_gap=None if samplegap is None else enclosure(samplegap),continuous_peak_minus_station_peak=None if samplegap is None else enclosure(gap-samplegap),every_station_below_continuous_peak_exactly=True,sampling='Approximately equal normalized-arclength midpoints; station parameters then evaluated as exact rationals. Station counts are never used as length coverage.'),segments=[m['data']for m in models])
  alllength+=total;allexcluded+=excluded;allmissing+=missing;allmodels+=models
 peak=max((m for m in allmodels if m['domain']is not None),key=lambda m:m['peak_key']);maxgap=ar(max(F(0),peak['peak_u']))*ar(peak['L2']).sqrt()
 report=dict(status='PASS_CONTINUOUS_DIGITIZED_GEOMETRY_AUDIT_NOT_NEW_PHYSICAL_CERTIFICATE',input_polygon=str(source.relative_to(ROOT)),cut_count=len(cuts),polygon_audit=poly,normalization=dict(origin_exact=list(map(str,origin)),axis_ranges_exact=list(map(str,scale)),definition='Ranges of observed upper/lower dual polylines; isolated right-tip marker and missing cap excluded.'),method='Exact rational normal-line projection and piecewise-affine upper-envelope breakpoints on every reference segment; exact algebraic comparison of squared signed gaps between segments. Arb192 encloses metric lengths and length fractions.',missing_cap='Unobserved upper/right-tip and lower/right-tip arcs have unknown length and are excluded, never connected or assigned zero gap. All percentages refer only to observed normalized branch length.',corner_rule='Each straight reference segment has its own constant outward normal. Both one-sided endpoint values are considered; no normal cone is swept across corners. Endpoints have zero arclength.',combined=dict(observed_normalized_length=enclosure(alllength),excluded_reference_length_fraction=enclosure(allexcluded/alllength),normal_line_missing_length_fraction=enclosure(allmissing/alllength),continuous_max_positive_normal_gap=enclosure(maxgap),peak_branch=peak['branch'],peak_segment=peak['data']['segment'],peak=peak['data']['maximum_signed_gap']),branches=branches,sha256_manifest=manifest,seconds=time.time()-tic,limitations='Digitized coordinates are exact inputs to this geometry calculation, not exact author bounds or rigorous digitization error bars. A reference point can be excluded while its normal-line outer endpoint is positive; exclusion and normal-line gap are reported separately.')
 assert sha(source)==source_hash,'Input polygon changed during audit; rerun on the new completed cuts'
 out=ROOT/args.out;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2))
 fmt=lambda v:f"{v['midpoint_float']:.9g}"
 lines=['# Continuous normalized-boundary geometry audit','',f"Input: `{args.input}`; SHA256 `{source_hash}`.",'','This audits every point of each observed straight reference segment, including between the plotted stations. It produces no new physical inequality.','', '| Observed branch | Maximum positive normal gap | Excluded reference length | Normal lines missing polygon |', '|---|---:|---:|---:|']
 for side,d in branches.items():lines.append(f"| {side} | {fmt(d['continuous_max_positive_normal_gap'])} | {100*d['excluded_reference_length_fraction']['midpoint_float']:.6f}% | {100*d['normal_line_missing_length_fraction']['midpoint_float']:.6f}% |")
 c=report['combined'];lines+=['',f"Combined excluded observed length: **{100*c['excluded_reference_length_fraction']['midpoint_float']:.6f}%**. The global positive normal-gap maximum is **{fmt(c['continuous_max_positive_normal_gap'])}**, on the {c['peak_branch']} branch, segment {c['peak_segment']}.",'','The axes are divided by their observed dual ranges. For each segment, an unnormalized metric normal makes every halfplane intersection endpoint affine in the segment parameter. The script clips the exact domain where the normal line meets the polygon, evaluates endpoints and all active affine-envelope crossings, and checks the complete resulting piecewise envelope. Signed squared gaps allow exact comparisons across segments; square-root arclength totals and fractions have directed Arb enclosures in the JSON.','',f"All {args.stations} tested midpoint stations per branch lie below their continuous branch maximum. Their sampling deficit is recorded separately; coverage percentages use clipped segment lengths, never station or normal-angle counts.",'','The missing right cap is not reconstructed. Percentages cover only the observed digitized polylines. At corners the two adjacent segment normals are treated separately, with no normal sweep. A missing normal-line intersection is explicitly recorded and is not silently converted to a gap value. A point can be excluded by one cut while its outward normal-line endpoint remains positive, so exclusion and positive gap are distinct diagnostics.','',f"Reproduce: `python research/boundary_length_20261007/geometry_continuous_audit.py --input {args.input} --out {args.out} --md {args.md}`."]
 (ROOT/args.md).write_text('\n'.join(lines)+'\n');print(json.dumps(dict(combined=c,branch_summaries={k:{q:v[q]for q in['continuous_max_positive_normal_gap','excluded_reference_length_fraction','normal_line_missing_length_fraction','station_comparison']}for k,v in branches.items()},seconds=report['seconds']),indent=2))
if __name__=='__main__':main()
