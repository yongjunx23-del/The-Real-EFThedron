"""Graphical support functions for the massive 4D (g0,g2) comparison.

Convention: h(theta)=sup(g0*cos(theta)+g2*sin(theta)); angles default to
DEGREES. No rescaling of either coupling is applied. At 0,90,180,270 degrees
the functions give max g0, max g2, minus min g0, and minus min g2.

Sources are PDF vector outlines, not published optimizing dual witnesses.
Taking max over vertices exactly computes the support of their convex hull
without inventing an interpolated or closing boundary arc. The published
dual polylines end near g0=2.6563, BEFORE the maximum-g0 tip. We retain the
separate quoted 2021 maximum-g0 coordinate 2.7272 at the approximate plotted
ordinate .06. This is a graphical anchor, NOT an attained primal point.

support_report() exposes the omitted-cap sensitivity and an explicit working
graphical uncertainty envelope. That envelope is not a statistical confidence
interval, numerical certificate, or bound on the unpublished optimizer data.
No function here proves that a newly optimized witness is continuum-valid.
"""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import argparse
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
EXTRACTION_TOL=np.array([5e-6,5e-6])
MARKER=np.array([2.7272,.06])
# The graphical marker x differs from the quoted numerical coordinate by
# 0.0005644; .001 in y is an explicitly chosen sensitivity allowance.
MARKER_WORKING_TOL=np.array([.0006,.001])
ALIASES={'dual':'published_dual','blue':'manuscript','published_primal':'primal'}


def _read(name):
    return np.loadtxt(ROOT/name,delimiter=',',skiprows=1)


@lru_cache(None)
def _data(benchmark):
    benchmark=ALIASES.get(benchmark,benchmark)
    if benchmark=='published_dual':
        files=['benchmarks_manuscript_dual_upper.csv','benchmarks_manuscript_dual_lower.csv']
    elif benchmark=='primal':
        files=['benchmarks_manuscript_primal_upper.csv','benchmarks_manuscript_primal_lower.csv']
    elif benchmark=='manuscript':
        files=['benchmarks_manuscript_blue.csv']
    else:
        raise ValueError('benchmark must be published_dual, manuscript, or primal')
    blocks=[_read(name) for name in files]
    labels=[f'{name}:vector_vertex_{i}' for name,pts in zip(files,blocks) for i in range(len(pts))]
    pts=np.vstack(blocks+[np.array([[0.,0.]])])
    labels.append('exact_free_point')
    if benchmark=='published_dual':
        pts=np.vstack([pts,MARKER]);labels.append('isolated_maxg0_graphical_anchor')
    return pts,np.array(labels),files


def _directions(theta,degrees):
    theta=np.asarray(theta,dtype=float)
    if np.any(~np.isfinite(theta)): raise ValueError('Angles must be finite')
    radians=np.deg2rad(theta) if degrees else theta
    directions=np.column_stack([np.cos(radians).ravel(),np.sin(radians).ravel()])
    directions[np.abs(directions)<1e-15]=0.
    return theta,directions


def _restore(values,shape):
    result=np.asarray(values).reshape(shape)
    return result.item() if shape==() else result


def support(theta,benchmark='published_dual',*,degrees=True):
    """Return the exact support of the extracted vertices plus stated anchors.

    This is a graphical benchmark estimate, NOT the support of a newly
    reconstructed continuum dual region. Supports scalar or array angles.
    """
    theta,directions=_directions(theta,degrees)
    points,_,_=_data(benchmark)
    return _restore(np.max(directions@points.T,axis=1),theta.shape)


def support_report(theta,benchmark='published_dual',*,degrees=True):
    """Support estimate, active source point, and graphical sensitivity data.

    For the missing right cap, the conditional high envelope puts an unknown
    point anywhere in the rectangle x=[last sampled x,2.7272],
    y=[lower endpoint y,upper endpoint y]. This assumes the omitted cap stays
    between those endpoint ordinates, as expected for the convex upper/lower
    boundaries. It is deliberately labeled conditional, not certified.
    Numerical optimizer and continuum errors of the source are UNKNOWN.
    """
    benchmark=ALIASES.get(benchmark,benchmark)
    theta,directions=_directions(theta,degrees)
    points,labels,files=_data(benchmark)
    products=directions@points.T
    indexes=np.argmax(products,axis=1)
    values=products[np.arange(len(directions)),indexes]
    reading_tol=np.abs(directions)@EXTRACTION_TOL
    low=values-reading_tol;high=values+reading_tol
    cap_effect=np.zeros_like(values)
    marker_dominates=np.zeros(len(values),dtype=bool)
    cap_sensitive=np.zeros(len(values),dtype=bool)
    without_marker=values.copy()
    if benchmark=='published_dual':
        without_marker=np.max(products[:,:-1],axis=1)
        marker_h=directions@MARKER
        marker_tol=np.abs(directions)@MARKER_WORKING_TOL
        marker_dominates=indexes==len(points)-1
        low=np.maximum(without_marker-reading_tol,marker_h-marker_tol)
        high=np.maximum(without_marker+reading_tol,marker_h+marker_tol)
        upper=_read('benchmarks_manuscript_dual_upper.csv')
        lower=_read('benchmarks_manuscript_dual_lower.csv')
        up=upper[np.argmax(upper[:,0])];lo=lower[np.argmax(lower[:,0])]
        # No arc is added to the vertex set. This box is uncertainty metadata.
        xmin=min(up[0],lo[0]);xmax=MARKER[0]
        ymin=min(up[1],lo[1],MARKER[1]);ymax=max(up[1],lo[1],MARKER[1])
        cap_h=(np.where(directions[:,0]>=0,xmax,xmin)*directions[:,0]
               +np.where(directions[:,1]>=0,ymax,ymin)*directions[:,1])
        cap_effect=np.maximum(cap_h-values,0)
        cap_sensitive=cap_effect>reading_tol
        high=np.maximum(high,cap_h+reading_tol)
    result={
        'benchmark':benchmark,
        'theta':theta.item() if theta.shape==() else theta.copy(),
        'angle_unit':'degrees' if degrees else 'radians',
        'support':_restore(values,theta.shape),
        'active_g0':_restore(points[indexes,0],theta.shape),
        'active_g2':_restore(points[indexes,1],theta.shape),
        'active_source':_restore(labels[indexes],theta.shape),
        'support_without_isolated_marker':_restore(without_marker,theta.shape),
        'isolated_marker_dominates':_restore(marker_dominates,theta.shape),
        'conditional_missing_cap_sensitive':_restore(cap_sensitive,theta.shape),
        'conditional_missing_cap_support_allowance':_restore(cap_effect,theta.shape),
        'extraction_working_halfwidth':_restore(reading_tol,theta.shape),
        'graphical_low':_restore(low,theta.shape),
        'graphical_high_conditional':_restore(high,theta.shape),
        'source_files':files,
        'warning':'Graphical support estimate only; source optimizer/continuum uncertainty unavailable. Right-cap envelope is conditional, not a certified interval.'
    }
    return result


def compare_support(theta,new_support,*,degrees=True):
    """Compare independently computed supports in EXACTLY the same direction.

    delta_vs_dual=new-published: negative means tighter numerically.
    gap_closed_fraction=(published-new)/(published-primal) is only supplied
    when the plotted primal-dual gap exceeds graphical uncertainty. Negative
    fractions mean weaker than the published reference. A fraction >1 flags
    crossing below the plotted primal support and calls for a validity audit.
    """
    theta=np.asarray(theta,dtype=float)
    new=np.broadcast_to(np.asarray(new_support,dtype=float),theta.shape)
    dual=support_report(theta,'published_dual',degrees=degrees)
    hm=np.asarray(support(theta,'manuscript',degrees=degrees))
    hp=np.asarray(support(theta,'primal',degrees=degrees))
    hd=np.asarray(dual['support'])
    width=np.asarray(dual['graphical_high_conditional'])-np.asarray(dual['graphical_low'])
    gap=hd-hp
    fraction=np.full(theta.shape,np.nan)
    np.divide(hd-new,gap,out=fraction,where=gap>np.maximum(width,1e-8))
    result={
        'theta':theta.item() if theta.shape==() else theta.copy(),
        'angle_unit':'degrees' if degrees else 'radians',
        'new_support':_restore(new,theta.shape),
        'published_dual_support':dual['support'],
        'manuscript_support':_restore(hm,theta.shape),
        'primal_support':_restore(hp,theta.shape),
        'delta_vs_dual':_restore(new-hd,theta.shape),
        'improvement_vs_dual':_restore(hd-new,theta.shape),
        'delta_vs_manuscript':_restore(new-hm,theta.shape),
        'improvement_vs_manuscript':_restore(hm-new,theta.shape),
        'gap_closed_fraction':_restore(fraction,theta.shape),
        'delta_vs_dual_graphical_low':_restore(new-np.asarray(dual['graphical_high_conditional']),theta.shape),
        'delta_vs_dual_graphical_high':_restore(new-np.asarray(dual['graphical_low']),theta.shape),
        'conditional_missing_cap_sensitive':dual['conditional_missing_cap_sensitive'],
        'isolated_marker_dominates':dual['isolated_marker_dominates'],
        'below_plotted_primal':_restore(new<hp-np.maximum(width,1e-8),theta.shape),
        'warning':'Same-angle graphical comparison, not certification of new witnesses or numerical source bounds.'
    }
    return result


def provenance():
    discrepancies={}
    for name in ['dual_upper','dual_lower','primal_upper','primal_lower']:
        discrepancies[name]=np.max(abs(_read(f'benchmarks_{name}.csv')
                    -_read(f'benchmarks_manuscript_{name}.csv')),axis=0).tolist()
    return {
        'definition':'h(theta)=max(g0*cos(theta)+g2*sin(theta)); theta degrees by default',
        'coupling_normalization':'M/(32*pi)=g0+g2*(sbar^2+tbar^2+ubar^2)+...; m^2=1',
        'published_2026_url':'https://arxiv.org/html/2605.06613v1',
        'published_vector_source':'https://arxiv.org/src/2605.06613v1 ; member img/fig1_final-crop.pdf',
        'quoted_maxg0_source':'https://arxiv.org/pdf/2106.10257 ; Appendix G quotes 2.7272',
        'isolated_marker':MARKER.tolist(),
        'isolated_marker_working_tolerance':MARKER_WORKING_TOL.tolist(),
        'extraction_working_tolerance':EXTRACTION_TOL.tolist(),
        'independent_source_copy_coordinate_discrepancies':discrepancies,
        'manuscript_file':'D:/Dropbox/Codex/FULL U/The_Real_EFThedron-7.pdf ; Fig.4 p.6',
        'right_cap_note':'The outer published polylines end near g0=2.6563. A separate max-g0 marker is included. No closing/interpolating arc is invented. Missing-cap sensitivity is reported separately.',
        'uncertainty_note':'Working tolerances quantify graphical sensitivity only. They are not statistical errors or certified optimizer/continuum error bounds.',
        'restriction':'No assumption of saturation of the published exclusion curve; graphical markers are reference coordinates, not physical amplitudes.'
    }


def _json_ready(obj):
    if isinstance(obj,np.ndarray): return _json_ready(obj.tolist())
    if isinstance(obj,np.generic): return _json_ready(obj.item())
    if isinstance(obj,dict): return {k:_json_ready(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)): return [_json_ready(v) for v in obj]
    if isinstance(obj,float) and not np.isfinite(obj): return None
    return obj


def _self_test():
    angles=np.array([0.,90.,180.,270.,360.])
    hs=support(angles)
    np.testing.assert_allclose(hs,[2.7272,1.0816846732472978,8.014557983981451,0.,2.7272],atol=1e-9)
    np.testing.assert_allclose(support(np.deg2rad(angles),degrees=False),hs)
    np.testing.assert_allclose(support(angles+720),hs,atol=1e-12)
    theta=np.linspace(0,360,1441)
    dual=support(theta);primal=support(theta,'primal')
    assert np.min(dual-primal)>-1e-5
    details=support_report(theta)
    assert np.all(np.asarray(details['graphical_low'])<=dual+1e-12)
    assert np.all(np.asarray(details['graphical_high_conditional'])>=dual-1e-12)
    comparison=compare_support(angles,support(angles,'manuscript'))
    np.testing.assert_allclose(comparison['gap_closed_fraction'][1:3],[.1558807,.2060588],atol=1e-4)
    return {'cardinal_supports':hs.tolist(),'minimum_graphical_primal_dual_gap':float(np.min(dual-primal)),
            'largest_conditional_missing_cap_allowance':float(np.max(details['conditional_missing_cap_support_allowance'])),
            'manuscript_same_angle_comparison':comparison}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',help='Optional JSON metadata/check report destination')
    args=parser.parse_args()
    result={'provenance':provenance(),'checks':_self_test()}
    encoded=json.dumps(_json_ready(result),indent=2,allow_nan=False)
    if args.out: Path(args.out).write_text(encoded,encoding='utf8')
    print(encoded)
