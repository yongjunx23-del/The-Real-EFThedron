"""Normalize existing local normal gaps: dual=0, primal=1.

All ratios use exact rational line intersections with the supplied reference
coordinates. No ratios are clamped into [0,1], and no missing points are filled.
"""
from pathlib import Path
from fractions import Fraction as F
import csv, hashlib, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'output/primal_dual_fraction_20261008'
SOURCE = ROOT/'output/certified_points_ccw_20261008/plot_data.json'
GEO = ROOT/'output/six_targets_20261007/comparison/geometry_continuous_audit.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = json.loads(SOURCE.read_text(encoding='utf-8'))
    geo = json.loads(GEO.read_text(encoding='utf-8'))
    paths = [ROOT/f'research/benchmarks_manuscript_primal_{s}.csv' for s in ['lower','upper']]
    primal = {}
    for side,path in zip(['lower','upper'],paths):
        with path.open(newline='') as f:
            reader=csv.reader(f);next(reader)
            primal[side]=[tuple(map(F,row)) for row in reader]

    rows, omitted = [], []
    for marker in data['selected_markers']:
        p = tuple(map(F,marker['reference_point_exact']))
        ray = tuple(map(F,marker['physical_ray_unnormalized_exact']))
        uc = F(marker['normal_parameter_exact'])
        a0,a2 = map(F,marker['objective_exact'])
        C = F(marker['support_upper_exact'])
        assert a0*(p[0]+uc*ray[0])+a2*(p[1]+uc*ray[1]) == C
        hits=[]
        for j,(a,b) in enumerate(zip(primal[marker['branch']],primal[marker['branch']][1:])):
            edge=sub(b,a);offset=sub(a,p);det=cross(ray,edge)
            if det == 0:
                # No silently discarded collinear overlap at any displayed station.
                assert cross(offset,ray) != 0
                continue
            up=cross(offset,edge)/det;v=cross(offset,ray)/det
            if 0 <= v <= 1:
                hit=tuple(p[k]+up*ray[k] for k in range(2))
                assert hit==tuple(a[k]+v*edge[k] for k in range(2))
                hits.append((up,j,v,hit))
        if not hits:
            omitted.append(dict(marker, reason='No normal-line intersection with the corresponding primal branch; ratio undefined.'))
            continue
        up,j,v,hit=max(hits)
        if up >= 0:
            omitted.append(dict(marker, reason='Zero or reversed primal-dual separation; the displayed inner-to-outer ordering is undefined.',primal_parameter_exact=str(up)))
            continue
        eta=uc/up
        L2=F(geo['branches'][marker['branch']]['segments'][marker['segment']]['metric_length_squared_exact'])
        delta_primal=float(up)*math.sqrt(float(L2))
        assert abs(float(eta)-marker['gap']/delta_primal)<2e-14
        rows.append(dict(marker, s_over_L=marker['ccw_percent']/100,
            eta_exact=str(eta), eta=float(eta), primal_parameter_exact=str(up),
            primal_normal_distance=delta_primal, primal_segment=j,
            primal_segment_fraction_exact=str(v), primal_intersection_exact=list(map(str,hit)),
            primal_intersection_count=len(hits), exact_ratio_unclipped=True))

    assert len(rows)+len(omitted)==19 and len(rows)==15 and len(omitted)==4
    assert all(0<F(r['eta_exact'])<1 for r in rows)
    rows.sort(key=lambda r:r['s_over_L'])
    split=data['branch_join_percent']/100
    plt.rcParams.update({
        'font.family':'serif','font.serif':['STIXGeneral'],'mathtext.fontset':'stix',
        'font.size':9,'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,
        'legend.fontsize':9,'axes.linewidth':.65,'xtick.direction':'in','ytick.direction':'in',
        'xtick.major.size':3.2,'ytick.major.size':3.2,'xtick.minor.size':1.8,'ytick.minor.size':1.8,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(7.05,3.65))
    fig.subplots_adjust(left=.110,right=.984,bottom=.160,top=.890)
    ax.axhspan(-.18,0,color='#f4e5e3',zorder=0)
    ax.axhspan(1,1.18,color='#e5efe8',zorder=0)
    dual=ax.axhline(0,color='#252525',linewidth=.9,label='Dual',zorder=2)
    primal_line=ax.axhline(1,color='#3d6a9c',linewidth=1.1,linestyle=(0,(3,2)),label='Primal',zorder=2)
    pts=ax.scatter([r['s_over_L'] for r in rows],[r['eta'] for r in rows],
        s=18,color='#00785b',marker='o',edgecolors='white',linewidths=.3,label='Certified',zorder=4)
    ax.axvline(split,color='#a0a0a0',linewidth=.6,linestyle=(0,(3,3)),zorder=1)
    ax.text(.25,-.090,'Excluded',ha='center',va='center',fontsize=10,color='#743d38')
    ax.text(.25,1.090,'Allowed (primal side)',ha='center',va='center',fontsize=10,color='#315b41')
    for x,label,ha in [(0,r'$g_{0,\min}$','left'),(split,r'$g_{0,\max}^{\dagger}$','center'),(1,r'$g_{0,\min}$','right')]:
        ax.text(x,1.025,label,transform=ax.get_xaxis_transform(),ha=ha,va='bottom',fontsize=10)
    ax.text(.25,.94,'Lower',ha='center',va='center',fontsize=8.5,color='#686868')
    ax.text(.75,.94,'Upper',ha='center',va='center',fontsize=8.5,color='#686868')
    ax.set(xlim=(0,1),ylim=(1.18,-.18),xlabel=r'Boundary arclength, $s/L$',
        ylabel=r'Fraction of primal–dual gap, $\eta$')
    ax.set_xticks(np.linspace(0,1,6));ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_yticks([0,.25,.5,.75,1]);ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2));ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='both',top=True,right=True,width=.6)
    ax.legend(handles=[pts,dual,primal_line],loc='center right',bbox_to_anchor=(.96,.48),
        frameon=False,handlelength=2.2,labelspacing=.4,borderpad=.15)
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'boundary_gap_fraction.{ext}',**({'dpi':400} if ext=='png' else {}))
    plt.close(fig)

    caption=r'''Certified bounds expressed as a fraction of the local primal–dual separation. At each displayed boundary position, define $\eta(s)=\delta_C(s)/\delta_P(s)$, where $\delta_C$ and $\delta_P$ are signed displacements from the dual reference along its outward normal in axis-range-normalized coordinates. Thus the dual reference is $\eta=0$ and the primal reference is $\eta=1$. The vertical axis is inverted: the excluded side lies above the dual, and the primal-allowed side below the primal. Increasing $\eta$ measures the fraction of the reference gap closed by a certified bound; the interval between the two references is not established as allowed by this plot. All 15 defined certificate ratios lie strictly between zero and one. Four locations without a primal normal-line intersection are omitted; ratios are not clipped. Points retain the previous display thinning. Arclength increases counterclockwise from minimum $g_0$. The dagger denotes the join across the unresolved right-tip arc, excluded from $L$.
'''
    (OUT/'caption.txt').write_text(caption,encoding='utf-8')
    explanation='''LOCAL VERTICAL NORMALIZATION

Let delta_P(s)<0 be the signed normal displacement from the dual reference to the outermost same-branch primal intersection, and delta_C(s) the displacement to the displayed certified supporting line. Both use the same reference point, normal line and axis-range normalization.

    eta(s) = delta_C(s)/delta_P(s).

This is a local normalization, not a single constant rescaling of the old vertical axis. The dual is exactly0 and the primal exactly1 wherever the denominator exists and is nonzero. The plotted axis runs downward from0 to1, so the excluded side is visually above the dual and the primal-allowed side below the primal. Those labels describe the local sides of the reference boundaries, not an assertion that an entire infinite ray beyond the primal remains feasible.

For a point p and the saved unnormalized ray r, the certified inequality a.g<=C gives u_C=(C-a.p)/(a.r). Exact intersections p+u_P*r with each observed same-branch primal segment are computed using rational determinants. The largest valid u_P is used. The plotted ratio is exactly u_C/u_P; the metric normalization factor cancels. No near-zero denominator regularization, clipping into[0,1], invented endpoint, or shift toward the middle is used.

The horizontal reference lines indicate the normalization levels. They do not fill missing primal intersections. Four previously shown certificate locations have no such intersection and have no ratio. Their source data and omission reasons remain in plot_data.json. The source curves come from the published figure; the rational arithmetic is exact for those input coordinates, not an uncertainty enclosure for the original authors' data. No new physical bound or certificate is claimed.
'''
    (OUT/'normalization.txt').write_text(explanation,encoding='utf-8')
    report=dict(status='EXACT_LOCAL_GAP_RATIO_DISPLAY_OF_EXISTING_CERTIFICATES',
        displayed=rows,undefined=omitted,displayed_count=len(rows),undefined_count=len(omitted),
        minimum_ratio=min(r['eta'] for r in rows),maximum_ratio=max(r['eta'] for r in rows),
        all_defined_ratios_strictly_between_zero_and_one=True,clipped_ratios=0,
        vertical_axis='eta=delta_certified/delta_primal; inverted display:dual0 above primal1',
        inherited_thinning=True,new_certifications=0,
        sha256_manifest={p.relative_to(ROOT).as_posix():sha(p)for p in [SOURCE,GEO,*paths,Path(__file__)]},
        output_sha256={p.name:sha(p)for p in [OUT/'boundary_gap_fraction.pdf',OUT/'boundary_gap_fraction.svg',OUT/'boundary_gap_fraction.png',OUT/'caption.txt',OUT/'normalization.txt']})
    (OUT/'plot_data.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    with (OUT/'gap_fractions.csv').open('w',newline='',encoding='utf-8')as f:
        columns=['s_over_L','branch','percent','cut_index','eta','eta_exact','gap','primal_normal_distance','source']
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print(json.dumps({k:report[k]for k in ['displayed_count','undefined_count','minimum_ratio','maximum_ratio','all_defined_ratios_strictly_between_zero_and_one']},indent=2))


if __name__=='__main__':
    main()
