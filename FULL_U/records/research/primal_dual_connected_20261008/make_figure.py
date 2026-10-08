"""Display all existing boundary markers, with primal=0 and dual=1.

Connecting segments are visual guides, not new pointwise certificates.
Undefined or reversed reference gaps get a separate position-only strip.
"""
from pathlib import Path
from fractions import Fraction as F
import csv
import hashlib
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'output/primal_dual_connected_20261008'
SOURCE = ROOT/'output/certified_points_plot_20261007/plot_data.json'
GEO = ROOT/'output/six_targets_20261007/comparison/geometry_continuous_audit.json'
PREVIOUS = ROOT/'output/primal_dual_fraction_20261008/plot_data.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = json.loads(SOURCE.read_text(encoding='utf-8'))
    geo = json.loads(GEO.read_text(encoding='utf-8'))
    previous = json.loads(PREVIOUS.read_text(encoding='utf-8'))
    lengths = {s: geo['branches'][s]['observed_normalized_length']['midpoint_float']
               for s in ['lower', 'upper']}
    total = sum(lengths.values())
    split = lengths['lower']/total
    paths = [ROOT/f'research/benchmarks_manuscript_primal_{s}.csv'
             for s in ['lower', 'upper']]
    primal = {}
    for side, path in zip(['lower', 'upper'], paths):
        with path.open(newline='') as f:
            reader = csv.reader(f)
            next(reader)
            primal[side] = [tuple(map(F, row)) for row in reader]

    rows, undefined = [], []
    for marker in data['markers']:
        fraction = marker['percent']/100
        arc = (fraction*lengths['lower'] if marker['branch'] == 'lower'
               else lengths['lower']+(1-fraction)*lengths['upper'])/total
        row = dict(marker, s_over_L=arc)
        p = tuple(map(F, marker['reference_point_exact']))
        ray = tuple(map(F, marker['physical_ray_unnormalized_exact']))
        uc = F(marker['normal_parameter_exact'])
        a0, a2 = map(F, marker['objective_exact'])
        C = F(marker['support_upper_exact'])
        assert a0*(p[0]+uc*ray[0])+a2*(p[1]+uc*ray[1]) == C
        hits = []
        for j, (a, b) in enumerate(zip(primal[marker['branch']], primal[marker['branch']][1:])):
            edge = sub(b, a)
            offset = sub(a, p)
            det = cross(ray, edge)
            if det == 0:
                assert cross(offset, ray) != 0
                continue
            up = cross(offset, edge)/det
            v = cross(offset, ray)/det
            if 0 <= v <= 1:
                hit = tuple(p[k]+up*ray[k] for k in range(2))
                assert hit == tuple(a[k]+v*edge[k] for k in range(2))
                hits.append((up, j, v, hit))
        if not hits:
            undefined.append(dict(row, reason='No corresponding primal normal-line intersection.',
                                  ratio=None, display='Position-only strip'))
            continue
        up, j, v, hit = max(hits)
        L2 = F(geo['branches'][marker['branch']]['segments'][marker['segment']]['metric_length_squared_exact'])
        dp = float(up)*math.sqrt(float(L2))
        if up >= 0:
            undefined.append(dict(row, reason='Reference order reversed or coincident; no positive inward primal-dual separation.',
                                  ratio=None, primal_parameter_exact=str(up),
                                  primal_normal_distance=dp, display='Position-only strip'))
            continue
        ratio = 1-uc/up
        assert abs(float(ratio)-(1-marker['gap']/dp)) < 2e-14
        rows.append(dict(row, ratio_exact=str(ratio), ratio=float(ratio),
                         primal_parameter_exact=str(up), primal_normal_distance=dp,
                         primal_segment=j, primal_segment_fraction_exact=str(v),
                         primal_intersection_exact=list(map(str, hit)),
                         primal_intersection_count=len(hits), clipped=False))

    rows.sort(key=lambda r: r['s_over_L'])
    undefined.sort(key=lambda r: r['s_over_L'])
    assert len(rows) == 27 and len(undefined) == 5
    assert len({r['cut_index'] for r in rows+undefined}) == len(data['markers']) == 32
    assert sum(r['ratio'] > 1 for r in rows) == 9
    assert all(not r['analytic'] for r in rows)
    assert sum(r['analytic'] for r in undefined) == 1
    by_cut = {r['cut_index']: r for r in rows}
    for old in previous['displayed']:
        assert F(by_cut[old['cut_index']]['ratio_exact']) == 1-F(old['eta_exact'])
        assert abs(by_cut[old['cut_index']]['s_over_L']-old['s_over_L']) < 2e-15

    # A missing station or the omitted right-tip arc breaks a visual guide.
    # Straight segments interpolate only between actual displayed markers.
    groups = []
    for side in ['lower', 'upper']:
        group = []
        for row in sorted([r for r in rows+undefined if r['branch'] == side],
                          key=lambda r: r['s_over_L']):
            if row.get('ratio') is None:
                if group:
                    groups.append(group)
                    group = []
            else:
                group.append(row)
        if group:
            groups.append(group)

    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['STIXGeneral'], 'mathtext.fontset': 'stix',
        'font.size': 9, 'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
        'legend.fontsize': 9, 'axes.linewidth': .65,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.size': 3.2, 'ytick.major.size': 3.2,
        'xtick.minor.size': 1.8, 'ytick.minor.size': 1.8,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(7.05, 4.30))
    ax = fig.add_axes([.110, .235, .874, .650])
    strip = fig.add_axes([.110, .143, .874, .064], sharex=ax)
    green = '#00785b'

    def draw_points(axis, inset=False):
        for group in groups:
            axis.plot([r['s_over_L'] for r in group], [r['ratio'] for r in group],
                      color=green, lw=.70, ls=(0, (2.4, 2)), alpha=.72, zorder=3)
        axis.scatter([r['s_over_L'] for r in rows], [r['ratio'] for r in rows],
                     s=11 if inset else 12, color=green, marker='o',
                     edgecolors='white', linewidths=.3, zorder=4)

    ax.axhspan(1, 2.13, color='#f4e5e3', zorder=0)
    ax.axhspan(-.17, 0, color='#e5efe8', zorder=0)
    dual = ax.axhline(1, color='#252525', linewidth=.9, label='Dual', zorder=2)
    primal_line = ax.axhline(0, color='#3d6a9c', linewidth=1.1,
                           linestyle=(0, (3, 2)), label='Primal', zorder=2)
    draw_points(ax)
    ax.axvline(split, color='#a0a0a0', linewidth=.6, linestyle=(0, (3, 3)), zorder=1)
    ax.text(.17, 1.87, 'Excluded', ha='center', va='center', fontsize=10, color='#743d38')
    ax.text(.25, -.086, 'Allowed (primal side)', ha='center', va='center', fontsize=9, color='#315b41')
    for x, label, ha in [(0, r'$g_{0,\min}$', 'left'),
                         (split, r'$g_{0,\max}^{\dagger}$', 'center'),
                         (1, r'$g_{0,\min}$', 'right')]:
        ax.text(x, 1.025, label, transform=ax.get_xaxis_transform(), ha=ha,
                va='bottom', fontsize=10)
    ax.text(.27, .14, 'Lower', ha='center', va='center', fontsize=8.5, color='#686868')
    ax.text(.75, .14, 'Upper', ha='center', va='center', fontsize=8.5, color='#686868')
    ax.set(xlim=(0, 1), ylim=(-.17, 2.13), ylabel=r'Normalized bound, $q$')
    ax.set_yticks([0, .5, 1, 1.5, 2])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='both', top=True, right=True, width=.6, labelbottom=False)
    handle = Line2D([], [], color=green, lw=.70, ls=(0, (2.4, 2)), marker='o',
                    markersize=3.5, markeredgewidth=.3, markeredgecolor='white', label='Certified')
    ax.legend(handles=[handle, dual, primal_line], loc='center left',
              bbox_to_anchor=(.055, .205), ncol=3, frameon=False,
              handlelength=1.7, columnspacing=1.3, borderpad=.15)

    # Magnification adds space without jittering, averaging, or dropping dots.
    inset = ax.inset_axes([.650, .625, .320, .330])
    inset.axhspan(1, 2.05, color='#f4e5e3', zorder=0)
    inset.axhline(1, color='#252525', linewidth=.7, zorder=2)
    draw_points(inset, inset=True)
    inset.axvline(split, color='#a0a0a0', linewidth=.5, linestyle=(0, (3, 3)))
    inset.set(xlim=(.413, .535), ylim=(.63, 2.02), xticks=[.42, .47, .52], yticks=[1, 1.5, 2])
    inset.tick_params(which='both', labelsize=7, top=True, right=True, width=.5, length=2)
    inset.set_facecolor('white')

    # All five non-normalizable source markers remain visible at their true x.
    for r in undefined:
        strip.scatter([r['s_over_L']], [.5], marker='D' if r['analytic'] else 'o',
                      s=14, facecolors='white', edgecolors='#777777' if r['analytic'] else green,
                      linewidths=.65, zorder=4)
    strip.text(.67, .53, r'$q$ undefined', ha='center', va='center', fontsize=8, color='#666666')
    strip.set(ylim=(0, 1), yticks=[], xlabel=r'Boundary arclength, $s/L$')
    strip.set_xticks(np.linspace(0, 1, 6))
    strip.xaxis.set_minor_locator(AutoMinorLocator(2))
    strip.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    strip.spines[['top', 'left', 'right']].set_visible(False)
    strip.spines['bottom'].set_color('#aaaaaa')
    strip.tick_params(which='both', width=.6)
    strip.xaxis.labelpad = 5

    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT/f'boundary_all_points.{ext}', **({'dpi': 400} if ext == 'png' else {}))
    plt.close(fig)

    caption = r'''Certified bounds around the observed boundary, parameterized counterclockwise from minimum $g_0$. The ordinate is $q(s)=1-\delta_C(s)/\delta_P(s)$, with signed outward normal displacements from the dual reference measured in axis-range-normalized coordinates. Thus primal is $q=0$ and dual is $q=1$; smaller $q$ means a tighter outer bound. Values above one are weaker than the published dual. The excluded and allowed labels indicate the local outer and primal-inner sides, respectively. All 31 existing numerical certificate markers and the analytic $g_2\geq0$ marker are retained: 27 numerical markers have meaningful normalized ordinates, including 9 above the dual. The position-only strip shows four numerical markers without a primal intersection (open circles) and one analytic marker at a reversed reference separation (open diamond). Dashed segments merely connect neighboring displayed markers as visual guides; they are not certified interpolants. Connections stop at undefined stations and across the unresolved right-tip arc. The inset enlarges the crowded lower-boundary portion; no dots are averaged, shifted, or removed. The dagger denotes the join across the unresolved maximum-$g_0$ tip, whose length is excluded from $L$.
'''
    (OUT/'caption.txt').write_text(caption, encoding='utf-8')
    notes = '''DEFINITION AND DISPLAY CHECKS

q = 1 - delta_certified / delta_primal = 1 - u_certified/u_primal.
Distances are measured on the same outward normal line of the dual reference, after scaling each physical axis by its displayed dual range. The exact rational reference coordinates and normal rays are unchanged. Both reference lines are horizontal: primal0, dual1. This undoes the prior inverted axis while changing eta to1-eta. Smaller q means tighter. Values above1 are retained.

All32 markers from the original complete marker inventory are accounted for, with no sign or density filters:31 numerical inequalities and1 analytic inequality. The figure includes27 numerical ratios; the other4 numerical positions and the analytic position are shown in the separate undefined strip. The analytical point has positive outward primal separation in the supplied reference curves, so treating its negative formal q as a bound inside the primal would misrepresent the reference ordering. No missing or reversed denominator is repaired, extrapolated, or clipped.

The inset repeats the crowded portion at larger scale. It adds no data. The connecting dashed lines are piecewise-linear visual guides through numerical marker coordinates. They do not claim that interpolation is a certified boundary. Connections break at every undefined marker and the omitted maximum-g0 cap. Open circles in the lower strip denote missing primal intersections; the open diamond denotes the analytic constraint at the reversed local reference gap.

Each ratio is computed from exact rational segment intersections with the saved primal reference coordinates. Each own-cut equality is checked exactly. For the15 previously shown ratios, q=1-eta is checked exactly, and arclength is checked unchanged to floating-point rounding. Exact arithmetic here does not enclose uncertainty in the published reference coordinates. No optimization or physical certification was performed for this display revision. Previously sealed science data and earlier figures remain unchanged.
'''
    (OUT/'normalization.txt').write_text(notes, encoding='utf-8')
    report = dict(status='ALL_EXISTING_MARKERS_DISPLAY_WITH_PRIMAL_ZERO_DUAL_ONE',
                  displayed=rows, undefined=undefined, total_source_markers=len(data['markers']),
                  normalized_count=len(rows), undefined_count=len(undefined),
                  weak_count=sum(r['ratio'] > 1 for r in rows),
                  improved_count=sum(0 < r['ratio'] < 1 for r in rows),
                  ratio_range=[min(r['ratio'] for r in rows), max(r['ratio'] for r in rows)],
                  branch_join=split, arclengths=lengths,
                  vertical_axis='q=1-delta_certified/delta_primal; primal0, dual1; standard upward axis',
                  numerical_marker_count=31, analytic_marker_count=1,
                  omitted_markers=0, shifted_markers=0, averaged_markers=0, clipped_ratios=0,
                  connecting_lines='Piecewise-linear visual guides, not certified interpolants',
                  connected_groups=[[r['cut_index'] for r in group] for group in groups],
                  new_certifications=0,
                  sha256_manifest={p.relative_to(ROOT).as_posix(): sha(p)
                                   for p in [SOURCE, GEO, PREVIOUS, *paths, Path(__file__)]},
                  output_sha256={p.name: sha(p) for p in [OUT/f'boundary_all_points.{ext}'
                                                       for ext in ['png', 'pdf', 'svg']]})
    (OUT/'plot_data.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    with (OUT/'all_markers.csv').open('w', newline='', encoding='utf-8') as f:
        columns = ['s_over_L', 'branch', 'percent', 'cut_index', 'analytic', 'ratio', 'ratio_exact',
                   'gap', 'primal_normal_distance', 'reason', 'source']
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sorted(rows+undefined, key=lambda r: r['s_over_L']))
    print(json.dumps({k: report[k] for k in ['total_source_markers', 'normalized_count',
                     'undefined_count', 'weak_count', 'improved_count', 'ratio_range',
                     'omitted_markers', 'connected_groups']}, indent=2))


if __name__ == '__main__':
    main()
