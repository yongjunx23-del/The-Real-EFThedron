"""Point-only display of actual active certificates, with primal normal gaps.

Reads the frozen exact envelope audit. Does not optimize or certify new bounds.
"""
from pathlib import Path
from fractions import Fraction as F
import csv, hashlib, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/certified_points_plot_20261007'
BASE = ROOT / 'output/six_targets_20261007'
manifest = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    manifest[path.relative_to(ROOT).as_posix()] = sha(path)
    return json.loads(path.read_text(encoding='utf-8'))


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def cross(a, b):
    return a[0]*b[1] - a[1]*b[0]


def primal_intersection(p, ray, points):
    """Outermost intersection with the observed SAME branch; no connectors."""
    hits = []
    for j, (a, b) in enumerate(zip(points, points[1:])):
        edge = b-a
        det = cross(ray, edge)
        if abs(det) < 1e-14:
            continue
        offset = a-p
        t = cross(offset, edge)/det
        u = cross(offset, ray)/det
        if -1e-12 <= u <= 1+1e-12:
            assert np.linalg.norm(p+t*ray-(a+u*edge)) < 1e-10
            hits.append((float(t), j, float(u)))
    return max(hits) if hits else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    packet = read(BASE/'comparison/certified_polygon.json')
    audit = read(BASE/'comparison/geometry_continuous_audit.json')
    assert audit['status'] == 'PASS_CONTINUOUS_DIGITIZED_GEOMETRY_AUDIT_NOT_NEW_PHYSICAL_CERTIFICATE'
    assert audit['sha256_manifest'][audit['input_polygon']] == sha(BASE/'comparison/certified_polygon.json')
    selection = read(BASE/'CURRENT_CERTIFICATE_SELECTION.json')
    target_map = read(BASE/'target_map.json')
    targets = {t['tag']: t for t in target_map['targets']}
    new_sources = {str(Path(x['path'])): targets[x['target']] for x in selection['certificates']}
    cuts = [tuple(F(c[k]) for k in ['a_exact', 'b_exact', 'C_exact']) for c in packet['cuts']]
    scale = np.array([float(F(x)) for x in audit['normalization']['axis_ranges_exact']])
    markers, primal_series, omitted = [], {}, []
    active_indices = set()

    for side, branch in audit['branches'].items():
        segments = branch['segments']
        # Exact audited active intervals; arclength is used only for placement.
        by_cut = {}
        for seg in segments:
            start = seg['branch_arc_interval']['start']['midpoint_float']
            stop = seg['branch_arc_interval']['end']['midpoint_float']
            for piece in seg['upper_envelope_pieces']:
                lo, hi = map(F, piece['t_interval_exact'])
                if lo == hi:
                    continue
                k = piece['active_cut']
                by_cut.setdefault(k, []).append(dict(
                    start=start+float(lo)*(stop-start), end=start+float(hi)*(stop-start),
                    segment=seg, piece=piece, lo=lo, hi=hi))

        for k, pieces in by_cut.items():
            active_indices.add(k)
            pieces.sort(key=lambda p: p['start'])
            groups = []
            for piece in pieces:
                if groups and abs(groups[-1][-1]['end']-piece['start']) < 1e-12:
                    groups[-1].append(piece)
                else:
                    groups.append([piece])
            longest = max(groups, key=lambda g: g[-1]['end']-g[0]['start'])
            position = (longest[0]['start']+longest[-1]['end'])/2
            source = packet['cuts'][k]['source']
            target = new_sources.get(str(Path(source)))
            placement = 'Midpoint of longest contiguous active interval'
            if target is not None:
                assert target['branch'] == side
                position = target['target_percent']/100
                placement = 'Requested target, verified inside this certificate active interval'
            containing = [p for p in pieces if p['start'] < position < p['end']]
            if not containing:
                # A midpoint may coincide with a reference corner. Choose a
                # definite interior point, never invent an interpolated normal.
                assert target is None
                chosen = min(longest, key=lambda p: abs((p['start']+p['end'])/2-position))
                position = (chosen['start']+chosen['end'])/2
                placement += '; adjusted to interior of a reference segment'
            else:
                chosen = containing[0]
            seg = chosen['segment']
            start = seg['branch_arc_interval']['start']['midpoint_float']
            stop = seg['branch_arc_interval']['end']['midpoint_float']
            t = F((position-start)/(stop-start))
            assert chosen['lo'] < t < chosen['hi']
            p0, p1 = [tuple(map(F, seg[q])) for q in ['start_exact', 'end_exact']]
            p = tuple(x+t*(y-x) for x, y in zip(p0, p1))
            ray = tuple(map(F, seg['physical_ray_unnormalized_exact']))
            a, b, C = cuts[k]
            den = dot((a, b), ray)
            assert den > 0
            u = (C-dot((a, b), p))/den
            assert u == F(chosen['piece']['intercept_exact'])+t*F(chosen['piece']['slope_exact'])
            hit = tuple(x+u*y for x, y in zip(p, ray))
            # The plotted ordinate is computed from this cut alone. This
            # exact check shows why it also lies on the existing envelope.
            assert all(aa*hit[0]+bb*hit[1] <= cc for aa, bb, cc in cuts)
            gap = float(u)*math.sqrt(float(F(seg['metric_length_squared_exact'])))
            analytic = not source.endswith('.json')
            if not analytic:
                file = ROOT/source
                assert file.is_file()
                manifest[file.relative_to(ROOT).as_posix()] = sha(file)
            markers.append(dict(branch=side, percent=100*position, gap=gap, cut_index=k,
                source=source, analytic=analytic, recent_target=None if target is None else target['tag'],
                placement=placement, active_interval_percent=[100*longest[0]['start'], 100*longest[-1]['end']],
                objective_exact=[str(a), str(b)], support_upper_exact=str(C),
                segment=seg['segment'], segment_t_exact=str(t), reference_point_exact=list(map(str, p)),
                physical_ray_unnormalized_exact=list(map(str, ray)), normal_parameter_exact=str(u),
                intersection_exact=list(map(str, hit)), own_cut_equals_envelope_exactly=True,
                all_cuts_satisfied_exactly=True))

        path = ROOT/f'research/benchmarks_manuscript_primal_{side}.csv'
        primal = np.loadtxt(path, delimiter=',', skiprows=1)
        manifest[path.relative_to(ROOT).as_posix()] = sha(path)
        traces = []
        for seg in segments:
            p0, p1 = [np.array([float(F(x)) for x in seg[q]]) for q in ['start_exact', 'end_exact']]
            ray = np.array([float(F(x)) for x in seg['physical_ray_unnormalized_exact']])
            ray /= math.sqrt(float(F(seg['metric_length_squared_exact'])))
            assert abs(np.linalg.norm(ray/scale)-1) < 1e-12
            start = seg['branch_arc_interval']['start']['midpoint_float']
            stop = seg['branch_arc_interval']['end']['midpoint_float']
            rows = []
            for t in np.linspace(0, 1, 41):
                p = p0+t*(p1-p0)
                hit = primal_intersection(p, ray, primal)
                rows.append(dict(percent=100*(start+t*(stop-start)), gap=None if hit is None else hit[0],
                    primal_segment=None if hit is None else hit[1], primal_segment_fraction=None if hit is None else hit[2]))
            traces.append(dict(dual_segment=seg['segment'], points=rows))
        primal_series[side] = traces

    for k, c in enumerate(packet['cuts']):
        if k not in active_indices:
            omitted.append(dict(cut_index=k, source=c['source'], reason='No positive-length active normal-gap interval on either observed branch; includes superseded cuts and unobserved corner/cap cuts.'))
    assert len(markers) == 32 and sum(not r['analytic'] for r in markers) == 31
    assert len([r for r in markers if r['recent_target']]) == 6

    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2))
    for ax, side in zip(axes, ['upper', 'lower']):
        first = True
        for trace in primal_series[side]:
            rr = trace['points']
            ax.plot([r['percent'] for r in rr], [np.nan if r['gap'] is None else r['gap'] for r in rr],
                color='#4777b6', lw=2, ls=':', label='Published primal (digitized)' if first else None, zorder=2)
            first = False
        ax.axhline(0, color='#303030', lw=1.2, label='Published dual = 0 (digitized)', zorder=1)
        selected = [r for r in markers if r['branch'] == side]
        for kind, color, shape, size, label in [
            ('earlier', '#14845d', 'o', 42, 'Completed directional certificates'),
            ('recent', '#d17728', '*', 155, 'Six latest certified targets'),
            ('analytic', '#777777', 'D', 42, r'Analytic bound $g_2\geq0$')]:
            rr = [r for r in selected if ('analytic' if r['analytic'] else 'recent' if r['recent_target'] else 'earlier') == kind]
            if rr:
                ax.scatter([r['percent'] for r in rr], [r['gap'] for r in rr], color=color,
                    marker=shape, s=size, edgecolors='white', linewidths=.65, label=label, zorder=4)
        for r in selected:
            if r['recent_target']:
                ax.annotate(f"{r['percent']:.0f}%", (r['percent'], r['gap']), xytext=(0, 10),
                    textcoords='offset points', ha='center', fontsize=9, color='#9a5217')
        count = sum(not r['analytic'] for r in selected)
        ax.set(title=f'{side.capitalize()} boundary — {count} active certificates', xlim=(0, 101),
            xlabel='Normalized branch length, left to right (%)', ylabel='Signed outward normal gap from dual')
        ax.grid(alpha=.2)
        ax.margins(y=.18)
        ax.legend(fontsize=8.5, loc='lower right')
    fig.suptitle('Certified bounds as individual points, with the published primal boundary', fontsize=15, y=.98)
    fig.text(.5, .082, 'Certified markers below zero improve on the digitized dual; above zero are weaker.', ha='center', fontsize=10)
    fig.text(.5, .049, 'One marker per active completed inequality; representative positions, with the six requested targets retained.', ha='center', fontsize=9)
    fig.text(.5, .019, 'No lines connect certified markers. Primal gaps use the same local normals; missing intersections and unobserved caps are omitted.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .115, 1, .94))
    fig.savefig(OUT/'certified_points_with_primal.png', dpi=180)
    fig.savefig(OUT/'certified_points_with_primal.pdf')
    plt.close(fig)
    manifest[Path(__file__).relative_to(ROOT).as_posix()] = sha(Path(__file__))
    data = dict(status='POINT_ONLY_GRAPHICAL_DISPLAY_OF_EXISTING_CERTIFICATES', markers=markers,
        primal_normal_gap_traces=primal_series, omitted_cuts=omitted, sha256_manifest=manifest,
        counts=dict(completed_directional_markers=31, separate_analytic_markers=1, new_certifications=0),
        limitation='These are representative displays of completed halfplanes, not physical feasible points or newly optimized stations. Published curves are digitized, with no rigorous extraction-error bars. Primal branches are not closed or extrapolated.')
    (OUT/'plot_data.json').write_text(json.dumps(data, indent=2), encoding='utf-8')
    with (OUT/'certified_markers.csv').open('w', newline='', encoding='utf-8') as f:
        cols = ['branch', 'percent', 'gap', 'cut_index', 'source', 'analytic', 'recent_target', 'placement']
        writer = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        writer.writeheader(); writer.writerows(markers)
    (OUT/'README.md').write_text('''# Certified points and published primal comparison

The two panels retain the previous axis-range-normalized branch-length coordinate and signed outward normal-line gap from the digitized published dual curve. Zero is the dual reference. Negative certificate values lie locally inward and improve the dual outer bound; positive certificate values are weaker locally. The primal curve is a separate published reference, not an additional certified outer bound.

Each circle or star represents one completed directional inequality that actually controls a positive-length portion of the observed boundary comparison. There are 7 upper and 24 lower certificate markers. A separate gray diamond displays the analytic inequality g2 >= 0. The six latest targets remain at upper 10%, 55%, 95% and lower 15%, 45%, 60%. Other markers use the midpoint of that inequality's longest active interval, adjusted to an interior segment if necessary. These positions represent a completed halfplane; they are not additional optimization locations or physical feasible points.

At each marker, its own inequality a.g <= C gives u=(C-a.p)/(a.r), where p is an exact rational reference point and r is that segment's unnormalized metric normal in physical coordinates. Its plotted gap is u times the normalized segment length. Exact rational checks establish that this inequality is active and that its intersection satisfies every saved cut. No dense samples or polygon vertices are promoted to extra certificate markers. Superseded cuts and cuts with no positive-length active interval on the observed branches are listed in plot_data.json and omitted from the figure.

The blue dotted primal reference intersects each local dual normal line with the observed primal polyline of the same branch. The outermost observed intersection is used. Each reference segment is drawn separately, preserving one-sided normals at corners. Missing intersections are gaps; no endpoints are extrapolated or missing caps connected. The published curves remain digitized graphical inputs, not exact author optimizer data.

The exact source halfplanes, representative rational coordinates, individual gap calculations, source hashes, omitted-cut inventory and primal traces are in plot_data.json. Marker values are also in certified_markers.csv. No new optimization or physical certification was performed. The previously sealed files and archive are unchanged.

Reproduce from the workspace root:

```text
python research/certified_points_plot_20261007/make_plot.py
```
''', encoding='utf-8')
    print(json.dumps(dict(counts=data['counts'], per_branch={s:sum(not r['analytic'] and r['branch']==s for r in markers) for s in ['upper','lower']}, figure=str(OUT/'certified_points_with_primal.png')), indent=2))


if __name__ == '__main__':
    main()
