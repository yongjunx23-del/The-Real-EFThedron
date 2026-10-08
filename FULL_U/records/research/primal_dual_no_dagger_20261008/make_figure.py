"""Add a declared approximate left-tip continuation to the saved comparison."""
from pathlib import Path
from fractions import Fraction as F
import csv
import hashlib
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'output/primal_dual_connected_20261008/plot_data.json'
OUT = ROOT/'output/primal_dual_no_dagger_20261008'
LEFT = ROOT/'output/competitive_20261007/left/complete_left_certificate.json'
DUAL = ROOT/'research/benchmarks_manuscript_dual_upper.csv'
PRIMAL = ROOT/'research/benchmarks_manuscript_primal_upper.csv'
POLYGON = ROOT/'output/six_targets_20261007/comparison/certified_polygon.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = json.loads(SOURCE.read_text(encoding='utf-8'))
    weak = sorted([r for r in data['displayed'] if r['ratio'] > 1], key=lambda r: r['ratio'])
    keep_weak = {r['cut_index'] for r in weak[:3]}
    rows = [r for r in data['displayed'] if r['ratio'] <= 1 or r['cut_index'] in keep_weak]
    kept = {r['cut_index'] for r in rows}
    omitted = [dict(r, display_omission='Weaker than the three retained above-dual points.')
               for r in data['displayed'] if r['cut_index'] not in kept]
    omitted += [dict(r, display_omission='Undefined-position strip removed at user request.')
                for r in data['undefined']]
    assert len(rows) == 21 and len(omitted) == 11
    assert sum(r['ratio'] > 1 for r in rows) == 3
    assert len(kept | {r['cut_index'] for r in omitted}) == 32

    certificate = json.loads(LEFT.read_text(encoding='utf-8'))
    assert certificate['status'] == 'PASS_COMPLETE_CONTINUUM_LEFT_BOUND'
    assert certificate['objective_exact'] == ['-1', '0']
    C = F(certificate['support_upper_outward_6dp'])
    assert C >= F(certificate['support_upper_exact'])
    with DUAL.open(newline='') as f:
        reader = csv.reader(f); next(reader)
        dmin = min(F(r[0]) for r in reader)
    with PRIMAL.open(newline='') as f:
        reader = csv.reader(f); next(reader)
        pmin = min(F(r[0]) for r in reader)
    qmin_exact = (-C-pmin)/(dmin-pmin)
    qmin = float(qmin_exact)
    assert 0 < qmin < 1
    low = sorted([r for r in rows if r['branch'] == 'lower'], key=lambda r: r['s_over_L'])
    upper = sorted([r for r in rows if r['branch'] == 'upper'], key=lambda r: r['s_over_L'])
    # Unwrap the seam at s/L=1 and interpolate a local normalized ordinate.
    # This does NOT reconstruct an unknown primal boundary or add certificates.
    control_x = [r['s_over_L'] for r in upper[-2:]] + [1.] + [1+r['s_over_L'] for r in low[:2]]
    control_q = [r['ratio'] for r in upper[-2:]] + [qmin] + [r['ratio'] for r in low[:2]]
    continuation = PchipInterpolator(control_x, control_q, extrapolate=False)
    assert np.allclose(continuation(control_x), control_q, rtol=0, atol=2e-15)
    assert abs(float(continuation.derivative()(1.))) < 1e-14
    continuation_x = [np.linspace(upper[-1]['s_over_L'], 1., 201),
                      np.linspace(0., low[0]['s_over_L'], 201)]
    continuation_q = [continuation(continuation_x[0]), continuation(1+continuation_x[1])]
    assert all(np.all((y > 0) & (y < 1)) for y in continuation_q)
    endpoint = dict(s_over_L=[0, 1], ratio=qmin, ratio_exact=str(qmin_exact),
                    cut_index=21, source=str(LEFT.relative_to(ROOT)),
                    certified_lower_g0=str(-C), dual_min_g0=str(dmin), primal_min_g0=str(pmin),
                    definition='q_min=(g0_cert_min-g0_primal_min)/(g0_dual_min-g0_primal_min)',
                    caveat='Endpoint uses normalized extremal g0 support; interior markers use local normal distance.')

    # Join every retained marker in boundary order. These segments are visual
    # guides, including across omitted markers and the unobserved right-tip join.
    rows.sort(key=lambda r: r['s_over_L'])
    groups = [rows]
    assert all(a['s_over_L'] < b['s_over_L'] for a, b in zip(rows, rows[1:]))

    packet = json.loads(POLYGON.read_text(encoding='utf-8'))
    assert packet['bounded_all_directions']
    vertices_exact = [tuple(map(F, v)) for v in packet['vertices_exact']]
    cuts_exact = [tuple(F(c[k]) for k in ['a_exact', 'b_exact', 'C_exact']) for c in packet['cuts']]
    assert all(a*x+b*y <= c for x, y in vertices_exact for a, b, c in cuts_exact)
    for p, q in zip(vertices_exact, vertices_exact[1:]+vertices_exact[:1]):
        assert any(a*p[0]+b*p[1] == c == a*q[0]+b*q[1] for a, b, c in cuts_exact)
    vertices = np.array([[float(x), float(y)] for x, y in vertices_exact])
    closed_vertices = np.vstack([vertices, vertices[0]])
    reference_paths = {kind: {side: ROOT/f'research/benchmarks_manuscript_{kind}_{side}.csv'
                             for side in ['upper', 'lower']} for kind in ['dual', 'primal']}
    reference = {kind: {side: np.loadtxt(path, delimiter=',', skiprows=1)
                        for side, path in paths.items()} for kind, paths in reference_paths.items()}
    # The isolated published maximum-g0 reference was already part of the saved
    # comparison. The two straight connections to it are graphical approximations.
    tip = np.array([2.7272, .06])
    cap = np.vstack([reference['dual']['lower'][-1], tip, reference['dual']['upper'][-1]])
    dual_outline = np.vstack([reference['dual']['lower'], tip, reference['dual']['upper'][::-1]])
    primal_outline = np.vstack([reference['primal']['lower'], reference['primal']['upper'][::-1]])

    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['STIXGeneral'], 'mathtext.fontset': 'stix',
        'font.size': 9, 'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
        'legend.fontsize': 9, 'axes.linewidth': .65,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.size': 3.2, 'ytick.major.size': 3.2,
        'xtick.minor.size': 1.8, 'ytick.minor.size': 1.8,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(7.05, 5.65))
    upper_ax = fig.add_axes([.110, .555, .874, .380])
    ax = fig.add_axes([.110, .165, .874, .275])
    green = '#00785b'
    split = data['branch_join']

    # Match the lower panel: red is outside the published dual reference,
    # green is inside the primal reference, and their gap is left white.
    upper_ax.set_facecolor('#f4e5e3')
    upper_ax.fill(dual_outline[:, 0], dual_outline[:, 1], color='white', lw=0, zorder=0)
    upper_ax.fill(primal_outline[:, 0], primal_outline[:, 1], color='#e5efe8', lw=0, zorder=1)
    for side in ['upper', 'lower']:
        d = reference['dual'][side]
        p = reference['primal'][side]
        upper_ax.plot(d[:, 0], d[:, 1], color='#252525', lw=.95, zorder=3)
        upper_ax.plot(p[:, 0], p[:, 1], color='#3d6a9c', lw=1.1, ls=(0, (3, 2)), zorder=3)
    upper_ax.plot(cap[:, 0], cap[:, 1], color='#252525', lw=.9, ls=(0, (1.3, 1.7)), zorder=3)
    upper_ax.plot(closed_vertices[:, 0], closed_vertices[:, 1], color=green, lw=1.15, zorder=4)
    upper_ax.text(.1, .89, 'Dual Excluded', ha='center', va='center', fontsize=10, color='#743d38')
    upper_ax.text(-3.8, .54, 'Primal Allowed', ha='center', va='center', fontsize=9, color='#315b41')
    upper_ax.text(.975, .92, '(a)', transform=upper_ax.transAxes, ha='right', va='center', fontsize=9)
    upper_ax.set(xlim=(-8.4, 3.05), ylim=(-.035, 1.16), xlabel=r'$g_0$', ylabel=r'$g_2$')
    upper_ax.set_xticks([-8, -6, -4, -2, 0, 2])
    upper_ax.set_yticks(np.arange(0, 1.11, .2))
    upper_ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    upper_ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    upper_ax.tick_params(which='both', top=True, right=True, width=.6)

    def draw_points(axis, inset=False):
        for group in groups:
            axis.plot([r['s_over_L'] for r in group], [r['ratio'] for r in group],
                      color=green, lw=.70, ls=(0, (2.4, 2)), alpha=.72, zorder=3)
        axis.scatter([r['s_over_L'] for r in rows], [r['ratio'] for r in rows],
                     s=11 if inset else 13, color=green, marker='o',
                     edgecolors='white', linewidths=.3, zorder=4)

    ax.axhspan(1, 1.47, color='#f4e5e3', zorder=0)
    ax.axhspan(-.26, 0, color='#e5efe8', zorder=0)
    dual = ax.axhline(1, color='#252525', linewidth=.9, label='Dual (Ref.)', zorder=2)
    primal = ax.axhline(0, color='#3d6a9c', linewidth=1.1,
                       linestyle=(0, (3, 2)), label='Primal (Ref.)', zorder=2)
    draw_points(ax)
    for xx, yy in zip(continuation_x, continuation_q):
        ax.plot(xx, yy, color='#759888', lw=1.1, ls=(0, (1.3, 1.7)), zorder=3)
    ax.scatter([0, 1], [qmin, qmin], s=18, color=green, marker='o',
               edgecolors='white', linewidths=.3, zorder=5, clip_on=False)
    ax.axvline(split, color='#a0a0a0', linewidth=.6, linestyle=(0, (3, 3)), zorder=1)
    ax.text(.20, 1.28, 'Dual Excluded', ha='center', va='center', fontsize=10, color='#743d38')
    ax.text(.25, -.13, 'Primal Allowed', ha='center', va='center', fontsize=9, color='#315b41')
    for x, label, ha in [(0, r'$g_{0,\min}$', 'left'),
                         (split, r'$g_{0,\max}$', 'center'),
                         (1, r'$g_{0,\min}$', 'right')]:
        ax.text(x, 1.025, label, transform=ax.get_xaxis_transform(), ha=ha,
                va='bottom', fontsize=10)
    ax.text(.27, .085, 'Lower', ha='center', va='center', fontsize=8.5, color='#686868')
    ax.text(.77, .085, 'Upper', ha='center', va='center', fontsize=8.5, color='#686868')
    ax.text(.975, .90, '(b)', transform=ax.transAxes, ha='right', va='center', fontsize=9)
    ax.set(xlim=(0, 1), ylim=(-.26, 1.47), ylabel=r'Normalized bound, $q$',
           xlabel=r'Boundary arclength, $s/L$')
    ax.set_yticks([0, .25, .5, .75, 1, 1.25])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='both', top=True, right=True, width=.6)
    handle = Line2D([], [], color=green, lw=.70, ls=(0, (2.4, 2)), marker='o',
                    markersize=3.5, markeredgewidth=.3, markeredgecolor='white', label='This work, certified dual')
    approx = Line2D([], [], color='#759888', lw=1.1, ls=(0, (1.3, 1.7)),
                    label='Approx. continuation')
    fig.legend(handles=[handle, dual, primal], loc='lower center',
               bbox_to_anchor=(.547, .015), frameon=False, ncol=3,
               handlelength=1.7, columnspacing=1.4, borderpad=.15)

    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(OUT/f'boundary_comparison.{ext}', **({'dpi': 400} if ext == 'png' else {}))
    plt.close(fig)

    caption = r'''(a) Massive scalar bounds in four dimensions. The green boundary is the polygon obtained by intersecting all 60 stored completed inequalities (39 vertices); every displayed edge lies on a saved certified inequality. Black and blue dashed curves are the published dual and primal references. Green shading marks the primal interior; red shading marks the exterior of the dual reference. At the right tip, the two dotted black segments joining the observed dual endpoints to its isolated maximum-$g_0$ reference are approximate graphical continuations. (b) Normalized boundary comparison, with an approximate continuation at minimum $g_0$. At interior markers the ordinate is $q=1-\delta_C/\delta_P$, with signed outward normal displacements from the dual reference in axis-range-normalized coordinates. At the minimum-$g_0$ endpoint, shown at both ends of the counterclockwise loop, use $q_{\min}=(g_{0,C}^{\min}-g_{0,P}^{\min})/(g_{0,D}^{\min}-g_{0,P}^{\min})=0.72559745$, from the certified bound $g_0\geq-7.751502$. The endpoint uses an extremal support comparison, not the interior normal-distance definition. A shape-preserving cubic interpolation in $q$ joins this endpoint to the nearest interior markers; these segments are explicitly approximate and reconstruct neither the primal boundary nor a newly certified bound. Primal remains $q=0$ and dual $q=1$. Smaller values indicate tighter bounds. The earlier display filter retains 18 defined below-dual markers and three above-dual markers; the minimum-$g_0$ certificate is also represented at the endpoint. Every displayed marker is connected in boundary order. Dashed connections in panel (b) guide the eye without certifying interpolation, including across omitted stations and the unresolved right-tip join. Unlike the interpolated guides in panel (b), the green polygon edges in panel (a) follow the saved exact half-planes. Published references and the top-panel shading are graphical comparisons, not new certified reference data.
'''
    (OUT/'caption.txt').write_text(caption, encoding='utf-8')
    # The prior minimum-g0 position stays omitted; the same certificate now has
    # an explicitly different, valid endpoint comparison rather than a fake ratio.
    for r in omitted:
        if r['cut_index'] == 21:
            r['display_omission'] = 'Old normal-distance position undefined; same certificate now represented by extremal-support endpoint.'
    report = dict(status='COMBINED_COMPARISON_NO_DAGGER',
                  source=str(SOURCE.relative_to(ROOT)), displayed=rows, omitted=omitted,
                  displayed_count=len(rows), improved_count=18, weak_count=3,
                  retained_weak_values=[r['ratio'] for r in weak[:3]],
                  minimum_g0_endpoint=endpoint,
                  unique_displayed_certificates=len(kept | {21}),
                  visual_marker_count=len(rows)+2,
                  approximate_continuation=dict(method='PCHIP of q around the unwrapped left-tip seam; no extrapolation',
                      control_s_over_L=control_x, control_q=control_q,
                      pieces=[dict(s_over_L=x.tolist(), q=y.tolist())
                              for x, y in zip(continuation_x, continuation_q)],
                      scope='Approximate graphical interpolation only. No new certificates or feasible primal points.'),
                  filter='Keep all defined q<=1 and the three smallest defined q>1; remove undefined strip.',
                  connected_groups=[[r['cut_index'] for r in group] for group in groups],
                  branch_join=split, normalization=data['vertical_axis'],
                  figure_inches=[7.05, 5.65], legend_location='Below both panels, outside axes',
                  legend_labels=['This work, certified dual', 'Dual (Ref.)', 'Primal (Ref.)'],
                  region_labels=['Dual Excluded', 'Primal Allowed'],
                  bottom_primal_strip_q_range=[-.26, 0],
                  top_panel=dict(certified_polygon=str(POLYGON.relative_to(ROOT)),
                      exact_vertex_feasibility_check=True, exact_edge_active_cut_check=True,
                      cuts=len(cuts_exact), vertices=len(vertices_exact),
                      vertices_exact=packet['vertices_exact'],
                      published_references={kind: {side: pts.tolist() for side, pts in sides.items()}
                                            for kind, sides in reference.items()},
                      approximate_dual_right_cap=cap.tolist(),
                      isolated_right_tip_source='Existing reference in research/all_direction_certification_20261007/plot_competitive_certificates.py',
                      shading='Red outside published dual outline, green inside primal; approximate closure at dual right tip.',
                      polygon_sha256=sha(POLYGON),
                      reference_sha256={str(p.relative_to(ROOT)): sha(p)
                                        for paths in reference_paths.values() for p in paths.values()}),
                  all_displayed_dots_connected=True, connecting_lines_certified=False,
                  new_certifications=0, input_sha256=sha(SOURCE),
                  reference_sha256={str(p.relative_to(ROOT)): sha(p) for p in [LEFT, DUAL, PRIMAL]},
                  script_sha256=sha(Path(__file__)),
                  output_sha256={p.name: sha(p) for p in [OUT/f'boundary_comparison.{ext}'
                                                       for ext in ['png', 'pdf', 'svg']]})
    (OUT/'plot_data.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    with (OUT/'displayed_points.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['s_over_L', 'branch', 'cut_index', 'ratio', 'ratio_exact', 'source'],
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    np.savetxt(OUT/'certified_boundary.csv', closed_vertices, delimiter=',',
               header='g0,g2', comments='', fmt='%.17g')
    with (OUT/'approximate_continuation.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['s_over_L', 'q_approximate'])
        for x, y in zip(continuation_x, continuation_q):
            writer.writerows(zip(x, y))
    (OUT/'continuation_method.txt').write_text(
        'Approximate display near the left tip. Existing normalized ordinates remain unchanged. '
        'The independently certified minimum-g0 inequality is compared with the minimum-g0 values '
        'of the saved dual and primal references, giving q_min=' + str(qmin) + '. '
        'This endpoint normalization differs from interior local normal distances. '
        'Unwrap boundary arclength across s/L=1. Build a shape-preserving PCHIP through the last two '
        'upper markers, the endpoint, and the first two lower markers. Draw only the two missing '
        'left-tip intervals. The interpolation is approximate and adds no certificate markers, '
        'physical constraints, or primal feasibility claims. The endpoint is one certificate '
        'drawn twice because the boundary parameter is cyclic. All interpolation controls, '
        'samples and hashes are saved in plot_data.json. All other displayed dots are joined in order, '
        'including across omitted stations and the right-tip join; these lines are graphical '
        'guides and do not reconstruct the missing physical geometry or its arclength.\n',
        encoding='utf-8')
    print(json.dumps({k: report[k] for k in ['displayed_count', 'improved_count', 'weak_count',
                                            'retained_weak_values', 'minimum_g0_endpoint',
                                            'unique_displayed_certificates']}, indent=2))


if __name__ == '__main__':
    main()
