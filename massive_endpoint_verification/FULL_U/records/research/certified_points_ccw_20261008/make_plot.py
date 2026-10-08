"""Presentation-only CCW remapping and declared thinning of certified markers."""
from pathlib import Path
from fractions import Fraction as F
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'output/certified_points_ccw_20261008'
OLD = ROOT/'output/certified_points_plot_20261007/plot_data.json'
GEOM = ROOT/'output/six_targets_20261007/comparison/geometry_continuous_audit.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    old = json.loads(OLD.read_text(encoding='utf-8'))
    geo = json.loads(GEOM.read_text(encoding='utf-8'))
    lengths = {s: geo['branches'][s]['observed_normalized_length'] for s in ['lower', 'upper']}
    lo, up = [lengths[s]['midpoint_float'] for s in ['lower', 'upper']]
    total = lo+up
    split = 100*lo/total
    split_lo = F(100)*F(lengths['lower']['lower_exact'])/(F(lengths['lower']['lower_exact'])+F(lengths['upper']['upper_exact']))
    split_hi = F(100)*F(lengths['lower']['upper_exact'])/(F(lengths['lower']['upper_exact'])+F(lengths['upper']['lower_exact']))

    def ccw(side, percent):
        return percent*lo/total if side == 'lower' else 100*(lo+(1-percent/100)*up)/total

    assert ccw('lower', 0) == 0 and abs(ccw('upper', 0)-100) < 1e-12
    assert abs(ccw('lower', 100)-split) < 1e-12 and abs(ccw('upper', 100)-split) < 1e-12

    # User-requested display filter. Proof data and old figures are untouched.
    selected, removed, crowded = [], [], []
    for marker in old['markers']:
        row = dict(marker, ccw_percent=ccw(marker['branch'], marker['percent']))
        if marker['gap'] > 0:
            removed.append(dict(row, display_omission='Weaker than digitized dual at this displayed point (positive gap).'))
        elif marker['branch'] == 'lower' and marker['percent'] > 80:
            crowded.append(row)
        else:
            selected.append(row)
    kept_crowded = []
    for row in sorted(crowded, key=lambda r: r['gap']):
        if all(abs(row['percent']-r['percent']) >= 5 for r in kept_crowded):
            kept_crowded.append(row)
        else:
            removed.append(dict(row, display_omission='Crowding: within 5 original lower-branch percentage points of an already retained stronger marker.'))
    selected.extend(kept_crowded)
    selected.sort(key=lambda r: r['ccw_percent'])
    assert len(selected)+len(removed) == len(old['markers'])
    assert len(selected) == 19 and len(kept_crowded) == 2
    assert all(r['gap'] <= 0 and not r['analytic'] for r in selected)
    assert len({r['cut_index'] for r in selected}) == len(selected)
    for a in kept_crowded:
        for b in kept_crowded:
            assert a is b or abs(a['percent']-b['percent']) >= 5

    traces = []
    for side in ['lower', 'upper']:
        pieces = old['primal_normal_gap_traces'][side]
        for trace in pieces if side == 'lower' else reversed(pieces):
            pts = trace['points'] if side == 'lower' else list(reversed(trace['points']))
            rows = [dict(p, ccw_percent=ccw(side, p['percent'])) for p in pts]
            assert all(a['ccw_percent'] <= b['ccw_percent']+1e-12 for a, b in zip(rows, rows[1:]))
            traces.append(dict(branch=side, dual_segment=trace['dual_segment'], points=rows))

    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(12.8, 6.4))
    first = True
    for trace in traces:
        p = trace['points']
        ax.plot([r['ccw_percent'] for r in p], [np.nan if r['gap'] is None else r['gap'] for r in p],
                color='#4777b6', ls=':', lw=2.2,
                label='Published primal (digitized)' if first else None, zorder=2)
        first = False
    ax.axhline(0, color='#303030', lw=1.2, label='Published dual = 0 (digitized)', zorder=1)
    ax.scatter([r['ccw_percent'] for r in selected], [r['gap'] for r in selected],
               s=48, color='#14845d', edgecolors='white', linewidths=.6,
               label='Certified bounds (selected markers)', zorder=4)
    ax.axvline(split, color='#777777', lw=1.1, ls=(0, (4, 4)), zorder=1)
    ax.text(0, 1.022, r'Minimum $g_0$', transform=ax.get_xaxis_transform(), ha='left', va='bottom', fontsize=11)
    ax.text(split, 1.022, r'Maximum $g_0$ — omitted tip arc', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=11)
    ax.text(100, 1.022, r'Minimum $g_0$ (return)', transform=ax.get_xaxis_transform(), ha='right', va='bottom', fontsize=11)
    ax.text(split/2, .968, 'Lower boundary (left → right)', transform=ax.get_xaxis_transform(), ha='center', va='top', color='#686868', fontsize=10)
    ax.text((split+100)/2, .968, 'Upper boundary (right → left)', transform=ax.get_xaxis_transform(), ha='center', va='top', color='#686868', fontsize=10)
    ax.set(xlim=(0, 100), ylim=(-.087, .010), xlabel='Counter-clockwise normalized observed boundary length (%)',
           ylabel='Signed outward normal gap from dual')
    ax.set_xticks(np.arange(0, 101, 10))
    ax.grid(alpha=.2)
    ax.legend(loc='lower left', fontsize=10, framealpha=.95)
    fig.suptitle('Certified bounds and published primal around the boundary', fontsize=15, y=.98)
    fig.text(.5, .090, 'Certified markers below zero improve on the digitized dual. Weaker markers are hidden; crowded lower-arc markers are thinned.', ha='center', fontsize=9)
    fig.text(.5, .059, 'The omitted right-tip arc has unknown length; its maximum-g₀ side is marked at the 53.15% observed-arc join.', ha='center', fontsize=9)
    fig.text(.5, .028, 'Primal blanks near 0% and 100% mean those normal lines miss the smaller primal region. No missing values are filled with zero.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .13, 1, .93))
    fig.savefig(OUT/'certified_points_ccw_with_primal.png', dpi=180)
    fig.savefig(OUT/'certified_points_ccw_with_primal.pdf')
    plt.close(fig)

    manifest = {p.relative_to(ROOT).as_posix(): sha(p) for p in [OLD, GEOM, Path(__file__)]}
    data = dict(status='PRESENTATION_ONLY_CCW_REMAP_WITH_EXPLICIT_DISPLAY_FILTER',
                selected_markers=selected, omitted_markers=removed, primal_traces_ccw=traces,
                lengths=lengths, branch_join_percent=split,
                branch_join_percent_enclosure=[str(split_lo), str(split_hi)],
                branch_order='Start at common observed minimum g0; lower left-to-right, omitted right cap, upper right-to-left; return to minimum g0.',
                maximum_g0_label='The true maximum is within an undigitized cap. The 53.15% junction is its side of the observed-arc parameter, not an exact full-perimeter location.',
                filter='Remove every positive-gap marker, then below-dual lower markers beyond original80% are greedily thinned at5-percentage-point separation, strongest gap first.',
                styling='All retained certified markers use the same green circular symbol. No latest distinction, stars, or connecting certificate line.',
                new_certifications=0, sha256_manifest=manifest)
    (OUT/'plot_data.json').write_text(json.dumps(data, indent=2), encoding='utf-8')
    with (OUT/'displayed_certified_points.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['ccw_percent', 'branch', 'percent', 'gap', 'cut_index', 'source'], extrasaction='ignore')
        writer.writeheader(); writer.writerows(selected)
    (OUT/'README.md').write_text(f'''# Counter-clockwise certified-point plot

This is a display revision of the existing proof-backed marker data. All retained certified points use identical circles. No special latest/star style remains, and no line connects certificate markers.

## Coordinate

Start at the common observed dual minimum g0, traverse the lower boundary left to right, then the upper boundary right to left, returning to the minimum. This is counter-clockwise. As before, each coordinate is divided by its observed dual range before measuring length. Observed lower and upper lengths are {lo:.16g} and {up:.16g}; their join is {split:.12f}% of the combined observed length. The two branches do not receive equal artificial lengths.

The observed dual minimum is the common left endpoint, approximately (-8.014558, 1.081685). The observed dual branch endpoints on the right both have g0 approximately 2.656293 but distinct g2 values. The true maximum-g0 side lies in the missing arc joining them. That missing arc has unknown length. The figure labels its side at the branch join; it does not claim that the true maximum occurs at exactly {split:.6f}% of the full, partly unobserved perimeter. Neither primal nor certificate traces are joined across the omitted cap.

## Why the primal does not extend to zero

The y-coordinate is a signed distance along the local outward normal line of the digitized dual boundary. The published primal is a smaller region, with its observed left tip approximately (-7.055908, 0.929694). Near the sharp dual left corner, the normal lines miss the entire observed primal boundary. The plotted distance is therefore undefined, even though both primal branches meet at a known left endpoint. These are geometric missing intersections, not a missing primal-left arc. They appear near both 0% and 100% in the counter-clockwise plot. They are left blank; zero would incorrectly assert equality with the dual. The primal also need not have zero vertical gap where it exists: it and the dual are different boundaries.

## Display filter

There are {len(selected)} displayed certificate markers. Every marker with positive gap (weaker than the digitized dual at that representative location) is hidden as requested. Among the remaining lower-branch markers beyond the original 80% position, retain the strongest first, then require at least five original lower-branch percentage points of separation. This keeps two markers in that crowded portion. The filter applies before converting to the combined counter-clockwise coordinate. All positions and ordinate values are inherited without alteration.

All {len(removed)} display omissions and their reasons are retained in plot_data.json. Hiding weak or crowded markers is only a presentation choice; it does not remove their certified inequalities, change the saved polygon, or establish dominance on unshown intervals. The published references remain digitized. No new bounds or certifications are claimed.

## Files

- plot_data.json: all selected and omitted markers, source hashes, remapped primal traces, length split and its inherited interval enclosure.
- displayed_certified_points.csv: the displayed numerical data.
- certified_points_ccw_with_primal.png and .pdf: figure exports.

Reproduce: `python research/certified_points_ccw_20261008/make_plot.py`.
''', encoding='utf-8')
    print(json.dumps(dict(selected_markers=len(selected), omitted_markers=len(removed), join_percent=split,
        retained_lower_after80=[r['percent'] for r in kept_crowded], figure=str(OUT/'certified_points_ccw_with_primal.png')), indent=2))


if __name__ == '__main__':
    main()
