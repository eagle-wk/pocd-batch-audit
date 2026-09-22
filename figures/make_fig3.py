"""Fig 3 (tool head-to-head) for the methodological-audit manuscript.

Rebuilt 2026-09-21 at printed size and re-sourced. Two problems with the v1
version are fixed here:

  1. Panel B was a hand-built `warned = [0, 0, 0, 0] / warned[3] = 1` bar chart.
     Nothing in the benchmark JSON records whether ComBat/sva/RUV warn, so the
     panel was asserting a software property as if it were a measurement. It is
     replaced by the audit's own per-scenario verdict distribution, which is a
     real recorded field (`verdict_dist`). The "none of the three tools has any
     design-admissibility output" statement belongs in the legend and in the
     text, not in a chart.
  2. The legend said the figure covered scenario C only, while panel A plotted
     all five scenarios. Panel A keeps all five; the legend is corrected in the
     manuscript to match.

Panel A marks scenarios A and E as n/a rather than NaN: no correction is
applied there, so a retention value does not exist and plotting 0 would read as
"signal destroyed".

Run:
  python make_fig3.py
"""
import argparse, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))        # <repo>/figures
REPO = os.path.dirname(HERE)                            # <repo>
DERIVED = os.path.join(REPO, 'results', 'derived')
sys.path.insert(0, HERE)
import figstyle as fs                                 # noqa: E402
fs.use()

import matplotlib.patches as mpatches                 # noqa: E402
import matplotlib.pyplot as plt                       # noqa: E402
import numpy as np                                    # noqa: E402

sys.stdout.reconfigure(encoding='utf-8')

SCEN_LABEL = {
    'A_no_batch':       'A\nno batch',
    'B_batch_crossed':  'B\ngenuine batch,\ncrossed',
    'C_batch_confound': 'C\nbatch \u2261 group',
    'D_pseudo_batch':   'D\npseudo-batch',
    'E_global_shift':   'E\nglobal shift',
}
TOOL_COLOR = {'ComBat': '#d1495b', 'sva (supervised)': '#edae49', 'RUV': '#00798c'}
VERDICT_COLOR = {'proceed': fs.GREEN, 'correct': fs.AMBER, 'stop': fs.RED,
                 'normalize': fs.VIOLET}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', default=os.path.join(
        DERIVED, 'benchmark', 'benchmark_tool_compare.json'))
    ap.add_argument('--outdir', default=os.path.join(HERE, 'output'))
    args = ap.parse_args()

    data = json.load(open(args.json, encoding='utf-8'))
    rows = data['rows']
    scen = [r['scenario'] for r in rows]
    n = len(scen)
    x = np.arange(n)
    w = 0.26

    fig = plt.figure(figsize=(fs.FIG_W, 5.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.55,
                          left=0.085, right=0.985, top=0.90, bottom=0.235)

    # ---- Panel A : true-signal retention after correction -----------------
    axA = fig.add_subplot(gs[0])
    tools = ['ComBat', 'sva (supervised)', 'RUV']
    keys = ['combat_retention', 'sva_retention', 'ruv_retention']
    for k, (tool, key) in enumerate(zip(tools, keys)):
        off = (k - 1) * w
        vals, xs = [], []
        for i, r in enumerate(rows):
            v = r.get(key)
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            vals.append(v * 100)
            xs.append(i + off)
            if v * 100 > 100:       # ComBat inflates in B; flag it
                axA.text(i + off, v * 100 + 2, f'{v * 100:.0f}%', ha='center',
                         va='bottom', fontsize=7.2, color='#111', rotation=90)
            else:
                axA.text(i + off, v * 100 + 2, f'{v * 100:.1f}%', ha='center',
                         va='bottom', fontsize=7.2, color='#111')
        axA.bar(xs, vals, w, color=TOOL_COLOR[tool], edgecolor='#333',
                linewidth=0.7, label=tool, zorder=3)

    # scenarios with no correction applied
    for i, r in enumerate(rows):
        if all((r.get(k) is None or np.isnan(r.get(k, np.nan))) for k in keys):
            axA.axvspan(i - 0.45, i + 0.45, color='#f2f2f2', zorder=0)
            axA.text(i, 62, 'n/a\nno correction\napplied', ha='center',
                     va='center', fontsize=7.2, color='#888', linespacing=1.35)

    axA.axhline(100, color='#999', lw=0.9, ls=':', zorder=1)
    axA.axhline(30, color='#666', lw=1.0, ls='--', zorder=1)
    axA.text(n - 0.52, 32, '30% retention reference', ha='right', va='bottom',
             fontsize=7.2, color='#666')
    axA.set_xticks(x)
    axA.set_xticklabels([SCEN_LABEL[s] for s in scen], fontsize=7.5)
    axA.set_xlim(-0.55, n - 0.45)
    axA.set_ylim(0, 132)
    axA.set_yticks([0, 25, 50, 75, 100])
    axA.set_ylabel('true signal retained\nafter correction  (%)')
    axA.set_title('A  Real-signal survival under batch correction', fontsize=8.8,
                  pad=6)
    for s in ('top', 'right'):
        axA.spines[s].set_visible(False)
    axA.grid(axis='y', alpha=0.22, lw=0.6, zorder=0)
    axA.legend(fontsize=7.2, loc='upper center', ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.02))

    # ---- Panel B : the audit's verdict for the same five scenarios --------
    axB = fig.add_subplot(gs[1])
    ORDER = ['proceed', 'correct', 'stop', 'normalize']
    bottoms = np.zeros(n)
    for v in ORDER:
        heights = np.array([r['verdict_dist'].get(v, 0) for r in rows], float)
        axB.bar(x, heights, 0.55, bottom=bottoms, color=VERDICT_COLOR[v],
                edgecolor='#333', linewidth=0.7, label=v, zorder=3)
        for i, h in enumerate(heights):
            if h >= 12:
                axB.text(i, bottoms[i] + h / 2, f'{v}\n{int(h)}/100', ha='center',
                         va='center', fontsize=7.2,
                         color='white' if v != 'correct' else '#3a2c00',
                         weight='bold', linespacing=1.3)
            elif h > 0:
                axB.text(i + 0.34, bottoms[i] + h / 2, f'{int(h)}', ha='left',
                         va='center', fontsize=7.2, color='#333')
        bottoms += heights
    axB.set_xticks(x)
    axB.set_xticklabels([SCEN_LABEL[s] for s in scen], fontsize=7.5)
    axB.set_xlim(-0.55, n - 0.45)
    axB.set_ylim(0, 118)
    axB.set_yticks([0, 25, 50, 75, 100])
    axB.set_ylabel('replicates  (%)')
    axB.set_title('B  Verdict returned by the audit for the same five scenarios',
                  fontsize=8.8, pad=6)
    for s in ('top', 'right'):
        axB.spines[s].set_visible(False)
    axB.grid(axis='y', alpha=0.22, lw=0.6, zorder=0)

    fig.suptitle('Fig. 3  Correction tools erase real signal under confounding '
                 'and say nothing', fontsize=9.8, weight='bold', y=0.985)
    fig.text(0.5, 0.012, fs.wrap(
        f'{data["reps"]} replicates per scenario, {data["p_genes"]:,} genes, '
        f'n = {data["n_per_group"]} per group, {data["perm"]} permutations. '
        'sva is run in supervised mode (group labels supplied). In B, ComBat '
        'returns >100% because it inflates rather than preserves the contrast. '
        'None of the three tools has any design-admissibility output, so none of '
        'them can report scenario C as unfit; only the audit refuses it.',
        148), ha='center', va='bottom', fontsize=7.2, color='#555',
        linespacing=1.55)

    problems = {'Fig3_tool_compare': fs.save(fig, args.outdir, 'Fig3_tool_compare')}
    fs.report(problems, 1)


if __name__ == '__main__':
    main()
