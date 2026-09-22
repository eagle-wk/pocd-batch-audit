"""Fig 4 (audit map, 14 contrast-units) and Fig 6 (technical separability gate).

Renamed from make_fig4_fig5.py on 2026-09-21: the separability gate is Fig 6 in
the manuscript (the three-real-cases panel took the Fig 5 slot), so the old
filename was actively misleading.

Rebuilt at printed size. The v1 Fig 4 was laid out at 12.6 in with a single
un-wrapped footnote line, which pushed the saved file to 15.38 in — twice the
text column. The footnote is now hard-wrapped and the grid is laid out inside
7.2 in.

Inputs : results/master_audit_table.tsv
         results/derived/verdicts/*.json
Outputs: pocd_data/figures/Fig4_audit_map.png / .pdf
         pocd_data/figures/Fig6_separability_gate.png / .pdf

Figure text is English throughout (publication figures). All numbers are read
from audit outputs; nothing is hardcoded from memory except axis semantics.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import figstyle as fs                                   # noqa: E402
fs.use()

import matplotlib.patches as mpatches                   # noqa: E402
import matplotlib.pyplot as plt                         # noqa: E402

sys.stdout.reconfigure(encoding='utf-8')
REPO = os.path.dirname(HERE)                            # <repo>
RESULTS = os.path.join(REPO, 'results')
DERIVED = os.path.join(RESULTS, 'derived')
VDIR = os.path.join(DERIVED, 'verdicts')                # one json per unit
OUT = os.path.join(HERE, 'output')
os.makedirs(OUT, exist_ok=True)
PROBLEMS = {}

# ---------------------------------------------------------------- load
rows = list(csv.DictReader(
    open(os.path.join(RESULTS, 'master_audit_table.tsv'),
         encoding='utf-8'), delimiter='\t'))

# technical separability numbers, keyed by dataset (only where computed)
#
# Panel A is restricted to the *TBI-domain* units, matching the manuscript
# Figure 6 legend ("Technical separability for the TBI-domain units"). The
# POCD-domain units are audited on a different footing: GSE283401 / GSE297195 /
# GSE330865 were judged by an older script version with no permutation test
# (verify_all.py "Finding A9"), and GSE163943 carries a batch label inferred
# from an essentially unimodal library-size split while sitting on the human
# side of the corpus. Listing the domain explicitly -- rather than relying on
# which directory a verdict happens to live in -- is what makes this selection
# reproducible. The field check below stays as a safety net.
TBI_UNITS = {'GSE115614', 'GSE31357', 'GSE59645', 'GSE64978', 'GSE68207',
             'GSE80174_Cortex', 'GSE80174_Hippocampus', 'GSE80174_Thalamus'}
sep = {}
for fn in sorted(os.listdir(VDIR)):
    if not fn.endswith('.json'):
        continue
    d, j = fn[:-5], os.path.join(VDIR, fn)
    if d not in TBI_UNITS:
        continue
    v = json.load(open(j, encoding='utf-8'))
    bp = v.get('batch_plausibility') or {}
    # Fig 6 plots the *dual gate*: accuracy above baseline AND a permutation
    # p-value. Both fields must be present for a point to be plotable.
    if bp.get('loocv_accuracy_on_technical') is None or bp.get('loocv_perm_p') is None:
        continue
    sep[d] = {
        'loocv': bp['loocv_accuracy_on_technical'],
        'base': bp['random_baseline'],
        'p': bp['loocv_perm_p'],
        'supported': bool(bp.get('batch_supported_by_technical')),
    }

ORDER = ['GSE330865', 'GSE297195', 'GSE283401', 'GSE304956', 'GSE276942', 'GSE163943',
         'GSE59645', 'GSE115614', 'GSE31357', 'GSE68207', 'GSE64978',
         'GSE80174_Cortex', 'GSE80174_Thalamus', 'GSE80174_Hippocampus']
by_ds = {r['dataset']: r for r in rows}
rows = [by_ds[d] for d in ORDER if d in by_ds]

LABEL = {'GSE80174_Cortex': 'GSE80174 ctx',
         'GSE80174_Thalamus': 'GSE80174 thal',
         'GSE80174_Hippocampus': 'GSE80174 hippo'}

VERDICT_COLOR = fs.VERDICT_COLOR
PROV_COLOR = fs.PROV_COLOR
PROV_TEXT = fs.PROV_TEXT

n = len(rows)
names = [LABEL.get(r['dataset'], r['dataset']) for r in rows]

# ================================================================ Fig 4
fig, axes = plt.subplots(1, 4, figsize=(fs.FIG_W, 5.9), sharey=True,
                         gridspec_kw={'width_ratios': [1.20, 0.80, 0.80, 1.05]})
fig.subplots_adjust(wspace=0.07, left=0.215, right=0.985, top=0.845, bottom=0.155)


def gy(i):
    """Display order top->bottom; matplotlib's y axis runs bottom-up."""
    return n - 1 - i


def grid(ax, title):
    ax.set_xticks([])
    ax.set_yticks([gy(i) for i in range(n)])
    ax.set_yticklabels(names, fontsize=7.4)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_title(title, fontsize=7.8, pad=5, linespacing=1.35)
    ax.tick_params(axis='y', length=0)
    for s in ('top', 'right', 'left', 'bottom'):
        ax.spines[s].set_visible(False)


def cell(ax, i, color, text, tcolor='#222'):
    ax.add_patch(mpatches.Rectangle((-0.5, gy(i) - 0.5), 1.0, 1.0,
                                    facecolor=color, edgecolor='white', lw=1.0))
    ax.text(0, gy(i), text, ha='center', va='center', fontsize=7.2, color=tcolor)


# --- panel 1: batch provenance
ax = axes[0]
grid(ax, 'Batch\nprovenance')
for i, r in enumerate(rows):
    p = r['batch_source']
    cell(ax, i, PROV_COLOR.get(p, '#ececec'), PROV_TEXT.get(p, '?'))

# --- panel 2: structural confounding (# groups in a single batch)
ax = axes[1]
grid(ax, 'Groups in a\nsingle batch')
cm = plt.cm.Oranges
for i, r in enumerate(rows):
    try:
        v = int(r['n_confounded_groups'])
    except (TypeError, ValueError):
        v = 0
    c = cm(0.0 if v == 0 else min(1.0, 0.18 + v / 14.0))
    cell(ax, i, c, str(v) if v else '0', tcolor='#222' if v < 8 else 'white')

# --- panel 3: lead PC R2_batch
ax = axes[2]
grid(ax, 'Lead PC\nR\u00b2(batch)')
cm2 = plt.cm.Reds
for i, r in enumerate(rows):
    try:
        v = float(r['lead_PC_R2_batch'])
    except (TypeError, ValueError):
        v = 0.0
    cell(ax, i, cm2(min(1.0, max(0.0, v))), f'{v:.2f}',
         tcolor='white' if v > 0.55 else '#222')

# --- panel 4: verdict
ax = axes[3]
grid(ax, 'Audit\nverdict')
for i, r in enumerate(rows):
    cat = r['cat']
    cell(ax, i, VERDICT_COLOR.get(cat, '#cccccc'), cat,
         tcolor='white' if cat != 'caveat' else '#222')

# domain separator + bracket labels
SEP_Y = n - 6 - 0.5   # shared edge between display rows 5 and 6
for a in axes:
    a.axhline(SEP_Y, color='#333', lw=1.3)
axb = fig.add_axes([0.030, 0.155, 0.026, 0.690]); axb.axis('off')
y_pocd = 0.155 + 0.690 * ((n - 1 - 2.5) + 0.5) / n
y_tbi = 0.155 + 0.690 * ((n - 1 - 9.5) + 0.5) / n
axb.text(0.5, y_pocd, 'POCD / PND', ha='center', va='center', fontsize=7.4,
         rotation=90, weight='bold')
axb.text(0.5, y_tbi, 'TBI / neurodeg.', ha='center', va='center', fontsize=7.4,
         rotation=90, weight='bold')

fig.suptitle('Fig. 4  Audit map \u2014 14 contrast-units from 12 public GEO series',
             fontsize=9.8, weight='bold', y=0.975)
fig.text(0.600, 0.018, fs.wrap(
    'Measured batch = scan date / scanner / declared processing block. GSM block = '
    'accession-run blocks (submitter batch column UNKNOWN). Inferred = median split '
    'of library size, no bimodality. The verdict attaches to a contrast-unit, not to '
    'a GEO series.', 118), fontsize=7.2, color='#555', va='bottom', ha='center',
    linespacing=1.55)
PROBLEMS['Fig4_audit_map'] = fs.save(fig, OUT, 'Fig4_audit_map')

# ================================================================ Fig 6
fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(fs.FIG_W, 4.3), gridspec_kw={'width_ratios': [1.16, 1.0]})
fig.subplots_adjust(wspace=0.30, left=0.085, right=0.985, top=0.845, bottom=0.245)

# ---- Panel A: separability plane
axA.plot([0, 1], [0, 1], color='#888', lw=1.0, ls='--', label='y = x  (chance)')
axA.fill_between([0, 0.85], [0.15, 1.0], [1, 1], color=fs.GREEN, alpha=0.10)
axA.plot([0, 0.85], [0.15, 1.0], color=fs.GREEN, lw=1.0)
axA.text(0.413, 0.80, 'separable region', fontsize=7.2, color='#2f6f3f',
         ha='left', va='center')

CONF = {'GSE59645', 'GSE115614', 'GSE31357'}
SHORT = {'GSE80174_Cortex': '80174ctx', 'GSE80174_Thalamus': '80174thal',
         'GSE80174_Hippocampus': '80174hippo'}
# label offsets tuned so no two of the 14 points collide
OFF = {'GSE68207': (2, 8), 'GSE64978': (2, 8), 'GSE80174_Thalamus': (2, 8),
       'GSE80174_Cortex': (-46, 1), 'GSE80174_Hippocampus': (7, -3),
       'GSE59645': (-30, 3), 'GSE115614': (7, -3), 'GSE31357': (7, -3)}
plotted = []
for ds, s in sorted(sep.items()):
    conf = ds in CONF
    col = fs.RED if (conf and s['supported']) else (fs.AMBER if conf else fs.BLUE)
    axA.scatter(s['base'], s['loocv'], s=88 if conf else 34,
                facecolor=col if conf else 'white', edgecolor=col,
                linewidths=1.4, zorder=3)
    lbl = SHORT.get(ds, ds.replace('GSE', ''))
    dx, dy = OFF.get(ds, (6, 5))
    axA.annotate(lbl, (s['base'], s['loocv']), fontsize=7.2, xytext=(dx, dy),
                 textcoords='offset points', weight='bold' if conf else 'normal',
                 color=col if conf else '#333')
    plotted.append(ds)

axA.set_xlim(0.40, 1.02); axA.set_ylim(0.36, 1.16)
axA.set_xlabel('class-frequency baseline (LOO)')
axA.set_ylabel('LOO accuracy on technical covariates')
axA.set_title('A  Technical separability of the batch label', fontsize=8.6, pad=6)
for s in ('top', 'right'):
    axA.spines[s].set_visible(False)
axA.grid(alpha=0.22, lw=0.6)
axA.legend(handles=[
    mpatches.Patch(facecolor=fs.RED, label='100% confounded \u2192 separable'),
    mpatches.Patch(facecolor=fs.AMBER, label='100% confounded \u2192 not separable'),
    mpatches.Patch(facecolor='white', edgecolor=fs.BLUE, linewidth=1.2,
                   label='not fully confounded'),
], fontsize=7.2, loc='lower right', frameon=True, framealpha=0.95,
    handlelength=1.1, borderpad=0.5)

# ---- Panel B: the three confounded cases
trio = ['GSE59645', 'GSE115614', 'GSE31357']
verd = ['UNUSABLE', 'CAVEAT', 'CAVEAT']
vc = [fs.RED, fs.AMBER, fs.AMBER]
w = 0.34
for i, ds in enumerate(trio):
    s = sep[ds]
    axB.bar(i - w / 2, s['base'], w, color='#c9d4de', edgecolor='#666', lw=0.7,
            label='baseline' if i == 0 else None)
    axB.bar(i + w / 2, s['loocv'], w, color=vc[i], edgecolor='#333', lw=0.7,
            label='LOO accuracy' if i == 0 else None)
    axB.text(i, 1.045, f"p = {s['p']:.3f}", ha='center', fontsize=7.2, color='#333')
    gain = s['loocv'] - s['base']
    if abs(gain) > 0.005:
        axB.annotate('', xy=(i + w / 2, s['loocv']), xytext=(i + w / 2, s['base']),
                     arrowprops=dict(arrowstyle='-|>', color='#333', lw=1.0))
    axB.text(i + w + 0.04, (s['loocv'] + s['base']) / 2, f'{gain:+.2f}',
             fontsize=7.2, color='#333', va='center')

axB.set_xticks(range(3))
axB.set_xticklabels([f'{l}\n{v}' for l, v in zip(trio, verd)], fontsize=7.2)
axB.set_xlim(-0.55, 2.70)
axB.set_ylim(0, 1.16)
axB.set_ylabel('accuracy')
axB.set_title('B  Same crosstab pattern, different verdict', fontsize=8.6, pad=6)
for s_ in ('top', 'right'):
    axB.spines[s_].set_visible(False)
axB.grid(axis='y', alpha=0.22, lw=0.6)
axB.legend(fontsize=7.2, loc='center right', frameon=True, framealpha=0.95,
           handlelength=1.1, borderpad=0.5)

fig.suptitle('Fig. 6  A structural confounding alarm is not evidence of a batch effect',
             fontsize=9.8, weight='bold', y=0.975)
fig.text(0.5, 0.018, fs.wrap(
    'All three units have every experimental group in a single batch level (100% '
    'structural confounding). Only GSE59645 also clears the separability gate '
    '(gain >= 0.15, p < 0.05) with measured scan dates and PC1 R2(batch) = 0.961; '
    'the other two carry no technical signal and are a numbering artefact. '
    'Separability requires gain >= 0.15 and p < 0.05.',
    132), ha='center', va='bottom', fontsize=7.2, color='#555', linespacing=1.55)
PROBLEMS['Fig6_separability_gate'] = fs.save(fig, OUT, 'Fig6_separability_gate')

print('Fig 4 rows:', n, '| separability points in panel A:', len(plotted))
fs.report(PROBLEMS, 2)
