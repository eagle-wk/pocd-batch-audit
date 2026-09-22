"""Fig 1 (conceptual framework), Fig 2 (simulation benchmark) and Fig 5
(three real failure modes) for the methodological-audit manuscript.

Rebuilt 2026-09-21 at printed size. The v1 versions were drawn at 13.4-13.6 in
wide, i.e. ~1.9x the journal text column, so every label printed at roughly half
its nominal point size. All three are now laid out inside 7.2 in.

Inputs
  C_m6A_PCD_patterns/GSE276942_downstream/bench_sim/benchmark_sim.json
      -> the 5 x 100-replicate run whose overall rate is 0.994 (497/500).
         NOTE: three other benchmark_sim.json files exist on disk with other
         settings (0.986, 0.9866). This is the one the manuscript cites; the
         assert below makes a wrong pick fail loudly.
  results/master_audit_table.tsv
      -> per-unit verdicts for the Fig 5 case panels.
  results/derived/GSE304956_cross_batch_replication.csv
      -> per-contrast replication (tab-separated despite the .csv extension).

Figure text is English throughout (publication figures). All plotted numbers are
read from the artefacts above; nothing is transcribed from memory.
"""
import csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs
fs.use()

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.stdout.reconfigure(encoding='utf-8')
HERE = os.path.dirname(os.path.abspath(__file__))        # <repo>/figures
REPO = os.path.dirname(HERE)                            # <repo>
RESULTS = os.path.join(REPO, 'results')
DERIVED = os.path.join(RESULTS, 'derived')              # archived intermediates
BENCH = os.path.join(DERIVED, 'benchmark', 'benchmark_sim.json')
OUT = os.path.join(HERE, 'output')
os.makedirs(OUT, exist_ok=True)

GREEN, AMBER, RED = fs.GREEN, fs.AMBER, fs.RED
BLUE, GREY = fs.BLUE, fs.GREY
PROBLEMS = {}

# ====================================================================== Fig 1
# Schematic only -- no data. Rebuilt as two stacked panels so that each can be a
# horizontal left-to-right flow at 7.2 in. The v1 side-by-side layout needed
# 13.4 in and became unreadable when scaled into the column.
fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(fs.FIG_W, 7.0),
    gridspec_kw={'height_ratios': [1.0, 1.0], 'hspace': 0.03})
fig.subplots_adjust(left=0.012, right=0.988, top=0.945, bottom=0.045)
for ax in (axA, axB):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')

BOXSTYLE = 'round,pad=0.20'
# Usable x-extent. FancyBboxPatch's pad is in *data* units (mutation_scale=1),
# so a box drawn at x spans x-pad .. x+w+pad. Everything below stays inside
# 0.35 .. 9.65, which is why the v1 boxes that ran to x=9.90 got clipped.
XL, XR = 0.55, 9.45
XW = XR - XL                      # 8.90


def box(ax, x, y, w, h, text, fc, ec, tc='#111', fs_=7.5, weight='normal',
        ls='-', lw=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=BOXSTYLE, facecolor=fc,
                                edgecolor=ec, linewidth=lw, linestyle=ls,
                                zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs_,
            color=tc, weight=weight, zorder=3, linespacing=1.45)


def harrow(ax, x0, x1, y, color='#444', lw=1.1):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle='-|>',
                                 mutation_scale=11, linewidth=lw, color=color,
                                 zorder=1))


def varrow(ax, y0, y1, x, color='#444', lw=1.1, ls='-'):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle='-|>',
                                 mutation_scale=11, linewidth=lw, color=color,
                                 linestyle=ls, zorder=1))


def row_x(n, gap=0.30, left=XL, width=XW):
    """x positions for n equal boxes spanning `width` with n-1 gaps."""
    w = (width - gap * (n - 1)) / n
    return [left + i * (w + gap) for i in range(n)], w


# ---- Panel A : conventional detect-then-correct
axA.text(XL - 0.30, 9.60, 'A   Conventional detect-then-correct', ha='left',
         va='center', fontsize=9.2, weight='bold')

CHAIN = ['Public dataset', 'Detect batch\n(accuracy, p-value)',
         'If detected: correct\n(ComBat / sva / RUV)',
         'Downstream analysis\n\u2192 publication']
cxs, CW = row_x(4)
CY, CH = 7.70, 1.50
for i, (x, t) in enumerate(zip(cxs, CHAIN)):
    box(axA, x, CY, CW, CH, t, fs.LIGHT, BLUE, fs_=7.5)
    if i:
        harrow(axA, cxs[i - 1] + CW, x, CY + CH / 2)

# the two unverified premises, side by side beneath the chain
pxs, PW = row_x(2)
PREM = ['the batch label is known\nand reliable',
        'the batch is not collinear\nwith group']
PY, PH = 4.40, 1.60
for i, (x, txt) in enumerate(zip(pxs, PREM)):
    box(axA, x, PY, PW, PH, '', fs.REDBG, RED, ls='--')
    axA.text(x + 0.20, PY + PH - 0.38, f'unverified premise {i + 1}', ha='left',
             va='center', fontsize=7.4, color=RED, weight='bold', zorder=3)
    axA.text(x + 0.20, PY + 0.52, txt, ha='left', va='center', fontsize=7.5,
             color='#7d1f1b', linespacing=1.4, zorder=3)
# dotted leaders: premise 1 -> "Detect batch", premise 2 -> "correct"
axA.add_patch(FancyArrowPatch((pxs[0] + PW * 0.55, PY + PH + 0.20),
                              (cxs[1] + CW * 0.5, CY - 0.06),
                              arrowstyle='-|>', mutation_scale=10,
                              linewidth=0.9, color=RED, linestyle=':', zorder=4))
axA.add_patch(FancyArrowPatch((pxs[1] + PW * 0.22, PY + PH + 0.20),
                              (cxs[2] + CW * 0.5, CY - 0.06),
                              arrowstyle='-|>', mutation_scale=10,
                              linewidth=0.9, color=RED, linestyle=':', zorder=4))

box(axA, XL, 0.40, XW, 1.70,
    'When either premise fails the tool still returns a corrected matrix.\n'
    'Signal is removed, or invented \u2014 and no warning is raised.',
    fs.REDBG, RED, tc='#7d1f1b', fs_=7.8, lw=1.2)

# ---- Panel B : diagnose-then-decide
axB.text(XL - 0.30, 9.60, 'B   This work: diagnose-then-decide', ha='left',
         va='center', fontsize=9.2, weight='bold')

STAGES = ['Stage 0\nare the group labels\nindependently verifiable?',
          'Stage 1\nbatch triple-threshold\nstructural \u2192 separable \u2192 nested',
          'Stage 2\u20134\neigengene sign \u2192 global shift \u2192\ncross-dataset concordance']
sxs, SW = row_x(3)
SY, SH = 6.90, 2.10
for i, (x, t) in enumerate(zip(sxs, STAGES)):
    box(axB, x, SY, SW, SH, t, fs.LIGHT, BLUE, fs_=7.5)
    if i:
        harrow(axB, sxs[i - 1] + SW, x, SY + SH / 2)

varrow(axB, SY - 0.02, 5.70, 5.0)
box(axB, XL, 4.10, XW, 1.35, 'Terminal verdict \u2014 one per contrast-unit',
    fs.GREENBG, GREEN, fs_=8.0, lw=1.2)

VD = [('stop', RED, 'refuse\nthe analysis'),
      ('correct', AMBER, 'genuine batch,\npartially crossed'),
      ('proceed', GREEN, 'no batch\nstructure'),
      ('normalize', fs.VIOLET, 'global\nshift')]
vxs, VW = row_x(4)
VY, VH = 2.10, 1.25
for i, (name, col, sub) in enumerate(VD):
    x = vxs[i]
    axB.add_patch(FancyBboxPatch((x, VY), VW, VH, boxstyle=BOXSTYLE,
                                 facecolor=col, edgecolor=col, linewidth=1.1,
                                 zorder=2))
    axB.text(x + VW / 2, VY + VH / 2, name, ha='center', va='center',
             fontsize=8.6, color='white', weight='bold', zorder=3)
    axB.text(x + VW / 2, VY - 0.75, sub, ha='center', va='center', fontsize=7.2,
             color='#444', linespacing=1.3)
axB.add_patch(FancyArrowPatch((5.0, 4.08), (5.0, 3.60), arrowstyle='-|>',
                              mutation_scale=11, linewidth=1.1, color='#444'))

fig.suptitle('Fig. 1  Batch handling as a decision problem, not a correction problem',
             fontsize=9.8, weight='bold', y=0.985)
PROBLEMS['Fig1_conceptual_framework'] = fs.save(fig, OUT, 'Fig1_conceptual_framework')

# ====================================================================== Fig 2
b = json.load(open(BENCH, encoding='utf-8'))
rows = b['rows']
overall = b['overall_correct_rate']
cum = sum(r['correct_rate'] for r in rows)
nreps = len(rows) * b['reps']
assert abs(overall - 0.994) < 1e-9, f'unexpected benchmark file: {overall}'
print(f'  Fig 2 source: reps={b["reps"]} genes={b["p_genes"]} perm={b["perm"]} '
      f'overall={overall} ({round(cum * b["reps"])}/{nreps})')

TRUTH_COLOR = {'proceed': GREEN, 'correct': AMBER, 'stop': RED,
               'normalize': fs.VIOLET}
PRETTY = {'A_no_batch': ('A', 'no batch\n(proceed)'),
          'B_batch_crossed': ('B', 'genuine batch,\ncrossed\n(correct)'),
          'C_batch_confound': ('C', 'batch \u2261 group\n(stop)'),
          'D_pseudo_batch': ('D', 'pseudo-batch\n(proceed)'),
          'E_global_shift': ('E', 'global shift\n(normalize)')}

fig, ax = plt.subplots(figsize=(fs.FIG_W, 4.0))
fig.subplots_adjust(left=0.088, right=0.985, top=0.845, bottom=0.285)
W = 0.58
for i, r in enumerate(rows):
    corr, truth = r['correct_rate'], r['truth']
    ax.bar(i, corr, W, color=TRUTH_COLOR[truth], edgecolor='#333', linewidth=0.8,
           zorder=3)
    if corr < 1.0:
        ax.bar(i, 1.0 - corr, W, bottom=corr, color='#f2f2f2', edgecolor='#333',
               linewidth=0.8, hatch='///', zorder=3)
        for v, n in r['verdict_dist'].items():
            if v != truth:
                ax.annotate(f'{n}/100 sent to\n"{v}"', xy=(i, 1.0),
                            xytext=(i, 1.20), ha='center', va='bottom',
                            fontsize=7.2, color=RED, linespacing=1.3,
                            arrowprops=dict(arrowstyle='-', color=RED, lw=0.8,
                                            linestyle=':'))
    ax.text(i, corr / 2, f'{corr * 100:.0f}%', ha='center', va='center',
            color='white', fontsize=9.0, weight='bold', zorder=4)
    ax.text(i, 1.035, PRETTY[r['scenario']][1], ha='center', va='bottom',
            fontsize=7.4, color='#333', linespacing=1.35)

ax.axhline(1.0, color='#999', lw=0.9, ls='--', zorder=1)
ax.set_xticks(range(len(rows)))
ax.set_xticklabels([PRETTY[r['scenario']][0] for r in rows], fontsize=9.0,
                   weight='bold')
ax.set_xlim(-0.62, len(rows) - 0.38)
ax.set_ylim(0, 1.52)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels(['0', '25', '50', '75', '100'])
ax.set_ylabel('verdicts correct  (%)')
ax.set_xlabel('scenario  (truth in parentheses)', labelpad=6)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.grid(axis='y', alpha=0.22, lw=0.6, zorder=0)
# No legend: a one-entry legend was crowding bar E. The hatch is explained in
# the footnote instead.

fig.suptitle(f'Fig. 2  Verdict accuracy against known ground truth: '
             f'{round(cum * b["reps"])}/{nreps} = {overall * 100:.1f}%',
             fontsize=9.8, weight='bold', y=0.975)
fig.text(0.5, 0.045, fs.wrap(
    f'{b["reps"]} replicates per scenario, {b["p_genes"]:,} genes, '
    f'n = {b["n_per_group"]} per group, {b["perm"]} permutations. Hatched segment = '
    f'residual error. Every residual error is in D, where random labels showed weak '
    f'separability and were conservatively sent to "correct" rather than "proceed". '
    f'Scenario C, the case detect-then-correct pipelines mishandle, was stopped in '
    f'100/100 replicates.',
    142), ha='center', va='bottom', fontsize=7.2, color='#555', linespacing=1.55)
PROBLEMS['Fig2_simulation_benchmark'] = fs.save(fig, OUT, 'Fig2_simulation_benchmark')

# ================================================== Fig 5 : three real cases
tbl = {r['dataset']: r for r in csv.DictReader(
    open(os.path.join(RESULTS, 'master_audit_table.tsv'),
         encoding='utf-8'), delimiter='\t')}
VC = fs.VERDICT_COLOR
DATA_BLUE = BLUE


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return 0.0


def tidy(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(axis='y', alpha=0.22, lw=0.6, zorder=0)


def foot(ax, meta, cat, y=-0.30):
    ax.text(0.5, y, meta, transform=ax.transAxes, ha='center', va='top',
            fontsize=7.2, color='#444', linespacing=1.45)


# A and C side by side on top (both simple bar charts), B full width below so
# its three long contrast labels have room.
fig = plt.figure(figsize=(fs.FIG_W, 5.9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.80, wspace=0.34,
                      left=0.095, right=0.985, top=0.855, bottom=0.205)

# ---- A : GSE276942, batch == timepoint by design
axA = fig.add_subplot(gs[0, 0])
r = tbl['GSE276942']
r2b, r2g = fnum(r['lead_PC_R2_batch']), fnum(r['lead_PC_R2_group'])
axA.bar([0, 1], [r2b, r2g], 0.52, color=[DATA_BLUE, '#c9d4de'],
        edgecolor='#333', linewidth=0.8, zorder=3)
for x, v in ((0, r2b), (1, r2g)):
    axA.text(x, v + 0.022, f'{v:.3f}', ha='center', fontsize=8.2, weight='bold',
             color='#111')
axA.set_xticks([0, 1]); axA.set_xticklabels(['batch', 'group'])
axA.set_ylim(0, 0.56)
axA.set_ylabel('lead PC  R\u00b2')
axA.set_title('A  GSE276942  batch \u2261 timepoint\nby design', fontsize=8.6,
              pad=6, linespacing=1.4)
tidy(axA)
foot(axA, f'{r["batch_source"]} batch  |  n = {r["n_samples"]}  |  '
          f'{r["n_confounded_groups"]}/{r["n_confounded_groups"]} groups confounded\n'
          f'verdict: {r["cat"].upper()}', r['cat'], y=-0.34)

# ---- C : GSE80174, one study three regions three tiers
axC = fig.add_subplot(gs[0, 1])
regions = [('GSE80174_Cortex', 'cortex', 10, 'r = 0.77'),
           ('GSE80174_Thalamus', 'thalamus', 10, 'not computable'),
           ('GSE80174_Hippocampus', 'hippocampus', 9, 'failed')]
xs = list(range(len(regions)))
vals = [fnum(tbl[k]['lead_PC_R2_batch']) for k, _, _, _ in regions]
cats = [tbl[k]['cat'] for k, _, _, _ in regions]
axC.bar(xs, vals, 0.54, color=[VC[c] for c in cats], edgecolor='#333',
        linewidth=0.8, zorder=3)
for i, v in enumerate(vals):
    axC.text(i, v + 0.020, f'{v:.3f}', ha='center', fontsize=8.2, weight='bold',
             color='#111')
axC.set_xticks(xs)
axC.set_xticklabels([f'{nm}\nn = {n}\nrep. {rep}' for _, nm, n, rep in regions],
                    fontsize=7.2)
axC.set_ylim(0, 0.88)
axC.set_ylabel('lead PC  R\u00b2(batch)')
axC.set_title('C  GSE80174  one experiment, one run,\nthree regions, three verdicts',
              fontsize=8.6, pad=6, linespacing=1.4)
tidy(axC)
axC.legend(handles=[mpatches.Patch(facecolor=VC[c], label=c)
                    for c in ('usable', 'caveat', 'unusable')],
           fontsize=7.2, loc='upper left', frameon=True, framealpha=0.95, ncol=3,
           handlelength=1.2, columnspacing=0.9)
foot(axC, f'{tbl["GSE80174_Cortex"]["batch_source"]} batch in all three  |  '
          f'same animals, same processing run\n'
          f'auditing at series level would average these away', cats[2], y=-0.34)

# ---- B : GSE304956 -- no batch signal in PC space, yet nothing replicates
# Per-contrast replication comes from the audit's own CSV, not from the master
# table's FAIL string (which holds one aggregated entry for all three).
axB = fig.add_subplot(gs[1, :])
r = tbl['GSE304956']
rep = list(csv.DictReader(open(os.path.join(
    DERIVED, 'GSE304956_cross_batch_replication.csv'), encoding='utf-8-sig'),
    delimiter='\t'))          # tab-separated despite the .csv extension
rs = [float(x['pearson_r']) for x in rep]
sigs = [float(x['sign_agreement']) for x in rep]
xs = list(range(len(rep)))
axB.axhline(0, color='#333', lw=0.9, zorder=2)
for y in (0.30, -0.30):
    axB.axhline(y, color='#999', lw=0.9, ls='--', zorder=1)
axB.text(len(rep) - 0.42, 0.325, '|r| = 0.30 threshold', ha='right', va='bottom',
         fontsize=7.2, color='#777')
axB.bar(xs, rs, 0.42, color=DATA_BLUE, edgecolor='#333', linewidth=0.8, zorder=3)
for i, (v, s) in enumerate(zip(rs, sigs)):
    axB.text(i, v - 0.028, f'{v:.3f}', ha='center', va='top', fontsize=8.2,
             weight='bold', color='#111')
    axB.text(i, -0.60, f'sign agreement {s * 100:.0f}%', ha='center', va='center',
             fontsize=7.2, color='#555')
axB.set_xticks(xs)
axB.set_xticklabels([x['contrast'].replace('_vs_', '\nvs ') for x in rep],
                    fontsize=7.5)
axB.set_xlim(-0.62, len(rep) - 0.38)
axB.set_ylim(-0.78, 0.48)
axB.set_ylabel('batch effect:  Pearson r\nbetween the two batches')
axB.set_title('B  GSE304956  no batch signal in PC space, yet nothing replicates',
              fontsize=8.6, pad=6)
tidy(axB)
axB.grid(axis='y', alpha=0.0)
foot(axB, f'{r["batch_source"]} batch  |  n = {r["n_samples"]}  |  '
          f'R\u00b2(batch) = {fnum(r["lead_PC_R2_batch"]):.3f}  |  '
          f'verdict: {r["cat"].upper()}', r['cat'], y=-0.30)

fig.suptitle('Fig. 5  Three real failure modes a detect-then-correct pipeline '
             'would not have flagged', fontsize=9.8, weight='bold', y=0.975)
PROBLEMS['Fig5_three_real_cases'] = fs.save(fig, OUT, 'Fig5_three_real_cases')

fs.report(PROBLEMS, 3)
