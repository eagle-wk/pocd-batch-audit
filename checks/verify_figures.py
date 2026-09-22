"""Pre-submission figure check.

Verifies every figure in <repo>/figures/output against the size contract:
  * physical width  <= FIG_W (7.2 in)
  * raster dpi      >= 300
  * a vector .pdf   exists alongside each .png
  * smallest font in the generator >= 7 pt (reported by figstyle.save)

Optionally cross-checks the manuscript drafts (embedded figure paths resolve,
legend numbering ascends, no leftover placeholders). The drafts are author-side
material and are not shipped in this repository, so point MANUSCRIPT_DIR at
them to enable that section:

  python verify_figures.py
  MANUSCRIPT_DIR=/path/to/pocd_data python verify_figures.py
"""
import os
import re
import sys

from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))        # <repo>/checks
REPO = os.path.dirname(HERE)                            # <repo>
sys.path.insert(0, os.path.join(REPO, 'figures'))       # figstyle lives here
import figstyle as fs                                   # noqa: E402

FIGDIR = os.path.join(REPO, 'figures', 'output')
MDIR = os.environ.get('MANUSCRIPT_DIR', '')
MANUSCRIPTS = ['方法学审计稿件_draft_v0.5.html',
               '方法学审计稿件_draft_v0.5_中文版.html']

sys.stdout.reconfigure(encoding='utf-8')

problems = []
print(f'Contract: width <= {fs.FIG_W} in, dpi >= {fs.DPI}, PNG + PDF\n')
print(f'{"figure":34s} {"width":>7s} {"height":>7s} {"dpi":>6s}  pdf  status')

for f in sorted(os.listdir(FIGDIR)):
    if not f.endswith('.png'):
        continue
    im = Image.open(os.path.join(FIGDIR, f))
    dpi = im.info.get('dpi', (72, 72))[0]
    w, h = im.size[0] / dpi, im.size[1] / dpi
    has_pdf = os.path.exists(os.path.join(FIGDIR, f[:-4] + '.pdf'))
    why = []
    if w > fs.FIG_W:
        why.append(f'width {w:.2f} > {fs.FIG_W}')
    # PNG stores resolution as integer pixels-per-metre, so 300 dpi round-trips
    # as 11811 ppm = 299.9995 dpi. Compare with a tolerance, not exactly.
    if dpi < fs.DPI - 0.5:
        why.append(f'dpi {dpi:.1f} < {fs.DPI}')
    if not has_pdf:
        why.append('no pdf')
    print(f'{f:34s} {w:6.2f}in {h:6.2f}in {dpi:5.0f}  {"y" if has_pdf else "n":3s}  '
          f'{"OK" if not why else "FAIL: " + "; ".join(why)}')
    if why:
        problems.append((f, why))

# ---- manuscript cross-check: every embedded figure resolves, numbering ascends
if not MDIR:
    print('\nManuscript cross-check: SKIPPED')
    print('  (manuscript drafts are not part of this repository; set '
          'MANUSCRIPT_DIR to the folder holding the .html drafts to enable)')
else:
    print('\nManuscript cross-check:')
    for mf in MANUSCRIPTS:
        p = os.path.join(MDIR, mf)
        if not os.path.exists(p):
            print(f'  {mf}: MISSING')
            problems.append((mf, ['file missing']))
            continue
        t = open(p, encoding='utf-8').read()
        imgs = re.findall(r'<img src="(figures/[^"]+)"', t)
        missing = [x for x in imgs if not os.path.exists(os.path.join(MDIR, x))]
        # Every legend number appears twice: once as the in-body caption and once
        # in the "Figure legends" table. Check order of first appearance only.
        caps = [int(x) for x in re.findall(r'<b>(?:Fig\.|图) ([0-9]+)\.?</b>', t)]
        first_seen = list(dict.fromkeys(caps))
        ascending = first_seen == sorted(first_seen)
        print(f'  {mf}')
        print(f'    embedded: {[os.path.basename(x) for x in imgs]}')
        print(f'    all resolve: {not missing}'
              + (f'  MISSING {missing}' if missing else ''))
        print(f'    legend order {first_seen} ascending: {ascending}')
        print(f'    leftover placeholders: '
              f'{t.count("to be generated") + t.count("待出图")}')
        if missing:
            problems.append((mf, [f'unresolved {missing}']))
        if not ascending:
            problems.append((mf, [f'legend order not ascending: {first_seen}']))

print()
if problems:
    print(f'{len(problems)} problem(s):')
    for k, v in problems:
        print(f'  - {k}: {"; ".join(v)}')
    sys.exit(1)
print('all checks passed')
