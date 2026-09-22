"""Publication figure style for the methodological-audit manuscript.

Why this module exists
----------------------
The v1 figures were laid out for a screen (10.99-15.38 in wide). When scaled
into a journal text column they shrink by 40-55%, so an 8 pt label prints at
~4 pt. Journals specify a *physical* width and reject anything that is not
legible at it, so figures have to be designed at the printed size.

Contract enforced here
----------------------
  * FIG_W  = 7.2 in   hard ceiling on figure width (double-column full width).
              BMC / most Springer journals accept up to ~7.5 in; 7.2 leaves
              margin for the bbox we do not control.
  * MIN_PT = 7.0 pt   hard floor for every text element.
  * 300 dpi raster + vector PDF, both emitted with the same geometry.
  * Fonts embedded as TrueType (pdf.fonttype 42) so the PDF stays editable.

`save()` asserts the two hard limits before writing, so a layout regression
fails loudly here rather than silently at the proofs stage.

Bbox policy
-----------
`bbox_inches='tight'` is used but the *result* is measured. tight-bbox can
grow a figure if an artist sits outside the axes, which is exactly the failure
that produced the oversized v1 files. save() therefore reports the true
physical size of what landed on disk, not the requested figsize.
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --------------------------------------------------------------- constants
FIG_W = 7.2          # in — width ceiling
MIN_PT = 7.0         # pt — font floor
DPI = 300

# Publication palette (colour-blind safe; also legible in greyscale because
# the three verdict colours differ in luminance, not only in hue).
GREEN = '#3f8f4f'    # usable / proceed
AMBER = '#d8a41f'    # caveat / correct
RED = '#b3352f'      # unusable / stop
BLUE = '#2e5f8a'     # neutral data
VIOLET = '#6b5b95'   # normalize
GREY = '#666666'
LIGHT = '#eef2f7'
REDBG = '#fdf3f2'
GREENBG = '#e8f1e9'

# The five batch-provenance classes used by Fig 4, ordered strong -> weak so
# the ramp reads as "how much do we actually know about this batch variable".
PROV_COLOR = {'measured': '#1f5c8b',
              'GSM block': '#7fa8cf',
              'inferred': '#c3d3e0',
              'design-derived': '#9b7fb5',
              'none': '#ececec'}
PROV_TEXT = {'measured': 'measured', 'GSM block': 'GSM block',
             'inferred': 'inferred', 'design-derived': 'design', 'none': 'none'}

VERDICT_COLOR = {'usable': GREEN, 'caveat': AMBER, 'unusable': RED}


def use():
    """Apply the publication rcParams. Call once per figure script."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        # Arial is what most journals ask for; DejaVu Sans is the bundled
        # fallback so the script also runs on a machine without Arial.
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 8.5,
        'axes.titlesize': 9.5,
        'xtick.labelsize': 7.5,
        'ytick.labelsize': 7.5,
        'legend.fontsize': 7.5,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'lines.linewidth': 1.1,
        'figure.dpi': 120,          # on-screen only; savefig.dpi governs files
        'savefig.dpi': DPI,
        'savefig.facecolor': 'white',
        'pdf.fonttype': 42,         # embed TrueType, keep text editable
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })


def _min_fontsize(fig):
    """Smallest fontsize (pt) among all text artists that are visible."""
    sizes = []
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible():
            continue
        if not (t.get_text() or '').strip():
            continue
        sizes.append(t.get_fontsize())
    return min(sizes) if sizes else None


def _measure(path):
    """Return (width_in, height_in) of a raster as written to disk."""
    from PIL import Image
    im = Image.open(path)
    dpi = im.info.get('dpi', (DPI, DPI))[0] or DPI
    return im.size[0] / dpi, im.size[1] / dpi


def save(fig, outdir, stem, pad=0.02, allow_tall=None):
    """Write <stem>.png (300 dpi) and <stem>.pdf, then verify the contract.

    allow_tall : optional height ceiling in inches (most journals want a
                 figure no taller than the text block, ~9.5 in).
    """
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, stem + '.png')
    pdf = os.path.join(outdir, stem + '.pdf')

    fig.savefig(png, dpi=DPI, bbox_inches='tight', pad_inches=pad,
                facecolor='white')
    fig.savefig(pdf, bbox_inches='tight', pad_inches=pad, facecolor='white')

    w, h = _measure(png)
    fs = _min_fontsize(fig)
    problems = []
    if w > FIG_W + 1e-6:
        problems.append(f'width {w:.2f} in exceeds {FIG_W} in')
    if fs is not None and fs < MIN_PT - 1e-6:
        problems.append(f'smallest font {fs:.2f} pt below {MIN_PT} pt')
    if allow_tall and h > allow_tall + 1e-6:
        problems.append(f'height {h:.2f} in exceeds {allow_tall} in')

    tag = 'OK  ' if not problems else 'FAIL'
    print(f'  [{tag}] {stem:34s} {w:5.2f} x {h:5.2f} in  '
          f'min font {fs:.1f} pt  {os.path.getsize(png) // 1024} KB')
    for p in problems:
        print(f'         !! {p}')
    plt.close(fig)
    return problems


def wrap(text, width):
    """Hard-wrap a caption string; keeps footnotes from widening a figure."""
    import textwrap
    return '\n'.join(textwrap.wrap(text, width))


def report(problems_all, n_figures):
    """Print the summary and exit non-zero if the contract was violated."""
    bad = {k: v for k, v in problems_all.items() if v}
    print()
    if bad:
        print(f'{len(bad)}/{n_figures} figures violate the size contract:')
        for k, v in bad.items():
            print(f'  - {k}: {"; ".join(v)}')
        sys.exit(1)
    print(f'all {n_figures} figures pass: '
          f'width <= {FIG_W} in, minimum font >= {MIN_PT} pt, {DPI} dpi, PNG + PDF')
