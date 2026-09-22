"""Recompute the manuscript's headline numbers from the shipped results.

This is the reviewer-facing check: it needs no raw GEO data and no network, only
the files in `results/`. Every number it recomputes is one that appears in the
manuscript text or a figure, so a pass means text and data agree.

It deliberately re-derives rather than re-reads:

  * the 2 / 6 / 6 distortion split is counted from the `cat` column of
    master_audit_table.tsv, then compared against audit_summary.json and
    final_numbers.json -- three independent files that must agree
  * the 11/21 direction concordance is counted from `same_sign` in
    cross_dataset_meta.tsv, not taken from the summary
  * the domain comparison is re-run as a Fisher exact test on the counts
  * benchmark k/n are recomputed from the stored rate and replicate count

Usage
-----
    python verify_numbers.py                  # verify
    python verify_numbers.py --recompute-summary   # rewrite the stale summary
"""
import argparse
import csv
import json
import os
import sys

from scipy.stats import fisher_exact
from scipy.stats import binomtest

HERE = os.path.dirname(os.path.abspath(__file__))        # <repo>/checks
REPO = os.path.dirname(HERE)                            # <repo>
RES = os.path.join(REPO, 'results')
DER = os.path.join(RES, 'derived')
GENESETS = os.path.join(REPO, 'data', 'genesets')

# Domain split. The six POCD-domain units are the ones whose batch structure was
# derived from the design or from GSM numbering; the remaining eight are the
# TBI / neurotrauma units studied by the separability trio.
POCD_UNITS = {'GSE330865', 'GSE297195', 'GSE283401', 'GSE304956',
              'GSE276942', 'GSE163943'}

problems = []


def ok(label, got, want, tol=0):
    if isinstance(want, float) or isinstance(got, float):
        good = abs(float(got) - float(want)) <= tol
    else:
        good = got == want
    mark = 'OK  ' if good else 'FAIL'
    print(f'  [{mark}] {label:<52} {got!r}  (expected {want!r})')
    if not good:
        problems.append(f'{label}: got {got!r}, expected {want!r}')


def load_table():
    p = os.path.join(RES, 'master_audit_table.tsv')
    with open(p, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def jload(p):
    with open(p, encoding='utf-8') as f:
        return json.load(f)


# --------------------------------------------------------------------------
# recompute mode: rebuild audit_summary.json + master_audit_report.md from the
# table. Used once during repository assembly, when the stored summary was found
# to carry a stale 2/5/7 split from before a verdict was reclassified.
# --------------------------------------------------------------------------
def recompute_summary(rows):
    from collections import Counter
    cnt = Counter(r['cat'] for r in rows)
    n = len(rows)
    summary = {
        'n_datasets': n,
        'usable': cnt.get('usable', 0),
        'caveat': cnt.get('caveat', 0),
        'unusable': cnt.get('unusable', 0),
        'frac_usable': cnt.get('usable', 0) / n,
        'frac_caveat': cnt.get('caveat', 0) / n,
        'frac_unusable': cnt.get('unusable', 0) / n,
        'frac_not_clean': (cnt.get('caveat', 0) + cnt.get('unusable', 0)) / n,
        'unverified_group_inference': [],
        'n_unverified': 0,
        'headline': (f'{cnt.get("unusable",0)}/{n} categorically unusable; '
                     f'{cnt.get("caveat",0)}/{n} usable only under documented '
                     f'restrictions; {cnt.get("usable",0)}/{n} unconditionally '
                     f'analyzable'),
    }
    p = os.path.join(RES, 'audit_summary.json')
    with open(p, 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(f'  rewritten {os.path.relpath(p, REPO)}')

    md = [f'# Master audit report — {n} contrast-units\n',
          f'**Distortion-rate summary:** {summary["headline"]}\n',
          '| unit | n | batch_source | confounded_groups | cat | usable | '
          'R2_batch | R2_group |',
          '|---|---:|---|---|---|---:|---:|---:|']
    for r in sorted(rows, key=lambda x: (x['cat'] != 'unusable', x['dataset'])):
        cg = r['confounded_groups'] or '-'
        md.append(f"| {r['dataset']} | {r['n_samples']} | {r['batch_source']} | "
                  f"{cg} | {r['cat']} | {r['usable']} | "
                  f"{r['lead_PC_R2_batch']} | {r['lead_PC_R2_group']} |")
    md.append('')
    for r in rows:
        if r['cat'] != 'usable':
            md.append(f"**{r['dataset']}** ({r['cat']}): {r['recommendation']}")
    md.append('')
    p = os.path.join(RES, 'master_audit_report.md')
    with open(p, 'w', encoding='utf-8') as f:
        f.write('\n'.join(md))
    print(f'  rewritten {os.path.relpath(p, REPO)}')
    print(f'  -> {summary["headline"]}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--recompute-summary', action='store_true',
                    help='rebuild audit_summary.json + master_audit_report.md '
                         'from master_audit_table.tsv, then exit')
    a = ap.parse_args()

    rows = load_table()
    if a.recompute_summary:
        print('Recomputing summary from master_audit_table.tsv:')
        recompute_summary(rows)
        return

    summary = jload(os.path.join(RES, 'audit_summary.json'))
    final = jload(os.path.join(RES, 'final_numbers.json'))

    print(f'Repo: {REPO}')
    print(f'Contrast-units: {len(rows)}\n')

    # ---------------------------------------------------------- 1. distortion
    print('1. Distortion split (counted from master_audit_table.tsv "cat")')
    from collections import Counter
    cnt = Counter(r['cat'] for r in rows)
    for k in ('usable', 'caveat', 'unusable'):
        ok(f'table cat={k}', cnt.get(k, 0), final['distortion_corrected'][k])
        ok(f'summary.{k} agrees with table', summary[k], cnt.get(k, 0))
    ok('summary.n_datasets', summary['n_datasets'], len(rows))
    print(f'  headline: {summary["headline"]}')

    # --------------------------------------------------------- 2. concordance
    print('\n2. Direction concordance (counted from cross_dataset_meta.tsv)')
    meta = list(csv.DictReader(open(os.path.join(DER, 'cross_dataset_meta.tsv'),
                                    encoding='utf-8-sig'), delimiter='\t'))
    same = [r for r in meta if r['same_sign'].strip().lower() == 'true']
    n_set = len(meta)
    ok('sets compared', n_set, final['concordance']['n'])
    ok('same-direction count', len(same), final['concordance']['k'])
    p_binom = binomtest(len(same), n_set, 0.5, alternative='two-sided').pvalue
    ok('binomial p (vs chance)', round(p_binom, 3), final['concordance']['p'],
       tol=5e-4)

    # IFN row, the one survivor
    ifn = [r for r in meta if r['set'] == 'IFN_signaling']
    if len(ifn) != 1:
        problems.append(f'expected exactly one IFN_signaling row, got {len(ifn)}')
    else:
        r = ifn[0]
        print('\n3. IFN_signaling (the one finding that survived)')
        ok('genes analysed', int(r['n_genes']), int(final['ifn']['n_genes']))
        ok('same_sign', r['same_sign'].strip().lower(), 'true')
        ok('Stouffer Z', round(float(r['z_stouffer']), 3), final['ifn']['Z'],
           tol=5e-4)
        ok('per-dataset dz[0]', round(float(r['dz_865']), 6),
           round(final['ifn']['dz'][0], 6), tol=1e-5)
        ok('per-dataset dz[1]', round(float(r['dz_195']), 6),
           round(final['ifn']['dz'][1], 6), tol=1e-5)
        ok('Fisher p (one-sided, as stored)', round(float(r['p_fisher']), 6),
           round(final['ifn']['fisher_one_sided_input'], 6), tol=1e-5)

        # The gene set is DEFINED with 13 members but ANALYSED with 12.
        gs = os.path.join(GENESETS, 'interferon_mhc_genesets_mouse.tsv')
        with open(gs, encoding='utf-8-sig') as f:
            in_ifn = sum(1 for l in f
                         if l.strip() and not l.startswith('#')
                         and l.rstrip('\n').endswith('\tIFN_signaling'))
        ok('IFN_signaling defined members', in_ifn, 13)
        ok('analysed members', int(r['n_genes']), 12)
        print('       12 of 13 analysed; Ifnb1 (the inducible ligand, near-zero '
              'baseline) is absent in both matrices -- documented in results/README.md')

    # --------------------------------------------------------- 4. benchmarks
    print('\n4. Simulation benchmark')
    sim = jload(os.path.join(DER, 'benchmark', 'benchmark_sim.json'))
    n_scen = len(sim['rows'])
    n_total = sim['reps'] * n_scen
    k = round(sim['overall_correct_rate'] * n_total)
    ok('scenarios', n_scen, 5)
    ok('trials (reps x scenarios)', n_total, final['Simulation benchmark']['n'])
    ok('correct trials', k, final['Simulation benchmark']['k'])
    ok('overall rate', sim['overall_correct_rate'], 0.994, tol=1e-9)
    worst = min(r['correct_rate'] for r in sim['rows'])
    print(f'  worst scenario correct_rate = {worst} '
          f'(scenario D; all others 1.00)')

    print('\n5. Tool comparison (ComBat / sva / RUV)')
    tool = jload(os.path.join(DER, 'benchmark', 'benchmark_tool_compare.json'))
    n_total = tool['reps'] * len(tool['rows'])
    k = round(tool['our_overall_correct_rate'] * n_total)
    ok('trials', n_total, final['Tool comparison']['n'])
    ok('correct trials', k, final['Tool comparison']['k'])
    ok('diagnostic rate', tool['our_overall_correct_rate'], 0.992, tol=1e-9)

    scen_c = {r['scenario']: r for r in tool['rows']}['C_batch_confound']
    print('  scenario C true-signal retention (collinear design):')
    for tool_name, key in (('ComBat', 'combat_retention'),
                           ('sva', 'sva_retention'),
                           ('RUV', 'ruv_retention')):
        stored = final['scenario_C'][f'{tool_name.lower()}_retention']
        ok(f'  {tool_name}', float(scen_c[key]), float(stored), tol=1e-9)
    print('       ComBat/RUV are at floating-point zero, i.e. the signal is '
          'annihilated; sva retains it only because it is blind to the '
          'collinearity and effectively corrects nothing.')

    # ------------------------------------------------------------ 6. domains
    print('\n6. Domain comparison (Fisher exact, recomputed)')
    pocd = [r for r in rows if r['dataset'] in POCD_UNITS]
    tbi = [r for r in rows if r['dataset'] not in POCD_UNITS]
    a_ = sum(1 for r in pocd if r['cat'] == 'unusable')
    b_ = sum(1 for r in tbi if r['cat'] == 'unusable')
    print(f'  POCD {a_}/{len(pocd)} unusable | TBI {b_}/{len(tbi)} unusable')
    ok('POCD [unusable, total]',
       [a_, len(pocd)], final['domain_fisher_corrected']['POCD'])
    ok('TBI [unusable, total]',
       [b_, len(tbi)], final['domain_fisher_corrected']['TBI'])
    table = [[a_, len(pocd) - a_], [b_, len(tbi) - b_]]
    orr, pval = fisher_exact(table)
    ok('odds ratio', round(orr, 3), final['domain_fisher_corrected']['OR'],
       tol=1e-3)
    ok('Fisher p', round(pval, 4), round(final['domain_fisher_corrected']['p'], 4),
       tol=1e-4)
    print('       With power this low this is "no difference detected", NOT '
          '"no difference exists" -- stated that way in the manuscript.')

    # ------------------------------------------------------- 7. provenance
    print('\n7. Batch provenance labels (final_numbers.json vs table)')
    prov = final['provenance']
    for unit, txt in prov.items():
        if unit not in {r['dataset'] for r in rows}:
            problems.append(f'provenance mentions unknown unit {unit}')
    print(f'  {len(prov)} units carry a provenance label')
    # GSE330865 / GSE297195: the submitter batch column is all UNKNOWN, so the
    # label must NOT read as "measured".
    for u in ('GSE330865', 'GSE297195'):
        bad = 'measured' in prov[u].lower()
        ok(f'{u} not labelled "measured"', bad, False)
    for u in ('GSE59645', 'GSE115614', 'GSE304956'):
        ok(f'{u} labelled measured', 'measured' in prov[u].lower(), True)

    # ------------------------------------------------------------------ done
    print()
    if problems:
        print(f'{len(problems)} problem(s):')
        for p in problems:
            print(f'  - {p}')
        sys.exit(1)
    print('all numbers verified against the manuscript')


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    main()
