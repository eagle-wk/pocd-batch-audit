# -*- coding: utf-8 -*-
"""把已下载的补充/processed 表达矩阵整理成 (expression_matrix.tsv + samplesheet.tsv),
供 batch_audit_generic.py 直接审计。只处理"列名自带分组标签"的数据集。

用法: python prep_downloaded_matrices.py --indir p0_1_downloads --outdir p0_1_downloads/prepared
"""
import argparse, csv, gzip, io, os, re, sys

def opener(p):
    return gzip.open(p, 'rt', encoding='utf-8', errors='replace') if p.endswith('.gz') else io.open(p, encoding='utf-8', errors='replace')

def read_table(path):
    with opener(path) as f:
        rows = [ln.rstrip('\n').rstrip('\r').split('\t') for ln in f]
    # 去掉尾部全空行
    while rows and all(c.strip() == '' for c in rows[-1]):
        rows.pop()
    return rows

CLASSIFY = [
    (r'_SHM_|_sham_|^sham|S[hS]am', 'control'),
    (r'_TRM_|_trauma_|^tbi|FPI|TBI|CCI', 'injured'),
]
def cls(name):
    for pat, lab in CLASSIFY:
        if re.search(pat, name):
            return lab
    return None

def build(name, path, gene_col, org, outdir, drop_unlabeled=True):
    rows = read_table(path)
    header = rows[0]
    # 跳过 value/status 交替列中的 status 列
    keep = [i for i, h in enumerate(header) if h.strip().lower() != 'status']
    header = [header[i] for i in keep]
    body = [[r[i] for i in keep] if len(r) > max(keep) else None for r in rows[1:]]
    body = [r for r in body if r is not None]
    gc = gene_col if isinstance(gene_col, int) else header.index(gene_col)
    sample_cols = [i for i in range(len(header)) if i != gc]
    samples, groups = [], []
    for i in sample_cols:
        lab = cls(header[i])
        if lab is None and drop_unlabeled:
            continue
        samples.append(header[i]); groups.append(lab or 'UNKNOWN')
    if len(samples) < 4:
        print(f'  x {name}: 可标注样本 {len(samples)} < 4,跳过')
        return None
    d = os.path.join(outdir, name)
    os.makedirs(d, exist_ok=True)
    with io.open(os.path.join(d, 'expression_matrix.tsv'), 'w', encoding='utf-8') as f:
        f.write('gene\t' + '\t'.join(samples) + '\n')
        for r in body:
            g = r[gc].strip()
            if not g:
                continue
            vals = []
            ok = True
            for i in sample_cols:
                if header[i] not in samples:
                    continue
                v = r[i] if i < len(r) else ''
                try:
                    vals.append(f'{float(v):.6f}')
                except Exception:
                    ok = False; break
            if ok:
                f.write(g + '\t' + '\t'.join(vals) + '\n')
    with io.open(os.path.join(d, 'samplesheet.tsv'), 'w', encoding='utf-8') as f:
        f.write('sample\tgroup\n')
        for s, gg in zip(samples, groups):
            f.write(f'{s}\t{gg}\n')
    print(f'  ✓ {name}: {len(samples)} 样本, 分组 {dict((g, groups.count(g)) for g in set(groups))}, org={org}')
    return {'name': name, 'dir': d, 'org': org}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', required=True)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    J = lambda f: os.path.join(a.indir, f)
    jobs = [
        ('GSE64978', 'GSE64978_tbi_fpkm_processed.txt.gz', 1, 'mouse'),
        ('GSE68207', 'GSE68207_tbi_ca1_fpkm_processed.txt.gz', 1, 'mouse'),
        ('GSE80174_Cortex', 'GSE80174_Cortex_DEseq2_raw_counts.txt.gz', 0, 'rat'),
        ('GSE80174_Hippocampus', 'GSE80174_Hippocampus_DEseq2_raw_counts.txt.gz', 0, 'rat'),
        ('GSE80174_Thalamus', 'GSE80174_Thalamus_DEseq2_raw_counts.txt.gz', 0, 'rat'),
    ]
    out = []
    for name, fn, gc, org in jobs:
        p = J(fn)
        if not os.path.exists(p):
            print(f'  - {name}: 文件缺失,跳过')
            continue
        r = build(name, p, gc, org, a.outdir)
        if r:
            out.append(r)
    print(f'\n整理完成: {len(out)} 个数据集 → {a.outdir}')
    with io.open(os.path.join(a.outdir, 'manifest.tsv'), 'w', encoding='utf-8') as f:
        f.write('name\tdir\torg\n')
        for r in out:
            f.write(f"{r['name']}\t{r['dir']}\t{r['org']}\n")

if __name__ == '__main__':
    main()
