# -*- coding: utf-8 -*-
"""
基因集逐基因 Δz 明细(可复现版)
================================================================================
为什么需要单独一个脚本
----------------------
`geneset_shift_diagnostic.py` 只输出集合层面的平均 Δz 与经验 p。但当某个集合
显著时,审稿人(和自己)一定会追问一句:

    "这是 12 个基因一起吃劲,还是某一个基因在拉着走?"

这个问题必须在论文里回答,否则单基因驱动的大效应会被当成通路效应。本脚本

  [1] 逐基因列出分层内组间 Δz(各层等权平均),原始矩阵与 QN 矩阵各一份
  [2] 给出"同向基因数 / 可检出基因数"与方向一致性比例
  [3] 标出反例基因(与集合主方向相反的),便于在 limitations 里点名

注意方向不要看错
----------------
Δz 的符号是"case − ctrl"经 z 标准化后的差值:
  Δz < 0 → 手术组相对 sham 组**下调**
  Δz > 0 → 上调
与 `geneset_shift_diagnostic.py` 的 mean_dz 完全同口径(都是先按基因 z 标准化
再求差),因此本脚本的单基因均值应当≈该脚本的 mean_dz,可互相验算。

用法
----
python geneset_per_gene_dz.py --matrix M.tsv --samples SS.tsv \
    --genesets G.tsv --set IFN_signaling \
    --group-col group --stratum-col age --case surgery --ctrl sham \
    --tag GSE297195_24h --out tsv/明细.tsv
"""
import argparse, csv, io, json, os, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--matrix', required=True)
ap.add_argument('--samples', required=True)
ap.add_argument('--genesets', required=True, help='TSV: gene<TAB>set')
ap.add_argument('--set', dest='setname', required=True, help='要展开的集合名')
ap.add_argument('--group-col', default='group')
ap.add_argument('--stratum-col', default='age')
ap.add_argument('--sample-col', default='sample')
ap.add_argument('--case', required=True)
ap.add_argument('--ctrl', required=True)
ap.add_argument('--qn', action='store_true', help='同时算分位数归一化后的 Δz')
ap.add_argument('--tag', default='')
ap.add_argument('--out', default='')
args = ap.parse_args()

# ---------------- 读数据 ----------------
srows = list(csv.DictReader(io.open(args.samples, encoding='utf-8'), delimiter='\t'))
samples = [r[args.sample_col].strip() for r in srows]
grp = [r[args.group_col].strip() for r in srows]
stg = [r[args.stratum_col].strip() for r in srows]

rows = list(csv.reader(io.open(args.matrix, encoding='utf-8'), delimiter='\t'))
hdr = rows[0][1:]
genes = [r[0] for r in rows[1:]]
X = np.array([[float(x) for x in r[1:]] for r in rows[1:]], dtype=float)
col = {s: j for j, s in enumerate(hdr)}
X = X[:, [col[s] for s in samples]]
print(f'矩阵: {X.shape[0]} 基因 x {X.shape[1]} 样本')

gi = {g: i for i, g in enumerate(genes)}
upper = {g.upper(): g for g in genes}

# ---------------- 读基因集 ----------------
want = []
with io.open(args.genesets, encoding='utf-8') as f:
    for line in f:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 2 and parts[1].strip() == args.setname:
            want.append(parts[0].strip())
want = sorted(set(want))
print(f'集合 {args.setname}: {len(want)} 个基因(表内)')

# ---------------- 分层 ----------------
strata = sorted(set(stg))
masks = []
for s in strata:
    ci = [i for i in range(len(samples)) if stg[i] == s and grp[i] == args.ctrl]
    ki = [i for i in range(len(samples)) if stg[i] == s and grp[i] == args.case]
    if len(ci) >= 2 and len(ki) >= 2:
        masks.append((ci, ki, s, len(ci), len(ki)))
print(f'层: {[(m[2], f"{m[3]}v{m[4]}") for m in masks]}')
if not masks:
    sys.exit('没有可用的分层(每组至少 2 个样本),退出')


def zscore_rows(M):
    return (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-12)


def per_gene_dz(M):
    Z = zscore_rows(M)
    acc, wsum = np.zeros(Z.shape[0]), 0.0
    for ci, ki, s, nc, nk in masks:
        acc += Z[:, ki].mean(axis=1) - Z[:, ci].mean(axis=1)
        wsum += 1
    return acc / max(wsum, 1)


def quantile_normalize(M):
    ranks = np.argsort(np.argsort(M, axis=0), axis=0)
    ref = np.sort(M, axis=0).mean(axis=1)
    out = np.empty_like(M)
    for j in range(M.shape[1]):
        out[:, j] = ref[ranks[:, j]]
    return out


dz_raw = per_gene_dz(X)
dz_qn = per_gene_dz(quantile_normalize(X)) if args.qn else None

# ---------------- 逐基因 ----------------
present, absent = [], []
for g in want:
    i = gi.get(g)
    if i is None and g.upper() in upper:
        i = gi[upper[g.upper()]]
    (present if i is not None else absent).append((g, i))

recs = []
for g, i in present:
    recs.append({
        'gene': genes[i], 'query': g,
        'dz_raw': float(dz_raw[i]),
        'dz_qn': float(dz_qn[i]) if dz_qn is not None else float('nan'),
    })
recs.sort(key=lambda r: r['dz_raw'])

v_raw = np.array([r['dz_raw'] for r in recs])
v_qn = np.array([r['dz_qn'] for r in recs]) if dz_qn is not None else None
n_neg = int((v_raw < 0).sum())
main_sign = -1 if n_neg * 2 >= len(v_raw) else 1

print(f'\n{"="*78}')
print(f'【{args.setname}】逐基因分层内 Δz   tag={args.tag}')
print(f'{"="*78}')
print(f'{"基因":<14}{"Δz(原始)":>11}{"Δz(QN)":>11}   方向')
for r in recs:
    d = r['dz_raw']
    flag = '一致' if (d < 0) == (main_sign < 0) else '**反例**'
    qs = f"{r['dz_qn']:>+11.3f}" if dz_qn is not None else f"{'-':>11}"
    print(f"{r['gene']:<14}{d:>+11.3f}{qs}   {flag}")

print(f'\n可检出 {len(recs)}/{len(want)} 个基因'
      + (f"(未检出: {', '.join(g for g, _ in absent)})" if absent else ''))
print(f'同向 {max(n_neg, len(recs)-n_neg)}/{len(recs)} = '
      f'{100*max(n_neg, len(recs)-n_neg)/max(len(recs),1):.1f}%  '
      f'(主方向 {"下调" if main_sign < 0 else "上调"})')
print(f'原始 Δz 均值 {v_raw.mean():+.3f}   中位数 {np.median(v_raw):+.3f}')
if v_qn is not None:
    nn_qn = int((v_qn < 0).sum())
    print(f'QN 后 Δz 均值 {v_qn.mean():+.3f}  同向 {max(nn_qn, len(v_qn)-nn_qn)}/{len(v_qn)}')

# 均值应当与 geneset_shift_diagnostic 的 mean_dz 对上,便于交叉验算
print(f'\n[验算] 本脚本 mean(dz_raw) 应 ≈ shift_diag 的 mean_dz: {v_raw.mean():+.4f}')

if args.out:
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with io.open(args.out, 'w', encoding='utf-8') as f:
        f.write('gene\tquery\tdz_raw\tdz_qn\tsame_as_main\n')
        for r in recs:
            same = (r['dz_raw'] < 0) == (main_sign < 0)
            f.write(f"{r['gene']}\t{r['query']}\t{r['dz_raw']:.6f}\t{r['dz_qn']:.6f}\t{same}\n")
    print(f'-> {args.out}')

    js = os.path.splitext(args.out)[0] + '.json'
    with io.open(js, 'w', encoding='utf-8') as f:
        json.dump({
            'tag': args.tag, 'set': args.setname,
            'n_query': len(want), 'n_present': len(recs), 'absent': [g for g, _ in absent],
            'strata': [m[2] for m in masks],
            'main_sign': int(main_sign),
            'n_same_sign': int(max(n_neg, len(recs) - n_neg)),
            'mean_dz_raw': float(v_raw.mean()),
            'mean_dz_qn': float(v_qn.mean()) if v_qn is not None else None,
            'per_gene': recs,
        }, f, ensure_ascii=False, indent=2)
    print(f'-> {js}')
