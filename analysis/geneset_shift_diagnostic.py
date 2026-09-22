# -*- coding: utf-8 -*-
"""
基因集 Δz 的"全局位移"诊断
================================================================================
动机
----
在 GSE297195 的 24h 数据上,同时有 5-6 个 m6A / 程序性细胞死亡子模块在手术组
"协同下调"(Δz 达 -2.6)。这种"一大批功能不相关的基因集一起动"的模式有两种解释:

  ① 真实生物学:手术后多个通路被共同抑制
  ② 技术假象:样本间表达谱整体位移(文库组成、3' 偏好、rRNA 残留比例不同),
     经 CPM 归一化后表现为"除了少数高表达基因,其余全部相对下降"

区分方法(本脚本依次做三件事)
------------------------------
  [1] 全基因组零分布
      先算每个基因的组间 Δz,看整个基因组的分布。若中位数明显偏离 0、
      或负值基因占绝大多数 → 存在全局位移。

  [2] 分位数归一化后重算
      对每个样本做 quantile normalization(抹掉整体分布差异),再重算 Δz。
      真实通路效应会保留;全局位移造成的效应会基本消失。

  [3] 表达水平偏倚检查
      检查基因集的平均表达水平是否异常高/低。高表达基因集(如含 Rpl/mt-)
      对文库组成最敏感,是假阳性的高发区。

用法
----
python geneset_shift_diagnostic.py --matrix M.tsv --samples SS.tsv \
    --genesets G.tsv --group-col group --stratum-col age \
    --case surgery --ctrl sham --tag LABEL --outdir DIR
"""
import argparse, csv, io, json, os, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--matrix', required=True)
ap.add_argument('--samples', required=True)
ap.add_argument('--genesets', required=True)
ap.add_argument('--group-col', default='group')
ap.add_argument('--stratum-col', default='age')
ap.add_argument('--sample-col', default='sample')
ap.add_argument('--case', required=True)
ap.add_argument('--ctrl', required=True)
ap.add_argument('--min-geneset', type=int, default=5)
ap.add_argument('--perm', type=int, default=5000)
ap.add_argument('--seed', type=int, default=42)
ap.add_argument('--tag', default='')
ap.add_argument('--outdir', default='shift_diag')
args = ap.parse_args()

rng = np.random.default_rng(args.seed)
os.makedirs(args.outdir, exist_ok=True)

# ---------------- 数据 ----------------
srows = list(csv.DictReader(io.open(args.samples, encoding='utf-8'), delimiter='\t'))
samples = [r[args.sample_col].strip() for r in srows]
grp = [r[args.group_col].strip() for r in srows]
stg = [r[args.stratum_col].strip() for r in srows]

rows = list(csv.reader(io.open(args.matrix, encoding='utf-8'), delimiter='\t'))
hdr = rows[0][1:]
genes = [r[0] for r in rows[1:]]
X = np.array([[float(x) for x in r[1:]] for r in rows[1:]], dtype=float)
col = {s: j for j, s in enumerate(hdr)}
order = [col[s] for s in samples]
X = X[:, order]
print(f'矩阵: {X.shape[0]} 基因 x {X.shape[1]} 样本')

gi = {g: i for i, g in enumerate(genes)}
upper = {g.upper(): g for g in genes}

sets = {}
with io.open(args.genesets, encoding='utf-8') as f:
    rd = csv.reader(f, delimiter='\t')
    next(rd, None)
    for r in rd:
        if len(r) >= 2 and r[0].strip() and r[1].strip():
            sets.setdefault(r[1].strip(), []).append(r[0].strip())

sets_idx = {}
for name, gs in sets.items():
    idx = sorted({gi[g] for g in gs if g in gi} | {gi[upper[g.upper()]] for g in gs if g.upper() in upper})
    if len(idx) >= args.min_geneset:
        sets_idx[name] = idx
print(f'基因集: {len(sets_idx)} 个')

# ---------------- 分层取样 ----------------
strata = sorted(set(stg))
masks = []                                   # (ctrl_idx, case_idx)
for s in strata:
    ci = [i for i in range(len(samples)) if stg[i] == s and grp[i] == args.ctrl]
    ki = [i for i in range(len(samples)) if stg[i] == s and grp[i] == args.case]
    if len(ci) >= 2 and len(ki) >= 2:
        masks.append((ci, ki, s, len(ci), len(ki)))
print(f'层: {[(m[2], f"{m[3]}v{m[4]}") for m in masks]}')


def zscore_rows(M):
    return (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-12)


def per_gene_dz(M):
    """每个基因的分层内组间 Δz(各层等权平均)"""
    Z = zscore_rows(M)
    acc, wsum = np.zeros(Z.shape[0]), 0.0
    for ci, ki, s, nc, nk in masks:
        acc += Z[:, ki].mean(axis=1) - Z[:, ci].mean(axis=1)
        wsum += 1
    return acc / max(wsum, 1), Z


def quantile_normalize(M):
    """标准分位数归一化: 把每个样本的分布对齐到平均分布, 抹平全局位移"""
    R = M.shape[0]
    ranks = np.argsort(np.argsort(M, axis=0), axis=0)      # 0..R-1
    sorted_M = np.sort(M, axis=0)
    ref = sorted_M.mean(axis=1)                            # 平均分位分布
    out = np.empty_like(M)
    for j in range(M.shape[1]):
        out[:, j] = ref[ranks[:, j]]
    return out


def report(M, label):
    dz, Z = per_gene_dz(M)
    neg = float((dz < 0).mean())
    print(f'\n{"="*78}')
    print(f'【{label}】全基因组 Δz 分布')
    print(f'  中位数 {np.median(dz):+.4f} | 均值 {dz.mean():+.4f} | '
          f'IQR [{np.percentile(dz,25):+.3f}, {np.percentile(dz,75):+.3f}]')
    print(f'  Δz<0 的基因占 {100*neg:.1f}%   (无全局位移应接近 50%)')

    gz = np.abs(dz)
    res = []
    for name, idx in sets_idx.items():
        obs = float(dz[idx].mean())
        B = args.perm
        null = np.abs(dz[rng.integers(0, dz.size, (B, len(idx)))].mean(axis=1))
        p = float((null >= abs(obs)).mean())
        res.append({
            'set': name, 'n_genes': len(idx), 'mean_dz': obs,
            'emp_p': p, 'mean_expr': float(M[idx].mean()),
            'median_expr_rank_pct': float((np.argsort(np.argsort(M.mean(axis=1)))[idx] / M.shape[0]).mean() * 100),
        })
    res.sort(key=lambda r: r['emp_p'])
    print(f'\n  {"基因集":<34}{"n":>4}{"平均Δz":>10}{"经验p":>9}{"表达秩%":>9}')
    for r in res[:14]:
        print(f"  {r['set']:<34}{r['n_genes']:>4}{r['mean_dz']:>+10.3f}"
              f"{r['emp_p']:>9.4f}{r['median_expr_rank_pct']:>9.1f}")
    return dz, res


print('\n' + '#' * 78)
print(f'# 原始矩阵 (CPM-log2)   tag={args.tag}')
print('#' * 78)
dz_raw, res_raw = report(X, '原始')

print('\n' + '#' * 78)
print(f'# 分位数归一化后 (全局位移被抹平)   tag={args.tag}')
print('#' * 78)
Xq = quantile_normalize(X)
dz_qn, res_qn = report(Xq, '分位数归一化')

# ---------------- 汇总 ----------------
print('\n' + '=' * 78)
print('汇总: 信号是否依赖全局位移?')
print('=' * 78)
qm = {r['set']: r for r in res_qn}
print(f'{"基因集":<34}{"原始Δz":>10}{"原始p":>10}{"QN后Δz":>10}{"QN后p":>10}{"结论":>12}')
out_rows = []
for r in res_raw:
    q = qm.get(r['set'])
    if not q:
        continue
    if r['emp_p'] < 0.05 and q['emp_p'] < 0.05:
        verdict = '[稳健-特异]'
    elif r['emp_p'] < 0.05 and q['emp_p'] >= 0.05:
        verdict = '[全局位移]'
    elif r['emp_p'] >= 0.05 and q['emp_p'] < 0.05:
        verdict = '[QN后显现]'
    else:
        verdict = '[无信号]'
    out_rows.append({**r, 'qn_mean_dz': q['mean_dz'], 'qn_emp_p': q['emp_p'], 'verdict': verdict})
    print(f"{r['set']:<34}{r['mean_dz']:>+10.3f}{r['emp_p']:>10.4f}"
          f"{q['mean_dz']:>+10.3f}{q['emp_p']:>10.4f}{verdict:>12}")

n_spec = sum(1 for r in out_rows if r['verdict'] == '[稳健-特异]')
n_glob = sum(1 for r in out_rows if r['verdict'] == '[全局位移]')
print(f'\n稳健-特异: {n_spec} | 疑似全局位移: {n_glob} | 总计 {len(out_rows)}')

with io.open(os.path.join(args.outdir, f'shift_diag{("_" + args.tag) if args.tag else ""}.json'),
             'w', encoding='utf-8') as f:
    json.dump({
        'tag': args.tag, 'strata': [m[2] for m in masks],
        'genome_median_dz_raw': float(np.median(dz_raw)),
        'genome_frac_negative_raw': float((dz_raw < 0).mean()),
        'genome_median_dz_qn': float(np.median(dz_qn)),
        'genome_frac_negative_qn': float((dz_qn < 0).mean()),
        'results': out_rows,
    }, f, ensure_ascii=False, indent=2)
print(f'\n-> {args.outdir}')
