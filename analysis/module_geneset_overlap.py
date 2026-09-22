# -*- coding: utf-8 -*-
"""
模块 × 预定义基因集 重叠检验(防"看 hub 基因猜身份"的过度推断)
================================================================================
动机(一次真实的错误)
--------------------
GSE330865 的 royalblue 模块,前 5 个 hub 是 Trim21 / Herc6 / Zc3hav1(ZAP) /
H2-Q7 / Mr1 —— 四个典型的 I 型干扰素刺激基因加一个 MHC-Ib。据此我在报告里
把它写成"干扰素 / MHC-I 抗原提呈模块",并由此推出"术后 ISG 上调"的假说。

后来用**文献驱动、与本研究数据无关**的预定义干扰素集合去核,发现:

    royalblue 的 121 个基因里,只有 5 个落入预定义干扰素集合 = 4.1%

而真正的免疫模块是 darkturquoise(13/34 = 38.2%)。也就是说,仅凭 hub 基因的
生物学标签去命名模块,在 hub 恰好是"知名基因"时会系统性产生误判。

教训:模块的身份必须用**外部的、预先定义的集合**做富集检验(超几何/Fisher),
而不是读基因名。本脚本把这一步固化下来。

判读
----
  fold > 1 且 p 小  → 该模块确实是这个集合对应的生物学模块
  fold ≈ 1、p 大    → 无关(即使 hub 名字看着像)
  fold < 1          → 该模块在此集合上反而"贫化",提示身份是别的东西

用法
----
python module_geneset_overlap.py \
    --modules  wgcna_results/module_genesets.tsv \
    --genesets ../platform_annot/interferon_mhc_genesets_mouse.tsv \
    --background wgcna_results/gene_module_assignment.csv \
    --out  tsv/模块_基因集重叠.tsv
"""
import argparse, csv, io, os, sys
from math import comb, lgamma, exp, log

ap = argparse.ArgumentParser()
ap.add_argument('--modules', required=True, help='TSV: gene<TAB>set (模块名在 set 列)')
ap.add_argument('--genesets', required=True, help='TSV: gene<TAB>set,可带 # 注释行')
ap.add_argument('--background', default='', help='可选: 参与检验的背景基因列表文件(取第一列)')
ap.add_argument('--fdr', type=float, default=0.05)
ap.add_argument('--out', default='')
args = ap.parse_args()


def read_sets(path):
    """读 gene<TAB>set 长表;跳过 # 注释行与 'gene<TAB>set' 表头"""
    d = {}
    with io.open(path, encoding='utf-8-sig') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) >= 2 and p[0].strip() and p[1].strip():
                if p[0].strip().lower() == 'gene' and p[1].strip().lower() == 'set':
                    continue                      # 表头
                d.setdefault(p[1].strip(), set()).add(p[0].strip())
    return d


modules = read_sets(args.modules)
genesets = read_sets(args.genesets)

if args.background:
    bg = set()
    with io.open(args.background, encoding='utf-8') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            if i == 0 and (line.lower().startswith('gene') or ',' in line and 'module' in line.lower()):
                continue
            bg.add(line.split('\t')[0].split(',')[0].strip().strip('"'))
else:
    bg = set().union(*[v for v in modules.values()], *[v for v in genesets.values()])
N = len(bg)
print(f'模块 {len(modules)} 个 | 基因集 {len(genesets)} 个 | 背景基因 {N} 个')

# —— 关键前置检查:集合成员是否真的在背景里 ——
# 若某个集合的基因几乎不在 WGCNA 网络内,那么"该模块不是这个集合的模块"这一结论
# 就不是生物学判断,而是"根本没纳入分析"。这种情况必须单独指出,否则会被误读。
print('\n[前置检查] 各集合有多少成员落在背景(WGCNA 已分配基因)内:')
untestable = []
for sn, sg in sorted(genesets.items(), key=lambda kv: -len(kv[1])):
    tin = len(sg & bg)
    print(f'  {sn:<18} {tin:>3}/{len(sg):<3} 可检验'
          + ('   ← 全部不在网络内,无法判断' if tin == 0 else ''))
    if tin == 0 and len(sg) >= 3:
        untestable.append(sn)
if untestable:
    print(f'  [!] {", ".join(untestable)} 完全不在 WGCNA 网络内 —— '
          f'对它们只能做基因集层面评分,不能用模块身份去推断(反之亦然)')



def _logcomb(n, k):
    """log C(n, k);n 可达 2000+,必须走 lgamma 否则 comb 溢出"""
    if k < 0 or k > n or n < 0:
        return float('-inf')
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def hypergeom_tail(k, n, K, N):
    """P(X >= k),X ~ Hypergeom(N, K, n):从 N 中抽 n 个,其中属于大小为 K 的类的个数 >= k

    用 log 空间做 log-sum-exp,避免 C(1913, 456) 这类天文数字溢出。
    """
    if N < n or K > N or n > N:
        return float('nan')
    hi = min(n, K)
    terms = []
    for x in range(k, hi + 1):
        if n - x < 0 or N - K < n - x:
            continue
        lg = _logcomb(K, x) + _logcomb(N - K, n - x) - _logcomb(N, n)
        terms.append(lg)
    if not terms:
        return 0.0
    mx = max(terms)
    if mx == float('-inf'):
        return 0.0
    return exp(mx) * sum(exp(t - mx) for t in terms) if mx > -700 else 0.0



rows = []
for mn, mg in sorted(modules.items(), key=lambda kv: -len(kv[1])):
    mg = mg & bg
    for sn, sg in sorted(genesets.items()):
        sg_bg = sg & bg
        K, n = len(sg_bg), len(mg)
        if K == 0 or n == 0:
            continue
        k = len(mg & sg_bg)
        expect = n * K / N if N else float('nan')      # 注意: 不能叫 exp,会遮蔽 math.exp
        fold = k / expect if expect else float('nan')
        p = hypergeom_tail(k, n, K, N)
        rows.append({'module': mn, 'geneset': sn, 'n_module': n, 'n_set': K,
                     'n_overlap': k, 'expected': expect, 'fold': fold, 'p': p,
                     'overlap_pct_of_module': 100 * k / n,
                     'members': sorted(mg & sg_bg)})

# BH
rows.sort(key=lambda r: r['p'])
m = len(rows)
prev = 1.0
for i in range(m - 1, -1, -1):
    prev = min(prev, rows[i]['p'] * m / (i + 1))
    rows[i]['padj'] = min(prev, 1.0)
rows.sort(key=lambda r: -r['fold'])

print(f'\n{"模块":<18}{"基因集":<18}{"n模块":>6}{"n集合":>6}{"重叠":>6}'
      f'{"期望":>8}{"fold":>7}{"p":>11}{"padj":>9}{"占模块%":>9}')
for r in rows:
    if r['n_overlap'] == 0 and r['fold'] < 1.5:
        continue
    print(f"{r['module']:<18}{r['geneset']:<18}{r['n_module']:>6}{r['n_set']:>6}"
          f"{r['n_overlap']:>6}{r['expected']:>8.2f}{r['fold']:>7.2f}"
          f"{r['p']:>11.3g}{r['padj']:>9.3g}{r['overlap_pct_of_module']:>9.1f}")

print('\n判定(FDR<%.2f):' % args.fdr)
hit = [r for r in rows if r['padj'] < args.fdr and r['fold'] > 1]
if hit:
    for r in hit:
        print(f"  [确证] {r['module']} ↔ {r['geneset']}: "
              f"{r['n_overlap']}/{r['n_module']} ({r['overlap_pct_of_module']:.1f}%), "
              f"fold={r['fold']:.1f}, padj={r['padj']:.3g}")
else:
    print('  (无)')

# 专门打印"名义显著但未过 FDR"的,以及单集合被多个模块命中的情况,便于判断
print('\n参考:各模块的最高 fold 命中')
best = {}
for r in rows:
    if r['n_overlap'] > 0:
        best.setdefault(r['module'], r)
for mn, r in sorted(best.items(), key=lambda kv: -kv[1]['fold']):
    print(f"  {mn:<18} → {r['geneset']:<18} {r['n_overlap']}/{r['n_module']} "
          f"({r['overlap_pct_of_module']:.1f}%) fold={r['fold']:.2f} p={r['p']:.3g}")

# 家族级(所有集合取并集)重叠 —— 回答"这个模块到底是不是干扰素模块"
fam = set().union(*[v for v in genesets.values()])
fam_bg = fam & bg
print(f'\n家族级重叠(所有 {len(genesets)} 个集合取并集, 共 {len(fam)} 基因 / {len(fam_bg)} 可检验):')
print(f'  {"模块":<18}{"n模块":>6}{"重叠":>6}{"期望":>8}{"fold":>7}{"p":>11}{"占模块%":>9}')
fam_rows = []
for mn, mg in sorted(modules.items(), key=lambda kv: -len(kv[1])):
    mg = mg & bg
    if not mg or not fam_bg:
        continue
    k = len(mg & fam_bg)
    expect = len(mg) * len(fam_bg) / N
    fold = k / expect if expect else float('nan')
    p = hypergeom_tail(k, len(mg), len(fam_bg), N)
    fam_rows.append((mn, len(mg), k, expect, fold, p, 100 * k / len(mg)))
fam_rows.sort(key=lambda r: -r[4])
for mn, n, k, expect, fold, p, pct in fam_rows:
    print(f'  {mn:<18}{n:>6}{k:>6}{expect:>8.1f}{fold:>7.2f}{p:>11.3g}{pct:>9.1f}')
if fam_rows:
    top = fam_rows[0]
    print(f'  → 重叠率最高的是 {top[0]}:{top[2]}/{top[1]} = {top[6]:.1f}%')

if args.out:
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with io.open(args.out, 'w', encoding='utf-8') as f:
        f.write('module\tgeneset\tn_module\tn_set\tn_overlap\texpected\tfold\tp\tpadj\tpct_of_module\tmembers\n')
        for r in rows:
            f.write(f"{r['module']}\t{r['geneset']}\t{r['n_module']}\t{r['n_set']}\t"
                    f"{r['n_overlap']}\t{r['expected']:.3f}\t{r['fold']:.4f}\t{r['p']:.6g}\t"
                    f"{r['padj']:.6g}\t{r['overlap_pct_of_module']:.2f}\t{','.join(r['members'])}\n")
    print(f'\n-> {args.out}')
