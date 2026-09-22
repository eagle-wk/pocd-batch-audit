# -*- coding: utf-8 -*-
"""
跨数据集 meta 分析:把 GSE330865 与 GSE297195 的集合级证据合并。

为什么需要:
  单看每个数据集,`IFN_signaling` 在 GSE330865 是 p=0.037、在 GSE297195 是 p=0.0076 ——
  都是"刚过线"。但这两个是**互相独立**的研究(不同实验室、不同动物、不同建库、不同批次结构),
  把两个独立证据合并,才能给出一个有量的结论,而不是停留在"两边方向一样"。

方法(两种,互为佐证,都不需要原始数据):
  Stouffer 加权:  Z = Σ(wᵢ·Zᵢ) / √(Σwᵢ²),  wᵢ = √nᵢ
  Fisher:         χ² = −2·Σ ln(pᵢ),  df = 2k

⚠️ 前提与限制(必须写进报告):
  - 两侧 p 均取自"随机等大小基因集"零分布(特异化检验),对相关基因偏乐观;
  - 两个数据集存在同一类设计缺陷(sham/surgery 各自 GSM 连续成块),
    因此这不是两个完全独立的证据 —— 只能算"两个独立队列 + 独立建库"的部分独立复现;
  - 多重检验:对全部候选集合做 BH 校正。
"""
import csv, io, json, math, os
import numpy as np
from scipy.stats import norm, chi2

HERE = os.path.dirname(os.path.abspath(__file__))        # <repo>/analysis
REPO = os.path.dirname(HERE)                            # <repo>
SD = os.path.join(REPO, 'results', 'derived', 'shift_diag')   # stage-2 inputs
OUT = os.path.join(REPO, 'results', 'derived')          # stage-4 outputs
# The two datasets' shift-diagnostic records are flattened into one directory
# in the repository; the filenames already carry the dataset prefix.
G865, G195 = SD, SD

N = {'GSE330865': 30, 'GSE297195(24h)': 23}     # 用于加权


def load(fp):
    if not os.path.exists(fp):
        return {}
    j = json.load(io.open(fp, encoding='utf-8'))
    return {r['set']: r for r in j['results']}


SRC = {
    'GSE330865': {
        'm6a_pcd_submodules': load(os.path.join(G865, 'shift_diag_GSE330865_m6aPCD.json')),
        'modules': load(os.path.join(G865, 'shift_diag_GSE330865_modules.json')),
        'interferon': load(os.path.join(G865, 'shift_diag_IFN_GSE330865.json')),
    },
    'GSE297195(24h)': {
        'm6a_pcd_submodules': load(os.path.join(G195, 'shift_diag_GSE297195_24h.json')),
        'modules': load(os.path.join(G195, 'shift_diag_GSE297195_24h_GSE330865modules.json')),
        'interferon': load(os.path.join(G195, 'shift_diag_IFN_GSE297195_24h.json')),
    },
}
FAMILY = {'m6a_pcd_submodules': 'm6A/PCD', 'modules': 'WGCNA模块', 'interferon': '干扰素(预定义)'}

rows = []
for fam, d865 in SRC['GSE330865'].items():
    d195 = SRC['GSE297195(24h)'].get(fam, {})
    for name in sorted(set(d865) | set(d195)):
        a, b = d865.get(name), d195.get(name)
        if not a or not b:
            continue
        # 用 QN 校正后的特异化 p(双侧) → 单侧(带方向)
        pa, pb = a['qn_emp_p'], b['qn_emp_p']
        da, db = a['qn_mean_dz'], b['qn_mean_dz']
        pa1 = max(min(pa / 2 if da < 0 else 1 - pa / 2, 1), 1e-12)
        pb1 = max(min(pb / 2 if db < 0 else 1 - pb / 2, 1), 1e-12)
        za, zb = norm.isf(pa1), norm.isf(pb1)
        wa, wb = math.sqrt(N['GSE330865']), math.sqrt(N['GSE297195(24h)'])
        z_st = (wa * za + wb * zb) / math.sqrt(wa ** 2 + wb ** 2)
        p_st = norm.sf(z_st)
        chi = -2 * (math.log(pa1) + math.log(pb1))
        p_fi = chi2.sf(chi, 4)
        rows.append({
            'family': FAMILY[fam], 'set': name, 'n_genes': a['n_genes'],
            'dz_865': da, 'p_865': pa, 'dz_195': db, 'p_195': pb,
            'same_sign': (da > 0) == (db > 0) and abs(da) > 1e-9 and abs(db) > 1e-9,
            'z_stouffer': z_st, 'p_stouffer_1s': p_st, 'p_fisher': p_fi,
        })

# BH 校正(对 Stouffer 单侧 p)
rows.sort(key=lambda r: r['p_stouffer_1s'])
m = len(rows)
prev = 1.0
for i in range(m - 1, -1, -1):
    prev = min(prev, rows[i]['p_stouffer_1s'] * m / (i + 1))
    rows[i]['p_bh'] = prev

print('=' * 118)
print('跨数据集 meta 分析 (GSE330865 n=30  ⊕  GSE297195 24h n=23)')
print('  合并的是"集合相对随机等大小基因集"的特异化检验(已做分位数归一化);单侧,方向由 Δz 决定')
print('=' * 118)
print(f"{'家族':<14}{'集合':<32}{'n':>4}{'865Δz':>9}{'865p':>7}"
      f"{'195Δz':>9}{'195p':>7}{'Stouffer Z':>11}{'合并p':>9}{'BH p':>8}{'同向':>5}")
for r in rows:
    print(f"{r['family']:<14}{r['set']:<32}{r['n_genes']:>4}{r['dz_865']:>+9.3f}{r['p_865']:>7.3f}"
          f"{r['dz_195']:>+9.3f}{r['p_195']:>7.3f}{r['z_stouffer']:>+11.2f}"
          f"{r['p_stouffer_1s']:>9.4f}{r['p_bh']:>8.3f}{'  是' if r['same_sign'] else '  否':>5}")

sig = [r for r in rows if r['p_bh'] < 0.05]
print(f'\nBH<0.05 的集合: {len(sig)} 个')
for r in sig:
    print(f"  ★ {r['set']:<32} n={r['n_genes']:<3} "
          f"GSE330865 {r['dz_865']:+.3f}(p={r['p_865']:.3f})  "
          f"GSE297195 {r['dz_195']:+.3f}(p={r['p_195']:.3f})  "
          f"Stouffer Z={r['z_stouffer']:+.2f}  合并p={r['p_stouffer_1s']:.5f}  BH={r['p_bh']:.4f}")

print('\n--- 方向一致且两侧都名义显著(QN 后 p<0.05)的集合 ---')
both = [r for r in rows if r['same_sign'] and r['p_865'] < 0.05 and r['p_195'] < 0.05]
if both:
    for r in both:
        print(f"  {r['family']:<14}{r['set']:<32} {r['dz_865']:+.3f} / {r['dz_195']:+.3f}")
else:
    print('  (无)')
print(f'\n全部候选中两侧都名义显著者: {len(both)}/{len(rows)}'
      f'  (若各集合独立、单侧命中率约 2.5%,随机期望 {0.025**2*len(rows):.3f} 个)')

json.dump(rows, io.open(os.path.join(OUT, 'cross_dataset_meta.json'), 'w', encoding='utf-8'),
          ensure_ascii=False, indent=2)
with io.open(os.path.join(OUT, 'cross_dataset_meta.tsv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader()
    for r in rows:
        w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
print(f'\n-> {os.path.join(OUT, "cross_dataset_meta.tsv")}')
