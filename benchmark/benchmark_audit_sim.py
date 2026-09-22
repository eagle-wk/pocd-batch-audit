#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
模拟数据基准:量化四条防线的敏感度与特异度,并与 ComBat 类校正方法做 head-to-head。

为什么需要这个
--------------
真实公共数据没有 ground truth,所以"我们的诊断对不对"无法在真实数据上验证。
方法学论文必须有这一步:在已知真值的模拟数据上回答
  (1) 每条防线的假阳性率与检出率各是多少;
  (2) 现有校正工具(ComBat 类)在同样的设计下会不会失效。

五个情景(每个都是公共数据里的真实常见情形)
------------------------------------------
  A no_batch       无批次标签、无位移          期望建议 proceed
  B batch_crossed  真批次,与分组交叉           期望建议 correct
  C batch_confound 批次与分组完全共线           期望建议 stop
  D pseudo_batch   伪批次标签(随机划分,无效应)  期望建议 proceed(无可分证据)
  E global_shift   无批次,但处理组整体位移      期望建议 normalize

判定逻辑(与 batch_audit_generic.py 一致)
---------------------------------------
  1. 全局位移  :全基因组 Δz 中位数 / 负值比例是否偏离零分布
  2. 共线      :批次分区与分组分区是否恒等
  3. 可区分    :用技术协变量留一分类 + 置换检验(不是只看准确率阈值)
  4. 决策      :位移 > 共线 > 可区分 > 默认

用法
----
  python benchmark_audit_sim.py --reps 100 --perm 200 --outdir bench_sim
"""

import argparse
import json
import os
import sys

import numpy as np

SCEN = ['A_no_batch', 'B_batch_crossed', 'C_batch_confound',
        'D_pseudo_batch', 'E_global_shift']
TRUTH = {'A_no_batch': 'proceed',
         'B_batch_crossed': 'correct',
         'C_batch_confound': 'stop',
         'D_pseudo_batch': 'proceed',
         'E_global_shift': 'normalize'}
CN = {'proceed': '直接分析', 'correct': '校正后分析',
      'stop': '不可用', 'normalize': '归一化后分析'}


# ---------------------------------------------------------------- 模拟生成
def simulate(scenario, n_per=15, p=2000, n_de=200, n_bat=300,
             d=1.0, b=1.2, shift=0.35, seed=0):
    """返回表达矩阵 X(p×n)、分组 g(n)、批次 batch(n 或 None)、技术协变量 tech(5×n)。"""
    rng = np.random.default_rng(seed)
    n = 2 * n_per
    g = np.array([0] * n_per + [1] * n_per)

    X = rng.normal(0.0, 1.0, (p, n))
    de_idx = rng.choice(p, n_de, replace=False)
    X[np.ix_(de_idx, np.where(g == 1)[0])] += d

    batch = None
    if scenario == 'B_batch_crossed':
        # 随机划分 → 自然与分组交叉
        batch = rng.permutation(np.repeat([0, 1], n // 2))
        bi = rng.choice(p, n_bat, replace=False)
        X[np.ix_(bi, np.where(batch == 1)[0])] += b
    elif scenario == 'C_batch_confound':
        batch = g.copy()                      # 与分组完全共线
        bi = rng.choice(p, n_bat, replace=False)
        X[np.ix_(bi, np.where(batch == 1)[0])] += b
    elif scenario == 'D_pseudo_batch':
        batch = rng.permutation(np.repeat([0, 1], n // 2))
        # 关键:不施加任何批次效应 —— 标签只是随机划分

    # 技术协变量(lib_size / n_genes / mt% / rpl% / dup%)
    tech = rng.normal(0.0, 0.5, (5, n))
    if scenario in ('B_batch_crossed', 'C_batch_confound'):
        for k in range(5):
            tech[k] += float(rng.normal(0, 1)) * b * (batch == 1)

    # 全局表达位移:处理组整体平移(文库组成差异在 CPM 归一化后的表现)
    if scenario == 'E_global_shift':
        X[:, g == 1] -= shift
        tech += float(rng.normal(0.3, 0.1)) * shift * 3.0 * (g == 1)

    return X, g, batch, tech, de_idx


# ------------------------------------------------------------ 防线 1:位移
def zscore_rows(M):
    return (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-12)


def detect_shift(X, g, samp_thr=0.20):
    """
    全局表达位移检测。

    [!] 为什么不用「逐基因 Δz 的中位数」作主判据:
    存在大量差异表达基因时,Δz 中位数会**系统性偏离 0**,偏离方向与 DE 基因的主导
    方向一致。模拟实测:1500 个基因里 200 个 DE(13%)就足以让 Δz 中位数偏到
    +0.06~+0.07,而它的标准误只有 ~0.012 —— 也就是 5–6 个 sigma。拿它当判据会在
    没有任何位移的情况下误报(本脚本第一版就是这个 bug,无位移情景误报 80%)。

    改用**样本层面的表达分布位置**(每个样本跨所有基因的中位数):
    要移动中位数需要 >50% 的基因改变,13% 的 DE 基因几乎推不动它,因此对 DE 比例
    不敏感,只对真正的全转录组位移敏感。Δz 中位数与负值比例仍计算,但仅作报告。
    """
    sm = np.median(X, axis=0)                    # 每个样本的表达分布位置
    d_samp = float(sm[g == 1].mean() - sm[g == 0].mean())
    Z = zscore_rows(X)
    dz = Z[:, g == 1].mean(axis=1) - Z[:, g == 0].mean(axis=1)
    med = float(np.median(dz))                   # 仅报告
    frac_neg = float((dz < 0).mean())            # 仅报告
    return bool(abs(d_samp) > samp_thr), d_samp, med, frac_neg


# ------------------------------------------------------------ 防线 2:共线
def is_confounded(batch, g):
    """批次分区与分组分区是否恒等(构造性共线)。"""
    if batch is None:
        return False
    bg = frozenset(frozenset(np.where(batch == b)[0]) for b in np.unique(batch))
    gg = frozenset(frozenset(np.where(g == v)[0]) for v in np.unique(g))
    return bg == gg


# -------------------------------------------------- 防线 3:可区分(置换检验)
def loocv_nc(F, y):
    """留一最近类中心准确率。F: d×n 标准化后的特征。"""
    n = F.shape[1]
    cls = np.unique(y)
    hit = 0
    for i in range(n):
        dist = []
        for c in cls:
            m = (y == c) & (np.arange(n) != i)
            if m.sum() == 0:
                dist.append(np.inf)
                continue
            cen = F[:, m].mean(axis=1)
            dist.append(np.sum((F[:, i] - cen) ** 2))
        hit += int(cls[int(np.argmin(dist))] == y[i])
    return hit / n


def separable(tech, batch, perm=200, seed=0, margin=0.15):
    """技术协变量能否区分批次。必须同时满足:准确率高于基线+margin,且置换 p<0.05。"""
    if batch is None:
        return False, float('nan'), float('nan')
    F = tech.copy()
    F = (F - F.mean(axis=1, keepdims=True)) / (F.std(axis=1, keepdims=True) + 1e-12)
    y = batch.astype(int)
    acc = loocv_nc(F, y)
    rng = np.random.default_rng(seed)
    n = len(y)
    null = np.empty(perm)
    for b in range(perm):
        null[b] = loocv_nc(F, rng.permutation(y))
    base = float((null > acc).mean())
    p = float((np.sum(null >= acc) + 1) / (perm + 1))
    return bool(acc > base + margin and p < 0.05), acc, p


# ------------------------------------------------------------------ 决策
def audit(X, g, batch, tech, perm=200, seed=0):
    shift, d_samp, med, frac_neg = detect_shift(X, g)
    conf = is_confounded(batch, g)
    sep, acc, p = separable(tech, batch, perm=perm, seed=seed)
    # 决策优先级(修正版):
    #  1. 结构共线(批次==分组)→ 不可用:任何校正都会把分组效应整个抹掉
    #  2. 存在「可检测批次」且非共线 → 位移应归因于该批次,建议校正
    #     (全局批次效应会让样本中位数偏移,但这正说明批次可检测、应校正而非归一化)
    #  3. 无批次标签 + 位移 → 归一化
    #  4. 有批次标签但既不可分也不共线 + 位移 → 位移无技术解释,归一化
    #  5. 其余 → 直接分析
    # 注:此优先级在「局部批次」设计上与旧版结果一致(98.6% 不变),
    #     但能正确处理「全局批次」设计下 B 情景的位移归因。
    if batch is not None and conf:
        verdict = 'stop'
    elif batch is not None and sep:
        verdict = 'correct'
    elif shift and batch is None:
        verdict = 'normalize'
    elif shift and not sep:
        verdict = 'normalize'      # 位移存在但批次不可分,无法安全校正
    else:
        verdict = 'proceed'
    return {'verdict': verdict, 'shift': shift, 'shift_d_sample': d_samp,
            'shift_med': med, 'shift_frac_neg': frac_neg, 'confounded': conf,
            'separable': sep, 'acc': acc, 'perm_p': p}


# ------------------------------------------------------- ComBat 类校正(简化)
def combat_center(X, batch):
    """按批次中心化 —— ComBat 的均值部分(无协变量版本)。共线时等价于减去分组均值。"""
    if batch is None:
        return X.copy()
    Xc = X.copy()
    grand = X.mean(axis=1, keepdims=True)
    for b in np.unique(batch):
        m = batch == b
        Xc[:, m] -= (X[:, m].mean(axis=1, keepdims=True) - grand)
    return Xc


def delta(X, g):
    return X[:, g == 1].mean(axis=1) - X[:, g == 0].mean(axis=1)


# ------------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=100)
    ap.add_argument('--perm', type=int, default=200)
    ap.add_argument('--n-per', type=int, default=15)
    ap.add_argument('--p', type=int, default=1500)
    ap.add_argument('--seed', type=int, default=20260918)
    ap.add_argument('--outdir', default='bench_sim')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows = []
    for si, sc in enumerate(SCEN):
        truth = TRUTH[sc]
        n_ok = 0
        accs, ps = [], []
        ret, infl = [], []          # ComBat: 信号保留率 / 噪声膨胀
        verdicts = {}
        for r in range(args.reps):
            seed = int(rng.integers(0, 2 ** 31 - 1))
            X, g, batch, tech, de_idx = simulate(
                sc, n_per=args.n_per, p=args.p, seed=seed)
            res = audit(X, g, batch, tech, perm=args.perm, seed=seed + 1)
            verdicts[res['verdict']] = verdicts.get(res['verdict'], 0) + 1
            n_ok += int(res['verdict'] == truth)
            if batch is not None:
                accs.append(res['acc'])
                ps.append(res['perm_p'])

            # ComBat 对照
            Xc = combat_center(X, batch)
            s0 = np.abs(delta(X, g)[de_idx]).mean()
            s1 = np.abs(delta(Xc, g)[de_idx]).mean()
            mask = np.ones(X.shape[0], bool)
            mask[de_idx] = False
            n0 = np.abs(delta(X, g)[mask]).mean()
            n1 = np.abs(delta(Xc, g)[mask]).mean()
            ret.append(s1 / (s0 + 1e-12))
            infl.append(n1 / (n0 + 1e-12))

        frac_ok = n_ok / args.reps
        rows.append({'scenario': sc, 'truth': truth, 'truth_cn': CN[truth],
                     'correct_rate': frac_ok,
                     'mean_acc': float(np.mean(accs)) if accs else float('nan'),
                     'mean_perm_p': float(np.mean(ps)) if ps else float('nan'),
                     'combat_signal_retention': float(np.mean(ret)),
                     'combat_noise_inflation': float(np.mean(infl)),
                     'verdict_dist': verdicts})
        print(f'[{sc:<17}] 真值={CN[truth]:<8} 判定正确率={frac_ok:5.1%}  '
              f'ComBat 信号保留={np.mean(ret):5.1%}  噪声膨胀={np.mean(infl):5.2f}x')

    # 汇总
    tot_ok = np.mean([r['correct_rate'] for r in rows])
    print()
    print(f'总体判定正确率 {tot_ok:.1%}  ({args.reps} 次重复 × {len(SCEN)} 情景)')

    out_json = os.path.join(args.outdir, 'benchmark_sim.json')
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump({'reps': args.reps, 'perm': args.perm, 'n_per_group': args.n_per,
                   'p_genes': args.p, 'overall_correct_rate': float(tot_ok),
                   'rows': rows}, f, ensure_ascii=False, indent=2)

    # 混淆矩阵(真值 × 判定)
    vkeys = ['proceed', 'correct', 'stop', 'normalize']
    print()
    print('混淆矩阵(行=真值, 列=我们的判定)')
    print('            ' + ''.join(f'{CN[k]:>10}' for k in vkeys))
    for r in rows:
        vd = r['verdict_dist']
        line = f'{CN[r["truth"]]:<10}  '
        for k in vkeys:
            line += f'{vd.get(k, 0):>10}'
        print(line)

    print()
    print(f'结果已写入 {out_json}')


if __name__ == '__main__':
    main()
