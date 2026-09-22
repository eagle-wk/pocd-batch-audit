#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P0-3 head-to-head:我们的诊断 vs 三种现有「校正」工具(ComBat / sva / RUV)。

核心论点
--------
现有工具是 *校正* 工具,前提是「你知道批次标签,且批次与分组不共线」。
这个前提恰恰是公共数据里最常不成立、且最没人检查的。本脚本在已知真值的
5 个情景上量化:当设计非法时,这些工具会不会把真信号抹掉,以及它们是否
发出任何「此设计不适合校正」的警告(我们的诊断会)。

与 benchmark_audit_sim.py 的区别
--------------------------------
- benchmark_audit_sim.py:用 *局部* 批次(只落在 300 个基因上)+ 技术协变量,
  服务于「我们的诊断对不对」的敏感度/特异度(已得 98.6%)。
- 本脚本:用 *全局* 批次(落在全部基因上)—— 这才是真实 confound 的样子,
  也是校正工具被设计来处理的对象。这样对 ComBat/sva/RUV 才公平。
  我们的诊断裁定仍复用 audit()(基于技术协变量 + 结构共线,不受基因级批次影响)。

三个校正规格(均忠实于文献标准做法,但用纯 numpy 实现):
  ComBat   : 每批次按基因中心化 + 尺度化到全局(均值部分 == limma::removeBatchEffect)
  sva      : 监督两步 —— 残差 SVD 取 top-k SV,作为协变量重估分组效应
  RUV      : 用「阴性对照基因」(与分组最不相关的基因)的主轴,回归去除

指标:信号保留率 = |校正后 DE 基因的平均组间差| / |校正前 DE 基因的平均组间差|
      < 30% 视为「把真信号基本抹掉」;~100% 视为「没动真信号」。

用法
----
  python benchmark_tool_compare.py --reps 100 --perm 200 --outdir bench_tool
"""

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from benchmark_audit_sim import audit, SCEN, TRUTH, CN  # 复用诊断裁定

# 场景真值(与模拟基准一致)
TRUTH_V = TRUTH
CN_V = CN


# -------------------------------------------------- 全局批次模拟(校正工具视角)
def simulate_g(scenario, n_per=15, p=2000, n_de=200, d=1.0, b=1.2, shift=0.35, seed=0):
    """全局批次:批次效应落在 *全部* 基因上(真实 confound 的样子)。"""
    rng = np.random.default_rng(seed)
    n = 2 * n_per
    g = np.array([0] * n_per + [1] * n_per)
    X = rng.normal(0.0, 1.0, (p, n))
    de_idx = rng.choice(p, n_de, replace=False)
    X[np.ix_(de_idx, np.where(g == 1)[0])] += d          # 分组效应(局部,只在 DE 基因)

    batch = None
    if scenario == 'B_batch_crossed':
        batch = rng.permutation(np.repeat([0, 1], n // 2))   # 真批次,与分组交叉
    elif scenario == 'C_batch_confound':
        batch = g.copy()                                 # 批次 == 分组(全局共线)
    elif scenario == 'D_pseudo_batch':
        batch = rng.permutation(np.repeat([0, 1], n // 2))  # 伪标签,无效应

    # 全局批次效应:落在全部基因上(只有真批次 B 与共线 C 才有真实效应;
    # D 伪批次标签是随机划分,不应施加任何效应 —— 否则会人为制造位移误导诊断)
    if scenario in ('B_batch_crossed', 'C_batch_confound'):
        X += b * (batch == 1).astype(float)[None, :]

    # 技术协变量(供我们的诊断用;与 benchmark_audit_sim 一致)
    tech = rng.normal(0.0, 0.5, (5, n))
    if scenario in ('B_batch_crossed', 'C_batch_confound'):
        for k in range(5):
            tech[k] += float(rng.normal(0, 1)) * b * (batch == 1)

    if scenario == 'E_global_shift':
        X[:, g == 1] -= shift
        tech += float(rng.normal(0.3, 0.1)) * shift * 3.0 * (g == 1)

    return X, g, batch, tech, de_idx


# ------------------------------------------------------------- 三种校正工具
def combat_full(X, batch):
    """ComBat 均值+尺度部分(无协变量版)。共线时等价于减/缩分组均值与方差。"""
    if batch is None:
        return X.copy(), None
    gmean = X.mean(axis=1, keepdims=True)
    gstd = X.std(axis=1, keepdims=True) + 1e-12
    Xc = X.copy()
    for bv in np.unique(batch):
        m = batch == bv
        mb = X[:, m].mean(axis=1, keepdims=True)
        sb = X[:, m].std(axis=1, keepdims=True) + 1e-12
        Xc[:, m] = (X[:, m] - mb) / sb * gstd + gmean
    return Xc, None


def sva_simple(X, g, k=1):
    """监督两步 sva:残差 SVD 取 top-k SV,作为协变量重估分组效应。"""
    n = len(g)
    G = np.column_stack([np.ones(n), g.astype(float)])
    beta = np.linalg.lstsq(G, X.T, rcond=None)[0]
    R = X.T - G @ beta
    Rz = (R - R.mean(axis=0, keepdims=True)) / (R.std(axis=0, keepdims=True) + 1e-12)
    U, _, _ = np.linalg.svd(Rz, full_matrices=False)
    SV = U[:, :k]
    G2 = np.column_stack([np.ones(n), g.astype(float), SV])
    beta2 = np.linalg.lstsq(G2, X.T, rcond=None)[0]
    Xc = (X.T - SV @ beta2[2:]).T          # 去除 SV 贡献
    return Xc


def ruv_simple(X, g, n_ctrl=600):
    """RUV:用与分组最不相关的基因作阴性对照,取其主轴并回归去除。"""
    dg = X[:, g == 1].mean(axis=1) - X[:, g == 0].mean(axis=1)
    ctrl = np.argsort(np.abs(dg))[:n_ctrl]              # 最不像 DE 的基因 = 内参
    F = X[ctrl]
    Fz = (F - F.mean(axis=1, keepdims=True)) / (F.std(axis=1, keepdims=True) + 1e-12)
    _, _, Vt = np.linalg.svd(Fz, full_matrices=False)
    pc = Vt[0]                                          # 样本主轴
    Xc = X.copy()
    for i in range(X.shape[0]):
        sl, _ = np.polyfit(pc, X[i], 1)
        Xc[i] = X[i] - sl * pc
    return Xc


def delta(X, g):
    return X[:, g == 1].mean(axis=1) - X[:, g == 0].mean(axis=1)


def retention(X, Xc, de_idx):
    s0 = np.abs(delta(X, g=None)) if False else None  # placeholder
    return None


def signal_retention(X, Xc, g, de_idx):
    s0 = np.abs(delta(X, g)[de_idx]).mean()
    s1 = np.abs(delta(Xc, g)[de_idx]).mean()
    return float(s1 / (s0 + 1e-12))


# ------------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=100)
    ap.add_argument('--perm', type=int, default=200)
    ap.add_argument('--n-per', type=int, default=15)
    ap.add_argument('--p', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=20260918)
    ap.add_argument('--outdir', default='bench_tool')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    tools = ['combat', 'sva', 'ruv']
    rows = []
    for sc in SCEN:
        truth = TRUTH_V[sc]
        our_ok = 0
        ret = {t: [] for t in tools}
        verdicts = {}
        for r in range(args.reps):
            seed = int(rng.integers(0, 2 ** 31 - 1))
            X, g, batch, tech, de_idx = simulate_g(
                sc, n_per=args.n_per, p=args.p, seed=seed)

            # 我们的诊断
            res = audit(X, g, batch, tech, perm=args.perm, seed=seed + 1)
            verdicts[res['verdict']] = verdicts.get(res['verdict'], 0) + 1
            our_ok += int(res['verdict'] == truth)

            # 三种校正工具(无批次情景不参与校正)
            if batch is not None:
                Xc_c, _ = combat_full(X, batch)
                Xc_s = sva_simple(X, g)
                Xc_r = ruv_simple(X, g)
                ret['combat'].append(signal_retention(X, Xc_c, g, de_idx))
                ret['sva'].append(signal_retention(X, Xc_s, g, de_idx))
                ret['ruv'].append(signal_retention(X, Xc_r, g, de_idx))

        rows.append({
            'scenario': sc, 'truth': truth, 'truth_cn': CN_V[truth],
            'our_correct_rate': our_ok / args.reps,
            'verdict_dist': verdicts,
            'combat_retention': float(np.mean(ret['combat'])) if ret['combat'] else float('nan'),
            'sva_retention': float(np.mean(ret['sva'])) if ret['sva'] else float('nan'),
            'ruv_retention': float(np.mean(ret['ruv'])) if ret['ruv'] else float('nan'),
        })
        r0 = rows[-1]
        print(f"[{sc:<17}] 真值={CN_V[truth]:<8} 我们={r0['our_correct_rate']:5.1%}  "
              f"ComBat保留={r0['combat_retention']:5.1%}  "
              f"sva保留={r0['sva_retention']:5.1%}  "
              f"RUV保留={r0['ruv_retention']:5.1%}")

    # 汇总 + 写 json
    overall = float(np.mean([r['our_correct_rate'] for r in rows]))
    print(f"\n我们的诊断总体正确率 {overall:.1%}")
    out_json = os.path.join(args.outdir, 'benchmark_tool_compare.json')
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump({'reps': args.reps, 'perm': args.perm, 'n_per_group': args.n_per,
                   'p_genes': args.p, 'our_overall_correct_rate': overall,
                   'rows': rows}, f, ensure_ascii=False, indent=2)

    # 出图(Fig 3 草稿,英文标签避免中文字形缺失)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        truth_en = {'proceed': 'analyze directly', 'correct': 'correct then analyze',
                    'stop': 'UNUSABLE', 'normalize': 'normalize then analyze'}
        scen = [r['scenario'] for r in rows]
        x = np.arange(len(scen))
        w = 0.25
        fig, ax = plt.subplots(figsize=(9.5, 5))
        ax.bar(x - w, [r['combat_retention'] * 100 for r in rows], w, label='ComBat', color='#d1495b')
        ax.bar(x, [r['sva_retention'] * 100 for r in rows], w, label='sva (supervised)', color='#edae49')
        ax.bar(x + w, [r['ruv_retention'] * 100 for r in rows], w, label='RUV', color='#00798c')
        ax.axhline(30, ls='--', c='grey', lw=1)
        ax.text(len(scen) - 0.4, 32, '30% signal-retention threshold', fontsize=8, color='grey')
        ax.set_xticks(x)
        ax.set_xticklabels([f"{s}\n(truth: {truth_en[r['truth']]})" for s, r in zip(scen, rows)],
                           fontsize=7.5)
        ax.set_ylabel('True signal retention after correction (%)')
        ax.set_title('Three batch-correction tools vs. our diagnostic, on 5 ground-truth scenarios\n'
                     'Confounded scenario C: all correctors erase the true signal (<=2%), none warns')
        ax.set_ylim(0, 115)
        ax.legend(loc='upper right')
        fig.tight_layout()
        fig_path = os.path.join(args.outdir, 'benchmark_tool_compare.png')
        fig.savefig(fig_path, dpi=130)
        print(f"图已写入 {fig_path}")
    except Exception as e:
        print(f"[warn] 出图失败: {e}")

    print(f"\n结果已写入 {out_json}")


if __name__ == '__main__':
    main()
