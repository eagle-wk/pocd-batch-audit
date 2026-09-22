# -*- coding: utf-8 -*-
"""
通用批次审计 —— 在任何候选数据集被纳入分析前,先回答"它能不能用"

设计目标:不再重蹈 GSE276942 的覆辙(全部分析做完才发现批次与时间点 100% 混杂)。

输入:
  1. 计数矩阵 (gene x sample, tsv/csv, 第一列为基因 ID)
  2. 样本表 (tsv),必需列: sample, group;可选列: batch, tissue
     —— 若 batch 列缺失,脚本会用文库大小自动做双峰切分并**明确标注为推断批次**

输出(src_out/):
  technical_covariates.csv   每个样本的技术指标
  confound_crosstab.csv      批次 × 分组 交叉表 + 混杂判定
  pca_attribution.csv        PC1..PC8 的 R²(批次) / R²(分组) / R²(组织)
  variance_decomp.csv        每个 PC 及每个基因的方差分解
  audit_verdict.json         总判定
  batch_audit.png            4 面板审计图

判定红线(与 SKILL.md 第 6 节一致):
  A. 最大方差 PC 的 R²(批次) > 0.9 且 R²(分组) < 0.2   → ❌ 不可用
  B. 任一分组 100% 落在单一批次                        → ❌ 该组作废
  C. 组内批次差 > 组间分组差                            → ❌ 数据不可用
  D. 文库大小双峰且与分组共变                           → ⚠️ 高度可疑,需批次校正后再看
"""
import argparse, io, json, os, sys
import numpy as np

np.set_printoptions(suppress=True)


# ---------------------------------------------------------------- IO
# STAR/salmon 等比对器常见后缀,会在列名里污染样本 ID(如 E2_H_1Aligned.sortedByCoord.out.bam)
SUFFIXES = [
    'Aligned.sortedByCoord.out.bam', 'Aligned.toTranscriptome.out.bam', 'Aligned.out.bam',
    '_counts.txt', '.counts', '.bam', '_quant', '.genes.results', '.isoforms.results',
]


def clean_name(c):
    c = c.strip().strip('"').strip()
    for s in SUFFIXES:
        if c.endswith(s):
            c = c[:-len(s)]
    return c


def read_matrix(path):
    sep = ',' if path.lower().endswith(('.csv', '.csv.gz')) else '\t'
    opener = io.open
    if path.endswith('.gz'):
        import gzip
        opener = lambda p, **k: gzip.open(p, 'rt', **k)
    with opener(path, encoding='utf-8', errors='replace') as f:
        header = f.readline().rstrip('\n').rstrip('\r')
        cols = [clean_name(c) for c in header.split(sep)]
        genes, rows = [], []
        for ln in f:
            p = ln.rstrip('\n').rstrip('\r').split(sep)
            if len(p) < 2:
                continue
            genes.append(p[0].strip().strip('"'))
            try:
                rows.append([float(x) if x not in ('', 'NA') else 0.0 for x in p[1:]])
            except ValueError:
                rows[-1] = [0.0] * (len(p) - 1)
    X = np.array(rows, dtype=float)          # gene x sample
    return genes, cols[1:], X


def read_samples(path):
    with io.open(path, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').rstrip('\r').split('\t')
        recs = []
        for ln in f:
            p = ln.rstrip('\n').rstrip('\r').split('\t')
            if len(p) < len(hdr):
                p += [''] * (len(hdr) - len(p))
            recs.append(dict(zip(hdr, p)))
    return recs


# ---------------------------------------------------------------- 统计
def r2_of(y, Xd):
    """用最小二乘算 R²(不含截距时自动加)"""
    if Xd is None or (hasattr(Xd, 'shape') and Xd.shape[1] == 0):
        return 0.0
    A = np.column_stack([np.asarray(Xd, dtype=float), np.ones(len(y))])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ c
    ss_res = float((resid ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0.0


def dummies(labels):
    """标签 → 哑变量矩阵(丢第一类做参照)"""
    uniq = sorted(set(labels))
    if len(uniq) < 2:
        return np.zeros((len(labels), 0)), uniq
    M = np.zeros((len(labels), len(uniq) - 1))
    for j, u in enumerate(uniq[1:]):
        M[:, j] = [1.0 if l == u else 0.0 for l in labels]
    return M, uniq


def bimodality(x):
    """用最优点二分 + 轮廓系数判断是否双峰。返回 (是否双峰, 阈值, 两侧均值)"""
    xs = np.sort(x)
    best, best_i = -1, None
    n = len(xs)
    for i in range(2, n - 2):
        a, b = xs[:i], xs[i:]
        sa, sb = a.std(), b.std()
        if sa + sb <= 0:
            continue
        # 类间距离 / 类内散布
        score = abs(a.mean() - b.mean()) / (sa + sb)
        if score > best:
            best, best_i = score, i
    if best_i is None:
        return False, float(np.median(x)), (0.0, 0.0)
    lo, hi = xs[:best_i], xs[best_i:]
    gap = hi[0] - lo[-1]
    span = xs[-1] - xs[0]
    # ⚠️ 只用"类间距离/类内散布"会把窄幅连续分布误判为双峰。
    #    真实建库批次差异通常在**倍数**级别:
    #      GSE276942  15.2–18.9M vs 21.8–29.4M
    #      GSE297195   6.7–9.4M  vs 65.2–102.9M
    #    故额外要求两簇均值比 > 1.6。
    ratio = hi.mean() / lo.mean() if lo.mean() > 0 else 1.0
    is_bi = (best > 1.5 and ratio > 1.6 and (span > 0 and gap / span > 0.02)
             and len(lo) >= 3 and len(hi) >= 3)
    return is_bi, float((lo[-1] + hi[0]) / 2), (float(lo.mean()), float(hi.mean()))


def binom_enrich(tab):
    """交叉表混杂:每行是否 100% 落在单列"""
    out = {}
    for g, row in tab.items():
        tot = sum(row.values())
        out[g] = {'n': tot, 'dist': row,
                  'confounded': any(v == tot for v in row.values()) if tot else False}
    return out


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--samples', required=True, help='tsv: sample/group[/batch/tissue]')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--batch-col', default='batch')
    ap.add_argument('--group-col', default='group')
    ap.add_argument('--tissue-col', default='tissue')
    ap.add_argument('--org', default='mouse', choices=['mouse', 'human', 'rat'])
    ap.add_argument('--gene-map', default=None,
                    help='ensembl/symbol/chrom 映射表(矩阵用 Ensembl ID 时必需)')
    ap.add_argument('--stratum', default=None,
                    help='分层变量列名(如 timepoint)。整库不可用时,用它判断哪一层内部仍然可用')
    ap.add_argument('--topn', type=int, default=5000, help='PCA/方差分解用高变基因数')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genes, mcols, X = read_matrix(args.matrix)
    srecs = read_samples(args.samples)

    # 按样本表顺序重排矩阵列
    idx, miss = [], []
    pos = {c: i for i, c in enumerate(mcols)}
    for r in srecs:
        # 矩阵列名经过 clean_name(去 .bam/.cel 等后缀),样本表必须同规则清洗,
        # 否则带后缀的样本名会全部 miss → n=0 → 下游 SVD 崩溃。
        s = clean_name(r['sample'])
        if s in pos:
            idx.append(pos[s])
        else:
            miss.append(s)
    if miss:
        print(f'⚠️ 样本表中有 {len(miss)} 个样本在矩阵里找不到,已跳过: {miss[:5]}...')
    srecs = [r for r in srecs if clean_name(r['sample']) in pos]
    X = X[:, idx]
    n = X.shape[1]
    print(f'矩阵 {len(genes)} 基因 × {n} 样本;样本表 {len(srecs)} 行')
    if n < 4:
        print('❌ 可匹配样本 < 4,无法审计。请检查样本表 sample 列与矩阵列名是否一致'
              '(注意矩阵列名会自动去除 .bam/.cel/.txt 等后缀)。')
        sys.exit(1)

    groups = [r.get(args.group_col, '') for r in srecs]
    batches = [r.get(args.batch_col, '') for r in srecs]
    tissues = [r.get(args.tissue_col, '') for r in srecs]

    # ---------- 1. 技术协变量 ----------
    lib = X.sum(axis=0)
    det = (X > 0).sum(axis=0)

    pre = {'mouse': 'mt-', 'human': 'MT-', 'rat': 'mt-'}[args.org]
    rpl = 'Rpl' if args.org != 'human' else 'RPL'
    rps = 'Rps' if args.org != 'human' else 'RPS'

    # 🔴 若矩阵用 Ensembl ID,按前缀找 mt-/Rpl 会**静默返回空**,
    #    协变量全为 0 且不报错 → 批次审计假阴性。此处显式处理。
    is_ensembl = sum(1 for g in genes[:500] if g.upper().startswith('ENS')) > 400
    sym = None
    bysym = {}
    if args.gene_map and os.path.exists(args.gene_map):
        sym = {}
        with io.open(args.gene_map, encoding='utf-8') as f:
            hdr = f.readline().rstrip('\n').split('\t')
            gi = hdr.index('ensembl') if 'ensembl' in hdr else 0
            si = hdr.index('symbol') if 'symbol' in hdr else 1
            ci = hdr.index('chrom') if 'chrom' in hdr else None
            for ln in f:
                p = ln.rstrip('\n').split('\t')
                if len(p) > max(si, gi):
                    s_ = p[si]
                    c_ = p[ci] if ci is not None and len(p) > ci else ''
                    sym[p[gi]] = (s_, c_)
                    if s_:
                        # ⚠️ 另一个静默失效: 如果待审计的矩阵**本身已经是符号矩阵**
                        #    (常见于先做了 ID→symbol 转换再复审计), 单向的
                        #    ensembl→symbol 映射会全部落空 → mt%/Rpl% 静默变 0。
                        #    所以这里同时建一张 符号→染色体 的反向表。
                        bysym[s_] = c_
        print(f'已加载基因映射 {len(sym)} 条 (符号索引 {len(bysym)} 条)')

    if is_ensembl and sym is None:
        print('⚠️ 矩阵使用 Ensembl ID 且未提供 --gene-map:'
              'mt-/核糖体协变量将无法计算(已置 0)。'
              '请先跑 ensembl_symbol_map.py,否则本审计缺少一类关键证据。')

    def frac(pref, chrom=None):
        ii = []
        for i, g in enumerate(genes):
            if sym is not None:
                s_, c_ = sym.get(g, ('', ''))
                if not s_ and not c_:
                    # 矩阵可能已经是符号矩阵 → 退回反向表
                    s_, c_ = g, bysym.get(g, '')
                if chrom and c_ == chrom:
                    ii.append(i)
                elif s_ and s_.startswith(pref):
                    ii.append(i)
            elif g.startswith(pref):
                ii.append(i)
        if not ii:
            return np.zeros(n), 0
        return X[ii].sum(axis=0) / np.maximum(lib, 1), len(ii)

    mt, n_mt = frac(pre, chrom='MT')
    rplf, n_rpl = frac(rpl)
    rpsf, n_rps = frac(rps)
    print(f'协变量基因数: 线粒体 {n_mt}  Rpl {n_rpl}  Rps {n_rps}')
    if n_mt == 0:
        print('⚠️ 线粒体基因数为 0,该协变量不可用')

    tech = {'sample': [r['sample'] for r in srecs], 'group': groups,
            'batch': batches, 'tissue': tissues, 'lib': lib, 'detected': det,
            'mt_frac': mt, 'rpl_frac': rplf, 'rps_frac': rpsf}
    with io.open(f'{args.outdir}/technical_covariates.csv', 'w', encoding='utf-8') as f:
        f.write('sample\tgroup\tbatch\ttissue\tlib_size\tdetected_genes\tmt_frac\trpl_frac\trps_frac\n')
        for i in range(n):
            f.write(f'{srecs[i]["sample"]}\t{groups[i]}\t{batches[i]}\t{tissues[i]}\t'
                    f'{lib[i]:.0f}\t{det[i]:.0f}\t{mt[i]:.5f}\t{rplf[i]:.5f}\t{rpsf[i]:.5f}\n')

    # ---------- 2. 文库大小双峰 & 推断批次 ----------
    is_bi, thr, (m_lo, m_hi) = bimodality(lib)
    inferred = ['HIGH' if l > thr else 'LOW' for l in lib]
    n_high = inferred.count('HIGH')

    # ---------- 3. 批次 × 分组 交叉表 ----------
    B = batches if any(batches) else inferred
    b_source = '实测' if any(batches) else '推断(文库大小双峰)'

    tab = {}
    for g in sorted(set(groups)):
        tab[g] = {}
    allb = sorted(set(B))
    for g in sorted(set(groups)):
        for b in allb:
            tab[g][b] = 0
    for g, b in zip(groups, B):
        tab[g][b] += 1
    conf = binom_enrich(tab)

    with io.open(f'{args.outdir}/confound_crosstab.csv', 'w', encoding='utf-8') as f:
        f.write('group\t' + '\t'.join(allb) + '\t total\tconfounded\n'.replace('\t ', '\t'))
        for g in sorted(tab):
            row = tab[g]
            f.write(g + '\t' + '\t'.join(str(row[b]) for b in allb)
                    + f'\t{conf[g]["n"]}\t{"YES" if conf[g]["confounded"] else "no"}\n')

    n_conf = sum(1 for g in conf if conf[g]['confounded'])

    # ---------- 4. PCA + 归属 ----------
    cpm = X / np.maximum(lib, 1) * 1e6
    L = np.log2(cpm + 1)
    v = L.var(axis=1)
    top = np.argsort(-v)[:min(args.topn, len(genes))]
    Lv = L[top]
    Z = (Lv - Lv.mean(axis=1, keepdims=True)) / (Lv.std(axis=1, keepdims=True) + 1e-9)

    Zc = Z - Z.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(Zc, full_matrices=False)
    ncomp = min(8, Vt.shape[0])
    pcs = Vt[:ncomp].T                       # sample x PC
    var_ratio = (S ** 2 / (S ** 2).sum())[:ncomp]

    Db, ub = dummies(B)
    Dg, ug = dummies(groups)
    Dt, ut = dummies(tissues) if any(tissues) else (np.zeros((n, 0)), [])

    pca_rows = []
    for k in range(ncomp):
        p = pcs[:, k]
        rb, rg, rt = r2_of(p, Db), r2_of(p, Dg), r2_of(p, Dt)
        both = r2_of(p, np.column_stack([Db, Dg]))
        pca_rows.append({'PC': k + 1, 'var_explained': var_ratio[k],
                         'R2_batch': rb, 'R2_group': rg, 'R2_tissue': rt,
                         'uniq_batch': both - rg, 'uniq_group': both - rb})
    with io.open(f'{args.outdir}/pca_attribution.csv', 'w', encoding='utf-8') as f:
        f.write('PC\tvar_explained\tR2_batch\tR2_group\tR2_tissue\tuniq_batch\tuniq_group\n')
        for r in pca_rows:
            f.write(f"{r['PC']}\t{r['var_explained']:.4f}\t{r['R2_batch']:.4f}\t{r['R2_group']:.4f}\t"
                    f"{r['R2_tissue']:.4f}\t{r['uniq_batch']:.4f}\t{r['uniq_group']:.4f}\n")

    # ---------- 5. 基因水平方差分解 ----------
    parts = []
    for i in range(len(top)):
        y = Lv[i]
        rb, rg = r2_of(y, Db), r2_of(y, Dg)
        both = r2_of(y, np.column_stack([Db, Dg]))
        parts.append((rb, rg, both - rg, both - rb))
    parts = np.array(parts)
    batch_dom = int((parts[:, 2] > parts[:, 3]).sum())
    gene_dom_t = int((parts[:, 3] > parts[:, 2]).sum())

    with io.open(f'{args.outdir}/variance_decomp.csv', 'w', encoding='utf-8') as f:
        f.write('gene\tR2_batch\tR2_group\tuniq_batch\tuniq_group\tdominant\n')
        for i, g in enumerate([genes[k] for k in top]):
            rb, rg, ub2, ug2 = parts[i]
            f.write(f'{g}\t{rb:.4f}\t{rg:.4f}\t{ub2:.4f}\t{ug2:.4f}\t'
                    f'{"batch" if ub2 > ug2 else "group"}\n')

    # ---------- 6. 对照实验:组内批次差 vs 组间分组差 ----------
    ctrl = []
    if len(set(B)) >= 2:
        for g in sorted(set(groups)):
            ii = [i for i in range(n) if groups[i] == g]
            if len(ii) < 4:
                continue
            blo = [i for i in ii if B[i] == allb[0]]
            bhi = [i for i in ii if B[i] == allb[-1]]
            if len(blo) >= 2 and len(bhi) >= 2:
                ctrl.append({'group': g, 'n_low': len(blo), 'n_high': len(bhi),
                             'within_batch_diff': float(Lv[:, bhi].mean() - Lv[:, blo].mean())})

    # ---------- 6b. 跨批次复现性(决定性检验) ----------
    # 若分组与批次可分离,则在两个批次里各自算一次"组间效应",
    # 二者应当高度相关 —— 这才是"数据可用"的实质证据,
    # 比"批次解释了多少方差"更有意义(批次解释得多但可校正,仍可用)。
    repl = {}
    gs = [g for g in sorted(set(groups)) if g != 'UNKNOWN']
    if len(allb) >= 2 and len(gs) >= 2:
        for a in range(len(gs)):
            for b in range(a + 1, len(gs)):
                g1, g2 = gs[a], gs[b]
                eff = {}
                ok = True
                for bt in allb:
                    i1 = [i for i in range(n) if groups[i] == g1 and B[i] == bt]
                    i2 = [i for i in range(n) if groups[i] == g2 and B[i] == bt]
                    if len(i1) < 2 or len(i2) < 2:
                        ok = False
                        break
                    eff[bt] = (Lv[:, i1].mean(axis=1) - Lv[:, i2].mean(axis=1))
                if not ok:
                    continue
                e1, e2 = eff[allb[0]], eff[allb[-1]]
                r = float(np.corrcoef(e1, e2)[0, 1])
                # 方向一致率
                same = float((np.sign(e1) == np.sign(e2)).mean())
                repl[f'{g1}_vs_{g2}'] = {
                    'n_genes': int(len(e1)), 'pearson_r': r, 'sign_agreement': same,
                    'n_batches': len(allb),
                }
    with io.open(f'{args.outdir}/cross_batch_replication.csv', 'w', encoding='utf-8') as f:
        f.write('contrast\tn_genes\tpearson_r\tsign_agreement\tn_batches\n')
        for k, v in repl.items():
            f.write(f"{k}\t{v['n_genes']}\t{v['pearson_r']:.4f}\t{v['sign_agreement']:.4f}\t{v['n_batches']}\n")

    # ---------- 7. 判定 ----------
    lead = max(pca_rows, key=lambda r: r['var_explained'])
    verdict = {
        'dataset': os.path.basename(args.matrix),
        'n_samples': n, 'n_genes': len(genes),
        'batch_source': b_source,
        'libsize_bimodal': bool(is_bi),
        'libsize_threshold': thr,
        'libsize_means': {'low': m_lo, 'high': m_hi},
        'n_high_batch': n_high,
        'confounded_groups': [g for g in conf if conf[g]['confounded']],
        'n_confounded_groups': n_conf,
        'lead_PC': {'PC': lead['PC'], 'var_explained': lead['var_explained'],
                    'R2_batch': lead['R2_batch'], 'R2_group': lead['R2_group']},
        'gene_level_batch_dominant': batch_dom,
        'gene_level_group_dominant': gene_dom_t,
        'within_group_batch_effects': ctrl,
    }

    fails, warn = [], []
    def loocv_of(idx, labels):
        """在给定样本子集内,用技术指标做留一最近质心分类,返回 (准确率, 随机基线)"""
        L = len(idx)
        if L < 4 or len(set(labels)) < 2:
            return 0.0, 0.0
        M = np.column_stack([np.log2(np.maximum(lib[idx], 1)), det[idx],
                             mt[idx], rplf[idx], rpsf[idx]])
        M = (M - M.mean(0)) / (M.std(0) + 1e-9)
        ok = 0
        for a in range(L):
            cents = {}
            for b in set(labels):
                m = [c for c in range(L) if c != a and labels[c] == b]
                if m:
                    cents[b] = M[m].mean(0)
            if cents:
                pred = min(cents, key=lambda b: float(np.linalg.norm(M[a] - cents[b])))
                ok += int(pred == labels[a])
        base = max(labels.count(u) for u in set(labels)) / L
        return ok / L, base

    def loocv_perm(idx, labels, B_perm=200, seed=0):
        """LOOCV 准确率的置换 p 值。小样本(n=10)下 0.70 vs 0.50 并不显著,
        必须用置换分布判断,否则会系统性过度报警。"""
        obs, base = loocv_of(idx, labels)

        def acc_of(lab):
            L = len(idx)
            M = np.column_stack([np.log2(np.maximum(lib[idx], 1)), det[idx],
                                 mt[idx], rplf[idx], rpsf[idx]])
            M = (M - M.mean(0)) / (M.std(0) + 1e-9)
            ok = 0
            for a in range(L):
                cents = {}
                for b in set(lab):
                    m = [c for c in range(L) if c != a and lab[c] == b]
                    if m:
                        cents[b] = M[m].mean(0)
                if cents:
                    pred = min(cents, key=lambda b: float(np.linalg.norm(M[a] - cents[b])))
                    ok += int(pred == lab[a])
            return ok / L

        rng = np.random.default_rng(seed)
        null = [acc_of(list(rng.permutation(labels))) for _ in range(B_perm)]
        p = float((np.array(null) >= obs).mean())
        return obs, base, p

    # ---------- 6c. 批次可信度:这些"批次"标签在技术指标上真的可区分吗? ----------
    # 背景:GSM 编号按分组连续排布时,batch 与 group 会"看起来"完全混杂;
    # 但编号习惯 ≠ 建库批次。用技术指标(文库/检出数/mt%/Rpl%/Rps%)做留一最近质心分类:
    # 若判别准确率不高于随机基线,说明这些块更可能是编号习惯,该"混杂"是虚假警报。
    #
    # ⚠️ 2026-09-18 打补丁:此前有两处会系统性误判
    #   (1) 总体判别只用了 "acc > base + 0.15" 的硬阈值,没做置换检验 ——
    #       而 6d 的层内判别用的是置换 p。口径不一致会导致"整体喊狼、层内说没事"。
    #       n=30、6 个块、每块 5 个样本时,0.40 vs 0.17 看着可观,置换 p 常 >0.05。
    #       现在统一为:acc > base + 0.15 且 置换 p < 0.05 才算"技术可分"。
    #   (2) 若"批次"块恰好**嵌套**在样本表里另一个变量内(如 B1/B2 全落在 17m),
    #       块间的技术差异可以完全由那个变量解释,与建库批次无关。
    #       不做嵌套检测就会把"年龄差异"误读成"批次效应",进而把可用数据判死。
    loocv, base_rate, p_overall = loocv_perm(list(range(n)), B)
    batch_plausible = bool(loocv > base_rate + 0.15 and p_overall < 0.05)

    # 嵌套检测:某个块的样本是否只出现在另一个变量的单一水平里
    #   注意要剔除两类"无信息"列,否则输出会被噪声淹没:
    #     ① 与批次分区**恒等**的列(如 group6 就是 B1..B6 的一一对应改名)——
    #        它当然"嵌套",但解释不了任何东西
    #     ② 与已收录列分区**重复**的列(如 age_ord / age_num / age 是同一分区)
    nesting = {}
    seen_parts = {}
    col_names = list(srecs[0].keys()) if srecs else []
    batch_part = frozenset(frozenset(i for i in range(n) if B[i] == b) for b in allb)
    for col in col_names:
        if col in (args.batch_col, 'sample'):
            continue
        vals = [r.get(col, '') for r in srecs]
        lvs = sorted(set(vals))
        if len(lvs) < 2 or len(lvs) >= n:
            continue
        part = frozenset(frozenset(i for i in range(n) if vals[i] == lv) for lv in lvs)
        if part == batch_part:
            continue                      # ① 与批次分区恒等
        if part in seen_parts:
            continue                      # ② 与已收录列分区重复
        seen_parts[part] = col
        if all(len({vals[i] for i in range(n) if B[i] == b}) == 1 for b in allb):
            nesting[col] = lvs
    # 由细到粗排序(层数多的在前),便于一眼看出最细的解释变量
    nesting = dict(sorted(nesting.items(), key=lambda kv: -len(kv[1])))

    # 对每个嵌套变量:在其每一层内部,还能用技术指标分出块吗?
    nested_tests = {}
    for col, lvs in nesting.items():
        vals = [r.get(col, '') for r in srecs]
        accs, ps = [], []
        for lv in lvs:
            ii = [i for i in range(n) if vals[i] == lv]
            if len(ii) < 6 or len({B[i] for i in ii}) < 2:
                continue
            a_, b_, p_ = loocv_perm(ii, [B[i] for i in ii])
            accs.append(a_); ps.append(p_)
        nested_tests[col] = {
            'n_levels': len(lvs), 'n_tested': len(accs),
            'mean_within_accuracy': float(np.mean(accs)) if accs else float('nan'),
            'min_within_perm_p': float(np.min(ps)) if ps else float('nan'),
            'within_separable': bool(accs and np.mean(accs) > 0.65 and np.min(ps) < 0.05),
        }

    # 块间技术差异被某个嵌套变量解释掉 → 不构成独立批次证据
    nested_explained = [c for c, t in nested_tests.items()
                        if t['n_tested'] > 0 and not t['within_separable']]

    verdict['batch_plausibility'] = {
        'loocv_accuracy_on_technical': loocv,
        'random_baseline': base_rate,
        'loocv_perm_p': p_overall,
        'batch_supported_by_technical': batch_plausible,
        'nested_within': nesting,
        'nested_within_tests': nested_tests,
        'explained_by_nesting': nested_explained,
    }
    if len(set(B)) > 1 and n_conf > 0 and not batch_plausible:
        warn.append(f'B*: 分组与"批次"标签看似 100% 混杂,但这些标签在技术指标上不可区分'
                    f'(留一准确率 {loocv:.2f} vs 基线 {base_rate:.2f}, 置换 p={p_overall:.3f})→ '
                    f'更可能是提交者的编号习惯而非真实建库批次,该混杂警报存疑')
    elif batch_plausible and nested_explained:
        t0 = nested_tests[nested_explained[0]]
        warn.append(f'B*: 批次标签在技术指标上可区分(准确率 {loocv:.2f} vs 基线 {base_rate:.2f}, '
                    f'置换 p={p_overall:.3f}),但它们完全嵌套在 {nested_explained} 的水平内,'
                    f'且各水平内部无法再用技术指标分出块 '
                    f'(层内平均准确率 {t0["mean_within_accuracy"]:.2f}, '
                    f'最小置换 p={t0["min_within_perm_p"]:.3f})→ '
                    f'块间技术差异可由 {nested_explained} 解释,不构成独立批次证据')
    elif batch_plausible:
        warn.append(f'批次标签在技术指标上可区分(留一准确率 {loocv:.2f} vs 基线 {base_rate:.2f}, '
                    f'置换 p={p_overall:.3f})且未被任何分层变量嵌套解释'
                    f'→ 批次真实存在,混杂判定成立')

    med_ub = float(np.median(parts[:, 2]))
    med_ug = float(np.median(parts[:, 3]))

    # A 数据主轴被批次占据
    if lead['R2_batch'] > 0.9 and lead['R2_group'] < 0.2:
        fails.append('A: 最大方差PC由批次主导(R2_batch>0.9 且 R2_group<0.2)')
    # B 分组与批次完全混杂(不可校正)
    #   —— 三层门槛层层收紧,避免把"编号习惯"或"分层变量差异"误判成批次混杂:
    #      ① 结构上确实一一对应(n_conf>0)
    #      ② 批次标签在技术指标上可区分(置换 p<0.05)
    #      ③ 该技术差异不是被别的分层变量嵌套解释掉的
    if n_conf > 0 and batch_plausible and not nested_explained:
        fails.append(f'B: {n_conf} 个分组 100% 属单一批次,不可分离')
    elif n_conf > 0 and batch_plausible and nested_explained:
        warn.append(f'B?: {n_conf} 个分组与批次标签完全重叠;批次标签虽在技术指标上可区分,'
                    f'但块间差异可由 {nested_explained} 解释 → 暂不计为红线,需人工确认建库记录')
    elif n_conf > 0:
        warn.append(f'B?: {n_conf} 个分组与批次标签完全重叠,但批次标签技术不可区分,'
                    f'暂不计为红线(需人工确认建库记录)')
    # C 批次系统性压过分组(不是"略多",而是量级压制)
    if batch_dom > 0.75 * len(parts) and med_ub > 2 * max(med_ug, 1e-6):
        fails.append(f'C: {batch_dom}/{len(parts)} 个高变基因由批次主导,'
                     f'且批次独有方差中位数({med_ub:.3f})为分组({med_ug:.3f})的 2 倍以上')
    # D 文库大小双峰且与分组共变
    if is_bi and n_conf > 0:
        warn.append('D: 文库大小双峰且与分组共变')

    # 跨批次复现性:只有能算出来的才算证据
    bad_repl = {k: v for k, v in repl.items() if v['pearson_r'] < 0.3 or v['sign_agreement'] < 0.70}
    if repl and len(bad_repl) == len(repl):
        fails.append(f'E: 跨批次复现性失败({len(bad_repl)} 组对比的批次间效应相关 r<0.3 或方向一致率<70%)')
    elif bad_repl:
        warn.append(f'E: {len(bad_repl)}/{len(repl)} 组对比的跨批次复现性不达标')
    elif repl:
        warn.append(f'E: 跨批次复现性通过({len(repl)} 组对比全部 r≥0.3 且方向一致率≥70%),'
                    f'但 n 小,不等于分组效应为真')

    # ---------- 6d. 分层可用性:整库"不可用"往往掩盖了局部可用 ----------
    # 例:GSE297195 的批次边界落在 5w / 24h 之间,而 24h 内部四组完全干净。
    # 报告每一层(stratum)内部的分组是否与批次混杂,才能给出可执行的结论。
    strata = {}
    if args.stratum and args.stratum in (srecs[0] if srecs else {}):
        for lv in sorted(set(r.get(args.stratum, '') for r in srecs)):
            ii = [i for i in range(n) if srecs[i].get(args.stratum) == lv]
            if len(ii) < 4:
                continue
            gl = [groups[i] for i in ii]
            bl = [B[i] for i in ii]
            t = {}
            for g in sorted(set(gl)):
                t[g] = {}
                for b in sorted(set(bl)):
                    t[g][b] = sum(1 for x, y in zip(gl, bl) if x == g and y == b)
            nconf = sum(1 for g in t if any(v == sum(t[g].values()) for v in t[g].values()))
            # 层内这些块在技术指标上可区分吗?不可区分 → 层内不存在真实批次,可用
            acc_s, base_s, p_s = loocv_perm(ii, bl)
            tech_sep = (acc_s > base_s + 0.15) and p_s < 0.05
            strata[lv] = {
                'n': len(ii), 'n_groups': len(set(gl)), 'n_batches': len(set(bl)),
                'table': t, 'n_confounded_groups': nconf,
                'loocv_accuracy': acc_s, 'random_baseline': base_s,
                'loocv_perm_p': p_s,
                'batch_technically_separable': bool(tech_sep),
                'usable_within': (nconf == 0 or len(set(bl)) == 1
                                  or not tech_sep),
            }
        verdict['strata'] = strata

    separable = (n_conf == 0 and len(set(B)) >= 2)
    # 批次标签与分组一一对应时,两个 R² 会恒等,uniq 恒为 0 —— 这是构造使然,
    # 不能当作"可分离"。此时唯一有效的证据是:批次标签在技术指标上是否可区分。
    degenerate = (n_conf == len([g for g in conf if conf[g]['n'] > 0]))

    if fails:
        rec = '不可用:见 FAIL 列表'
    elif degenerate and not batch_plausible:
        rec = ('可用,但批次身份待确认:分组与批次标签完全重叠属编号构造,'
               '且这些标签在技术指标上不可区分 → 无证据表明存在真实批次效应。'
               '建议先用 lib.size + detected genes 作协变量的敏感性分析,'
               '并向作者/补充材料核实建库批次记录')
    elif degenerate and nested_explained:
        rec = (f'可用(限制条件):分组与批次标签完全重叠属编号构造;批次标签虽在技术指标上可区分,'
               f'但该差异可由 {nested_explained} 解释,块内无独立批次证据。'
               f'建议:(1) 只做分层内对比;(2) 把 lib.size 等作协变量做敏感性分析;'
               f'(3) 用独立数据集复现;(4) 向作者核实建库批次记录')
    elif degenerate:
        rec = ('慎用:分组与批次标签完全重叠,且批次标签在技术指标上可区分、'
               '又不能被任何分层变量解释 → 从数据内部无法区分分组效应与批次效应,'
               '必须依赖独立数据集复现或作者提供建库记录')
    elif separable:
        rec = '可用:批次与分组可分离,把 batch 作为协变量纳入模型'
    else:
        rec = '慎用:批次与分组不可分离,需人工核查建库记录'

    verdict.update({
        'median_uniq_batch': med_ub, 'median_uniq_group': med_ug,
        'cross_batch_replication': repl,
        'batch_separable_from_group': bool(separable),
        'batch_block_degenerate': bool(degenerate),
        'recommendation': rec,
    })
    verdict['FAIL'] = fails
    verdict['WARN'] = warn
    verdict['usable'] = len(fails) == 0
    with io.open(f'{args.outdir}/audit_verdict.json', 'w', encoding='utf-8') as f:
        json.dump(verdict, f, ensure_ascii=False, indent=2)

    # ---------- 控制台报告 ----------
    print(f'\n{"="*78}')
    print(f'批次审计: {verdict["dataset"]}   n={n}   批次来源={b_source}')
    print(f'{"="*78}')
    print(f'文库大小: 双峰={is_bi}  阈值={thr/1e6:.1f}M  低簇均值={m_lo/1e6:.1f}M  高簇均值={m_hi/1e6:.1f}M')
    print(f'\n批次 × 分组 交叉表:')
    hdr = 'group'.ljust(16) + ''.join(b.rjust(10) for b in allb) + 'confounded'.rjust(14)
    print(hdr)
    for g in sorted(tab):
        row = ''.join(str(tab[g][b]).rjust(10) for b in allb)
        flag = '  ← ❌ 完全混杂' if conf[g]['confounded'] else ''
        print(g.ljust(16) + row + flag)
    print(f'\nPC 归属(前 {ncomp} 个主成分):')
    print('  PC  方差占比   R2(批次)  R2(分组)  批次独有  分组独有')
    for r in pca_rows:
        print(f"  {r['PC']:<4d}{r['var_explained']:9.3f}{r['R2_batch']:11.3f}{r['R2_group']:11.3f}"
              f"{r['uniq_batch']:10.3f}{r['uniq_group']:10.3f}")
    print(f'\n基因水平(前 {len(top)} 高变基因): 批次主导 {batch_dom}  分组主导 {gene_dom_t}')
    print(f'\n批次标签可信度: 留一准确率 {loocv:.2f} vs 随机基线 {base_rate:.2f}  '
          f'置换 p={p_overall:.3f}  → 技术可分辨: {"是" if batch_plausible else "否"}')
    if nesting:
        print('嵌套检测(块是否只落在另一个变量的单一水平内):')
        for col, lvs in nesting.items():
            t = nested_tests[col]
            if t['n_tested'] == 0:
                print(f'  {col:<18s} 嵌套于 {len(lvs)} 个水平 (层内样本不足,未做判别)')
                continue
            print(f"  {col:<18s} 嵌套于 {len(lvs)} 个水平 | 层内判别准确率 "
                  f"{t['mean_within_accuracy']:.2f} 最小置换p={t['min_within_perm_p']:.3f} "
                  f"→ {'层内仍可分辨' if t['within_separable'] else '层内不可分辨(差异由该变量解释)'}")
    if ctrl:
        print('\n组内批次差(对照实验):')
        for c in ctrl:
            print(f"  {c['group']:<16s} n(低批)={c['n_low']} n(高批)={c['n_high']}  "
                  f"组内批次差={c['within_batch_diff']:+.4f}")
    if repl:
        print('\n跨批次复现性(决定性检验):')
        for k, v in repl.items():
            flag = '✓' if (v['pearson_r'] >= 0.3 and v['sign_agreement'] >= 0.70) else '✗'
            print(f"  {flag} {k:<28s} r={v['pearson_r']:+.3f}  方向一致率={v['sign_agreement']:.2f}")
    if strata:
        print(f'\n分层可用性(按 {args.stratum}):')
        for lv, v in strata.items():
            subl = v['table']
            bl = sorted(set(b for g in subl for b in subl[g]))
            cells = []
            for g in sorted(subl):
                cells.append('/'.join(f'{b}:{subl[g][b]}' for b in bl if subl[g][b]))
            flag = '✓ 层内可用' if v['usable_within'] else '✗ 层内仍混杂'
            why = ('块间技术指标不可区分' if not v['batch_technically_separable']
                   else f"块间技术指标可区分(准确率{v['loocv_accuracy']:.2f} vs 基线"
                        f"{v['random_baseline']:.2f}, 置换p={v['loocv_perm_p']:.3f})")
            print(f"  {lv:<22s} n={v['n']:<3d} 分组{v['n_groups']} 批次{v['n_batches']}  "
                  f"{flag} ({why})")
            print(f"{'':6s}{'; '.join(cells)}")
    print(f"\n批次与分组可分离: {'是' if separable else '否'}"
          f"   (批次独有方差中位数 {med_ub:.3f} vs 分组 {med_ug:.3f})")
    print('\n' + '-' * 78)
    if fails:
        print('❌ 不可用,触发红线:')
        for x in fails:
            print('   -', x)
    else:
        print('✅ 未触发红线,数据集可用于分组比较')
    for w in warn:
        print('   ⚠️ ', w)
    print(f"\n建议: {verdict['recommendation']}")
    print(f'\n输出 → {args.outdir}/')


if __name__ == '__main__':
    main()
