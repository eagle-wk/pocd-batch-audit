# -*- coding: utf-8 -*-
"""
分层效应一致性检验 (stratified effect consistency)

动机:
  在 POCD 这类研究里,最有价值的信号不是"某个模块在整体上差异显著",
  而是**该效应在不同层之间是否可复现**。
  分层变量可以是: 建库批次 / 年龄 / 性别 / 时间点 / 中心。
  只在单一层出现的"显著"极可能是噪声;跨层同向且幅度递增才是真信号。

本脚本做什么:
  1. 对每个模块(或任意给定基因集)计算特征基因 ME
  2. 在**每个层(stratum)内部**做 处理组 vs 对照组 的效应量 Δ 与 Welch t 检验
  3. 汇总:
       n_same_sign      —— 多少层 Δ 与主方向一致
       pooled_delta     —— 层内标准化后的合并效应量
       p_pooled_perm    —— 层内置换检验 p(保持层结构不变,只打乱组标签)
       delta_trend_r    —— Δ 与**层顺序**(如年龄)的 Spearman 相关,
                           用置换检验判断"效应随层递增/递减"是否显著
  4. 判定档位:
       [稳健]   ≥2/3 层同向 且 合并效应置换 p<0.05
       [分层特异] 仅某一层显著(其余不显著)→ 不能外推,需明确标注
       [不一致] 方向相反 → 视为噪声

关键方法学要点:
  * 置换必须在**层内**做(permute group labels within stratum),
    否则会破坏层结构、产生假阳性。
  * 每层样本量小时(n<6),不要看单层 p 值,要看方向一致性。
  * Δ 用层内合并标准差标准化(Pooled SD),使不同层的效应量可比。

用法:
  python stratified_effect_consistency.py \
      --matrix cpm_log2.tsv --me me.tsv \
      --samples samplesheet.tsv --group-col group --stratum-col age \
      --case OP --ctrl Sham --stratum-order "3m,17m,27m" \
      --outdir consistency_results
"""
import argparse, csv, io, json, os
import numpy as np
from scipy.stats import ttest_ind, spearmanr

ap = argparse.ArgumentParser()
ap.add_argument('--matrix', required=True)
ap.add_argument('--me', default='', help='模块特征基因表(首列模块名, 其余列为各样本 ME);缺省则用 --genesets')
ap.add_argument('--genesets', default='', help='基因集文件(gene<TAB>set),对每个集合现算 ME')
ap.add_argument('--geneset-metric', choices=['eigengene', 'mean-z'], default='mean-z',
                help='用 --genesets 时如何把集合压成一个每样本分数。'
                     'eigengene=取第一主成分(基因数少时会过拟合噪声, 不推荐 <15 基因的集合); '
                     'mean-z=先按基因 z 标准化再取均值(默认, 对小集合稳健)')
ap.add_argument('--samples', required=True)
ap.add_argument('--group-col', default='group')
ap.add_argument('--stratum-col', default='age')
ap.add_argument('--sample-col', default='sample')
ap.add_argument('--case', required=True, help='处理组的 group 取值,如 OP / surgery')
ap.add_argument('--ctrl', required=True, help='对照组的 group 取值,如 Sham')
ap.add_argument('--stratum-order', default='', help='逗号分隔的层顺序,用于趋势检验(如 3m,17m,27m)')
ap.add_argument('--adjust', default='', help='逗号分隔的样本表数值列名。会把这些协变量(如文库量)'
                                            '从每个对象的 ME 里回归掉再做检验 —— '
                                            '这是判断"效应是不是技术因素驱动的"的关键敏感性分析')
ap.add_argument('--adjust-mode', choices=['within', 'global'], default='within',
                help='协变量校正方式。within=逐层各回归一次(默认, 推荐); '
                     'global=对全部样本统一回归一次。'
                     '注意: 当协变量本身与分层变量相关时(如"文库量随年龄变化"), '
                     'global 会用一个由组间变异主导的斜率去校正层内数据, 造成**过度校正**, '
                     '把真实的层内效应一起抹掉。除非有特殊理由, 一律用 within。')
ap.add_argument('--min-n', type=int, default=3, help='每层每组最少样本数')
ap.add_argument('--min-geneset', type=int, default=5, help='基因集最少匹配基因数')
ap.add_argument('--perm', type=int, default=5000)
ap.add_argument('--seed', type=int, default=42)
ap.add_argument('--outdir', default='consistency_results')
args = ap.parse_args()
rng = np.random.default_rng(args.seed)
os.makedirs(args.outdir, exist_ok=True)

# ---------- 样本表 ----------
srows = list(csv.DictReader(io.open(args.samples, encoding='utf-8'), delimiter='\t'))
samples = [r[args.sample_col].strip() for r in srows]
grp = np.array([r[args.group_col].strip() for r in srows])
strata = np.array([r[args.stratum_col].strip() for r in srows])
print(f'样本表: {len(samples)} 样本 | 组: {dict(zip(*np.unique(grp, return_counts=True)))}')
print(f'                   层: {dict(zip(*np.unique(strata, return_counts=True)))}')

# ---------- 表达矩阵 ----------
rows = list(csv.reader(io.open(args.matrix, encoding='utf-8'), delimiter='\t'))
hdr = rows[0][1:]
genes = [r[0] for r in rows[1:]]
X = np.array([[float(x) for x in r[1:]] for r in rows[1:]])
col = {s: j for j, s in enumerate(hdr)}
missing = [s for s in samples if s not in col]
if missing:
    raise SystemExit(f'矩阵缺少样本: {missing[:5]} ... 共 {len(missing)} 个')
order = [col[s] for s in samples]
X = X[:, order]
print(f'矩阵: {X.shape[0]} 基因 x {X.shape[1]} 样本')
gi = {g: i for i, g in enumerate(genes)}
upper = {g.upper(): g for g in genes}


def eigengene(idx_list):
    sub = X[idx_list]
    Z = (sub - sub.mean(axis=1, keepdims=True)) / (sub.std(axis=1, keepdims=True) + 1e-12)
    u, s, vt = np.linalg.svd(Z, full_matrices=False)
    pc = vt[0]
    if pc.sum() < 0:
        pc = -pc
    return pc


def mean_z(idx_list):
    """先按基因(行)z 标准化, 再对集合取均值。小集合(5~15 基因)比 eigengene 稳健得多。"""
    sub = X[idx_list]
    Z = (sub - sub.mean(axis=1, keepdims=True)) / (sub.std(axis=1, keepdims=True) + 1e-12)
    return Z.mean(axis=0)


def score_set(idx_list):
    return mean_z(idx_list) if args.geneset_metric == 'mean-z' else eigengene(idx_list)


# ---------- 待检验对象 ----------
targets = {}          # name -> ME 向量
if args.me:
    mrows = list(csv.reader(io.open(args.me, encoding='utf-8'), delimiter='\t'))
    mh = mrows[0][1:]
    mc = {s: j for j, s in enumerate(mh)}
    for r in mrows[1:]:
        targets[r[0]] = np.array([float(r[1 + mc[s]]) for s in samples])
elif args.genesets:
    sets = {}
    with io.open(args.genesets, encoding='utf-8') as f:
        rd = csv.reader(f, delimiter='\t')
        next(rd, None)
        for r in rd:
            if len(r) < 2:
                continue
            g, c = r[0].strip(), r[1].strip()
            if g and c:
                sets.setdefault(c, []).append(g)
    for name, gs in sets.items():
        idx = sorted({gi[g] for g in gs if g in gi} | {gi[upper[g.upper()]] for g in gs if g.upper() in upper})
        if len(idx) >= args.min_geneset:
            targets[name] = score_set(idx)
        else:
            print(f'   [跳过] {name}: 仅匹配 {len(idx)} 个基因 (低于 --min-geneset={args.min_geneset})')
else:
    raise SystemExit('必须给出 --me 或 --genesets')

print(f'待检验对象: {len(targets)} 个')

# ---------- 技术协变量校正 ----------
adj_cols = [c.strip() for c in args.adjust.split(',') if c.strip()]
A_global = None
if adj_cols:
    try:
        cols = [np.array([float(r[c]) for r in srows], float) for c in adj_cols]
    except KeyError as e:
        raise SystemExit(f'--adjust 指定的列不存在于样本表: {e}')
    print(f'技术协变量校正: {adj_cols}  方式 = {args.adjust_mode}')

    if args.adjust_mode == 'global':
        A_global = np.column_stack([np.ones(len(srows))] + cols)


def residualize(y):
    """按调整方式对 ME 去协变量。within 模式在每个层内部各拟合一次。"""
    if not adj_cols:
        return y
    if A_global is not None:
        beta, *_ = np.linalg.lstsq(A_global, y, rcond=None)
        return y - A_global @ beta
    out = y.copy()
    for st in set(strata):
        m = (strata == st)
        if m.sum() <= len(adj_cols) + 1:
            continue
        Xd = np.column_stack([np.ones(m.sum())] + [np.array([float(r[c]) for r in srows], float)[m]
                                                   for c in adj_cols])
        beta, *_ = np.linalg.lstsq(Xd, y[m], rcond=None)
        out[m] = y[m] - Xd @ beta
    return out


targets = {k: residualize(v) for k, v in targets.items()}
print()

# ---------- 层顺序 ----------
if args.stratum_order:
    SORDER = [s.strip() for s in args.stratum_order.split(',')]
else:
    SORDER = sorted(set(strata), key=lambda v: (len(v), v))
present = [s for s in SORDER if s in set(strata)]
print(f'层顺序: {present}\n')

# 数值层水平(用于趋势检验)
def level_of(s):
    try:
        return float(''.join(ch for ch in s if ch.isdigit() or ch == '.'))
    except Exception:
        return float(present.index(s))

# ---------- 核心 ----------
out = []
for name, me in targets.items():
    per = []
    for st in present:
        m = (strata == st)
        c = m & (grp == args.ctrl)
        p = m & (grp == args.case)
        if c.sum() < args.min_n or p.sum() < args.min_n:
            continue
        cv, pv = me[c], me[p]
        sd = np.sqrt((cv.var(ddof=1) + pv.var(ddof=1)) / 2)
        d = float(pv.mean() - cv.mean())
        dz = d / sd if sd > 0 else 0.0
        t, pval = ttest_ind(pv, cv, equal_var=False)
        per.append({'stratum': st, 'n_ctrl': int(c.sum()), 'n_case': int(p.sum()),
                    'mean_ctrl': float(cv.mean()), 'mean_case': float(pv.mean()),
                    'delta': d, 'delta_z': dz, 't': float(t), 'p': float(pval)})
    if not per:
        continue
    deltas = np.array([r['delta'] for r in per])
    dzs = np.array([r['delta_z'] for r in per])
    # 主方向 = 绝对值最大的效应方向
    main = np.sign(deltas[np.argmax(np.abs(deltas))])
    n_same = int(np.sum(np.sign(deltas) == main))
    # 合并效应(标准化后取均值 → 与层样本量无关)
    pooled_z = float(dzs.mean())

    # ---- 层内置换检验 ----
    # 保持每层的样本集合不变,只在层内重排"处理/对照"标签
    B = args.perm
    null = np.empty(B)
    stratum_masks = []
    for st in present:
        m = (strata == st)
        c = m & (grp == args.ctrl)
        p = m & (grp == args.case)
        if c.sum() < args.min_n or p.sum() < args.min_n:
            continue
        stratum_masks.append((np.where(c)[0], np.where(p)[0]))
    obs = abs(pooled_z)
    for b in range(B):
        dzb = []
        for ci, pi in stratum_masks:
            allidx = np.concatenate([ci, pi])
            perm = rng.permutation(allidx)
            nc = len(ci)
            a, bb = me[perm[:nc]], me[perm[nc:]]
            sd = np.sqrt((a.var(ddof=1) + bb.var(ddof=1)) / 2)
            dzb.append((bb.mean() - a.mean()) / sd if sd > 0 else 0.0)
        null[b] = abs(np.mean(dzb))
    p_perm = float((np.sum(null >= obs) + 1) / (B + 1))

    # ---- 趋势检验(效应是否随层水平递增) ----
    if len(per) >= 3:
        lv = np.array([level_of(r['stratum']) for r in per])
        rho, _ = spearmanr(lv, deltas)
        rho = float(rho) if np.isfinite(rho) else 0.0
        n_perm_t = min(B, 2000)
        cnt = 0
        for _ in range(n_perm_t):
            rr = spearmanr(lv, rng.permutation(deltas))[0]
            rr = float(rr) if np.isfinite(rr) else 0.0
            if abs(rr) >= abs(rho):
                cnt += 1
        p_trend = (cnt + 1) / (n_perm_t + 1)
    else:
        rho, p_trend = float('nan'), float('nan')

    # ---- 分档 ----
    n_sig = sum(1 for r in per if r['p'] < 0.05)
    if n_same == len(per) and p_perm < 0.05:
        tier = '稳健'
    elif n_same == len(per) and p_perm < 0.10:
        tier = '倾向一致'
    elif n_same >= max(2, len(per) - 1) and p_perm < 0.05:
        tier = '较稳健'
    elif n_sig >= 1 and n_sig < len(per):
        tier = '分层特异'
    elif n_same <= len(per) / 2 and len(per) >= 2:
        tier = '方向不一致'
    else:
        tier = '无信号'

    out.append({'target': name, 'n_strata': len(per), 'n_same_sign': n_same,
                'main_sign': int(main), 'pooled_delta_z': pooled_z,
                'p_pooled_perm': p_perm, 'n_strata_sig': n_sig,
                'trend_rho_vs_stratum': rho, 'trend_p': p_trend, 'tier': tier,
                'per_stratum': per})

out.sort(key=lambda r: (r['tier'] == '稳健' or r['tier'] == '较稳健',
                        abs(r['pooled_delta_z'])), reverse=True)

print('=' * 118)
print(f"{'对象':<22}{'层数':>5}{'同向':>5}{'合并Δz':>9}{'置换p':>9}{'显著层':>7}{'趋势rho':>9}{'趋势p':>8}  判定")
print('-' * 118)
for r in out:
    print(f"{r['target']:<22}{r['n_strata']:>5}{r['n_same_sign']:>5}{r['pooled_delta_z']:>+9.3f}"
          f"{r['p_pooled_perm']:>9.4f}{r['n_strata_sig']:>7}"
          f"{r['trend_rho_vs_stratum']:>+9.3f}{r['trend_p']:>8.3f}  [{r['tier']}]")

print('\n' + '=' * 118)
print('逐层效应明细:')
for r in out:
    if r['tier'] in ('无信号', '方向不一致') and abs(r['pooled_delta_z']) < 0.3:
        continue
    print(f"\n[{r['target']}]  {r['tier']}   (主方向 {'↑处理组' if r['main_sign']>0 else '↓处理组'})")
    for d in r['per_stratum']:
        star = ' *' if d['p'] < 0.05 else ''
        print(f"   {d['stratum']:<8} n={d['n_ctrl']}/{d['n_case']}  "
              f"对照={d['mean_ctrl']:+.3f} 处理={d['mean_case']:+.3f}  "
              f"Δ={d['delta']:+.3f} (Δz={d['delta_z']:+.2f})  p={d['p']:.4f}{star}")

# ---------- 输出 ----------
flat = [{k: v for k, v in r.items() if k != 'per_stratum'} for r in out]
if flat:
    with io.open(f'{args.outdir}/consistency_summary.csv', 'w', encoding='utf-8', newline='\n') as f:
        w = csv.DictWriter(f, fieldnames=list(flat[0].keys()))
        w.writeheader()
        for r in flat:
            w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    with io.open(f'{args.outdir}/consistency_per_stratum.csv', 'w', encoding='utf-8', newline='\n') as f:
        w = csv.writer(f)
        w.writerow(['target', 'tier', 'stratum', 'n_ctrl', 'n_case', 'mean_ctrl', 'mean_case', 'delta', 'delta_z', 't', 'p'])
        for r in out:
            for d in r['per_stratum']:
                w.writerow([r['target'], r['tier'], d['stratum'], d['n_ctrl'], d['n_case'],
                            f"{d['mean_ctrl']:.6g}", f"{d['mean_case']:.6g}",
                            f"{d['delta']:.6g}", f"{d['delta_z']:.6g}", f"{d['t']:.6g}", f"{d['p']:.6g}"])
json.dump({'case': args.case, 'ctrl': args.ctrl, 'strata': present,
           'n_perm': args.perm, 'results': out},
          io.open(f'{args.outdir}/consistency.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
print(f'\n已保存 -> {args.outdir}/consistency_summary.csv, consistency_per_stratum.csv, consistency.json')
