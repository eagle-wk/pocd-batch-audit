# -*- coding: utf-8 -*-
"""
P1-1 统一 CLI —— 一条命令把一批 GEO 数据集跑完四防线审计并出齐报告

设计目标:
    把"下载(可选) → 解析样本表 → 四防线审计 → 逐数据集报告 → 总表 + 失真率统计"
    收敛成**一个入口**。P0-1 的 12 个跨疾病矩阵一下完,直接:
        python audit_one.py --input <放满 *_series_matrix.txt.gz 的目录> --outroot <输出根>
    即可全自动产出 15 份报告 + 一张总表 + 一个失真率 JSON,无需再手工拼。

三种输入模式(可混用):
    1) 离线(首选): --input 一个或多个文件/目录,递归找 *_series_matrix.txt(.gz) 与
       expression_matrix.tsv,就地审计。适合"数据已下好/已在本机"的场景,不依赖网络。
    2) 在线: --candidates <tsv>,逐行下载 series_matrix 再审计(走 geo_download_and_audit)。
    3) 仅汇总: --verdict-dir 指向若干已审计目录(含 audit_verdict.json),只并表不出新审计。

依赖: Python>=3.9;审计核心 batch_audit_generic.py 仅需 numpy(用 --audit-python 指定带 numpy 的解释器)。

输出(<outroot>/):
    <GSE>/audit/audit_verdict.json        四防线判定(由 batch_audit_generic 产出)
    <GSE>/REPORT.md                       逐数据集人类可读报告(本脚本生成)
    master_audit_table.tsv                所有数据集汇总(机器可读)
    master_audit_report.md                所有数据集汇总(人类可读 + 失真率)
    audit_summary.json                    失真率三档统计(直接喂论文 §3.3)
"""
import argparse, csv, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import geo_download_and_audit as G   # 复用: parse_series_matrix / write_dataset / detect_org / download_series_matrix


# ------------------------------------------------------------- 分类与统计
def _as_bool(x):
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in ('true', '1', 'yes')


def normalize_record(src, kind=None):
    """把 audit_verdict.json 或 seed tsv 行统一成带 'cat' 三档的记录。

    ⚠️ 判别键不能用 `'FAIL' in src`: seed tsv **也有 FAIL 列**,会把它误判成
       verdict 分支,导致 confounded_groups 被当字符串逐字符拆分、R2 列取空。
       改用 seed 独有的 n_WARN 列判别;调用方也可显式传 kind。
    """
    if kind is None:
        kind = 'seed' if 'n_WARN' in src else 'verdict'
    if kind == 'verdict':
        rec = {
            'dataset': src.get('dataset', ''),
            'n_samples': src.get('n_samples', ''),
            'n_genes': src.get('n_genes', ''),
            'batch_source': src.get('batch_source', ''),
            'confounded_groups': src.get('confounded_groups', []) or [],
            'n_confounded_groups': src.get('n_confounded_groups',
                                           len(src.get('confounded_groups', []) or [])),
            'batch_separable': src.get('batch_separable_from_group', ''),
            'degenerate': src.get('batch_block_degenerate', ''),
            'lead_PC_R2_batch': (src.get('lead_PC') or {}).get('R2_batch', ''),
            'lead_PC_R2_group': (src.get('lead_PC') or {}).get('R2_group', ''),
            'usable': _as_bool(src.get('usable', False)),
            'recommendation': src.get('recommendation', ''),
            'FAIL': src.get('FAIL', []) or [],
            'n_WARN': len(src.get('WARN', []) or []),
        }
    else:                                              # seed tsv 行
        row = src
        cg = [x for x in (row.get('confounded_groups') or '').split(';') if x]
        rec = {
            'dataset': row.get('dataset', ''),
            'n_samples': row.get('n_samples', ''),
            'n_genes': row.get('n_genes', ''),
            'batch_source': row.get('batch_source', ''),
            'confounded_groups': cg,
            'n_confounded_groups': int(row.get('n_confounded_groups') or 0),
            'batch_separable': row.get('batch_separable', ''),
            'degenerate': row.get('degenerate', ''),
            'lead_PC_R2_batch': row.get('lead_PC_R2_batch', ''),
            'lead_PC_R2_group': row.get('lead_PC_R2_group', ''),
            'usable': _as_bool(row.get('usable', False)),
            'recommendation': row.get('recommendation', ''),
            'FAIL': [x for x in (row.get('FAIL') or '').split(';') if x],
            'n_WARN': int(row.get('n_WARN') or 0),
        }
    # 三档: usable(干净) / caveat(可用但受限) / unusable(红线)
    r = rec['recommendation']
    if not rec['usable']:
        rec['cat'] = 'unusable'
    elif any(k in r for k in ['限制', '慎用', '待确认', '需人工', '只做分层', '不可分离', '天花板']):
        rec['cat'] = 'caveat'
    else:
        rec['cat'] = 'usable'
    return rec


# ------------------------------------------------------------- 逐数据集报告
def render_report_md(gse, verdict_path, audit_dir):
    with open(verdict_path, encoding='utf-8') as f:
        v = json.load(f)
    rec = normalize_record(v, kind='verdict')
    # 分组推断不可信 → 本报告的"不可用/可用"不成立,明确标 UNVERIFIED,
    # 避免读者把它当成真实裁定(否则会像 GSE163943 那样凭空多出一个假 unusable)。
    gi_path = os.path.join(os.path.dirname(audit_dir), 'group_inference.json')
    if os.path.exists(gi_path):
        try:
            with open(gi_path, encoding='utf-8') as f:
                gi = json.load(f)
            if not gi.get('verified', True):
                rec['cat'] = 'unverified'
        except Exception:
            pass
    L = []
    L.append(f'# Audit report — {gse}\n')
    L.append(f'- **Samples:** {rec["n_samples"]}  **Genes:** {rec["n_genes"]}  '
             f'**Batch source:** {rec["batch_source"]}')
    if rec['cat'] == 'unverified':
        L.append(f'- **Verdict:** `UNVERIFIED` — group labels could not be inferred '
                 f'(see group_inference.json); **do not use the FAIL/PASS below**.')
    else:
        L.append(f'- **Verdict:** `{rec["cat"].upper()}`  '
                 f'(usable flag = {rec["usable"]}, n_WARN = {rec["n_WARN"]})')
    if rec['FAIL'] and rec['cat'] != 'unverified':
        L.append(f'- **FAIL:** ' + '; '.join(rec['FAIL']))
    L.append('')

    crosstab = os.path.join(audit_dir, 'confound_crosstab.csv')
    if os.path.exists(crosstab):
        L.append('## Batch x group crosstab')
        with open(crosstab, encoding='utf-8') as f:
            for i, line in enumerate(f):
                L.append('| ' + ' | '.join(line.rstrip('\n').split('\t')) + ' |')
                if i == 0:
                    L.append('|' + '|'.join('---' for _ in line.rstrip('\n').split('\t')) + '|')
        L.append('')

    pca = os.path.join(audit_dir, 'pca_attribution.csv')
    if os.path.exists(pca):
        L.append('## PCA attribution (leading PCs)')
        with open(pca, encoding='utf-8') as f:
            for i, line in enumerate(f):
                L.append('| ' + ' | '.join(line.rstrip('\n').split('\t')) + ' |')
                if i == 0:
                    L.append('|' + '|'.join('---' for _ in line.rstrip('\n').split('\t')) + '|')
        L.append('')

    bp = v.get('batch_plausibility', {})
    if bp:
        L.append('## Batch-label plausibility')
        L.append(f'- LOOCV accuracy on technical covariates: '
                 f'{bp.get("loocv_accuracy_on_technical"):.3f} '
                 f'(random baseline {bp.get("random_baseline"):.3f}, '
                 f'permutation p = {bp.get("loocv_perm_p"):.3f})')
        L.append(f'- Batch supported by technical signal: {bp.get("batch_supported_by_technical")}')
        nesting = bp.get('nested_within', {})
        if nesting:
            L.append(f'- Nested within: ' + ', '.join(f'{k}({len(vals)} lv)' for k, vals in nesting.items()))
        expl = bp.get('explained_by_nesting', [])
        if expl:
            L.append(f'- Block difference explained by: {", ".join(expl)}')
        L.append('')

    repl = os.path.join(audit_dir, 'cross_batch_replication.csv')
    if os.path.exists(repl):
        L.append('## Cross-batch replication')
        with open(repl, encoding='utf-8') as f:
            for i, line in enumerate(f):
                L.append('| ' + ' | '.join(line.rstrip('\n').split('\t')) + ' |')
                if i == 0:
                    L.append('|' + '|'.join('---' for _ in line.rstrip('\n').split('\t')) + '|')
        L.append('')

    if rec['n_WARN']:
        L.append('## Warnings')
        for w in v.get('WARN', []):
            L.append(f'- CAUTION: {w}')
        L.append('')

    L.append('## Recommendation')
    L.append(rec['recommendation'])
    L.append('')
    out = os.path.join(os.path.dirname(audit_dir), 'REPORT.md')
    with open(out, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L))
    return out


# ------------------------------------------------------------- 输入发现 / 审计
def gse_of(path):
    m = re.search(r'(GSE\d+)', os.path.basename(path))
    if not m:                       # 退化命名(如 series_matrix.txt.gz 无前缀)→ 取父目录名
        m = re.search(r'(GSE\d+)', os.path.basename(os.path.dirname(path)))
    return m.group(1) if m else None


def discover_matrices(inputs):
    """按 GSE 去重: 同一数据集的 .txt 与 .txt.gz 只取一个(优先 .gz),避免重复审计。"""
    cand = {}
    for p in inputs:
        paths = []
        if os.path.isfile(p) and (p.endswith('.gz') or p.endswith('.txt') or p.endswith('.tsv')):
            paths.append(p)
        elif os.path.isdir(p):
            for root, _, files in os.walk(p):
                for fn in files:
                    if (('series_matrix' in fn and (fn.endswith('.txt') or fn.endswith('.txt.gz')))
                            or fn == 'expression_matrix.tsv'):
                        paths.append(os.path.join(root, fn))
        for x in paths:
            g = gse_of(x)
            if not g:
                continue
            if g not in cand or x.endswith('.gz'):
                cand[g] = x
    return sorted(cand.values())


def check_group_inference(sheet_path, min_known_frac=0.5):
    """分组推断可信度: UNKNOWN 过多或有效分组<2 → 判定不可信(审计结论不成立)。"""
    groups = []
    with open(sheet_path, encoding='utf-8') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for row in rd:
            groups.append((row.get('group') or '').strip())
    n = len(groups)
    unknown = sum(1 for g in groups if g == '' or g.upper() == 'UNKNOWN')
    real = [g for g in groups if g and g.upper() != 'UNKNOWN']
    n_groups = len(set(real))
    known_frac = (n - unknown) / n if n else 0.0
    verified = bool(n_groups >= 2 and known_frac >= min_known_frac)
    return {'n_samples': n, 'n_unknown': unknown, 'n_groups': n_groups,
            'known_frac': known_frac, 'verified': verified,
            'reason': ('分组推断可信' if verified else
                       f'UNKNOWN {unknown}/{n}, 有效分组 {n_groups} 个(<2 或已知比例<{min_known_frac})')}


def run_audit(gse, matrix_path, sheet_path, outroot, audit_script, org, py):
    audit_dir = os.path.join(outroot, gse, 'audit')
    os.makedirs(audit_dir, exist_ok=True)
    cmd = [py, audit_script, '--matrix', matrix_path, '--samples', sheet_path,
           '--outdir', audit_dir, '--org', org]
    subprocess.run(cmd, check=True, capture_output=True, text=True,
                   encoding='utf-8', errors='replace')
    return os.path.join(audit_dir, 'audit_verdict.json')


def process_matrix(path, outroot, audit_script, py, org_force=None, skip_done=False):
    gse = gse_of(path)
    if not gse:
        print(f'  x 无法解析 GSE: {path}')
        return None
    d = os.path.join(outroot, gse)
    verdict_path = os.path.join(d, 'audit', 'audit_verdict.json')
    if skip_done and os.path.exists(verdict_path):
        print(f'. 跳过(已审计) {gse}')
        return d
    print(f'\n=== {gse} ===')
    try:
        genes, sids, meta_list, mlines = G.parse_series_matrix(path)
    except Exception as e:
        print(f'  x 解析失败: {e}')
        return None
    if not sids:
        print(f'  x 未解析到样本列,跳过')
        return None
    if len(mlines) < 500:
        print(f'  x 表达数据行过少({len(mlines)} 行),可能为 RNA-seq 无 processed matrix,跳过审计 {gse}')
        return None
    matrix_path, sheet_path = G.write_dataset(gse, genes, sids, meta_list, mlines, outroot)
    # 🔴 分组推断可信度闸门:
    #    启发式分组常把样本全判成 UNKNOWN(如 GSE163943 元数据无对照关键词)。
    #    此时 R2_group 恒为 0,审计会把全部高变基因归因于批次,触发假 FAIL
    #    ("C: 5000/5000 个高变基因由批次主导")—— 这是推断失败,不是真批次问题。
    #    不拦住就会让一个假 unusable 污染总表与论文失真率。
    gi = check_group_inference(sheet_path)
    with open(os.path.join(d, 'group_inference.json'), 'w', encoding='utf-8') as f:
        json.dump(gi, f, ensure_ascii=False, indent=2)
    if not gi['verified']:
        print(f"  !! 分组推断不可信({gi['n_unknown']}/{gi['n_samples']} UNKNOWN, "
              f"有效分组 {gi['n_groups']} 个) → 本轮判定为 UNVERIFIED,不会覆盖人工核定的结论")
    org = org_force or G.detect_org(genes)
    try:
        vp = run_audit(gse, matrix_path, sheet_path, outroot, audit_script, org, py)
    except subprocess.CalledProcessError as e:
        print(f'  x 审计失败 {gse}: {e.stderr[:500]}')
        return None
    render_report_md(gse, vp, os.path.join(outroot, gse, 'audit'))
    print(f'  样本表: {sheet_path}  (请人工核对 group 列!)')
    return d


# ------------------------------------------------------------- 汇总
def aggregate(verdict_dirs, seed_tsv, outroot):
    from collections import Counter
    recs2 = []
    seen = set()
    unverified = []
    # 1) verdict 优先(信息更完整)
    for d in verdict_dirs:
        vj = os.path.join(d, 'audit', 'audit_verdict.json')
        if not os.path.exists(vj):
            vj = os.path.join(d, 'audit_verdict.json')
        if not os.path.exists(vj):
            continue
        with open(vj, encoding='utf-8') as f:
            r = normalize_record(json.load(f), kind='verdict')
        # 🔴 batch_audit_generic 写进 verdict 的 dataset 值是**矩阵文件名**
        #    (如 expression_matrix.tsv),不是 GSE 号。若照单全收,去重会失效
        #    (同数据集算两次),且与 seed tsv 的 dataset 对不上。这里强制取目录名。
        if not re.match(r'^GSE\d+$', str(r['dataset'])):
            holder = os.path.dirname(vj)
            r['dataset'] = (os.path.basename(os.path.dirname(holder))
                            if os.path.basename(holder) == 'audit'
                            else os.path.basename(holder))
        # 分组推断不可信 → 本轮判定不采用,留待人工核定的 seed 补位
        gi_path = os.path.join(d, 'group_inference.json')
        if os.path.exists(gi_path):
            try:
                with open(gi_path, encoding='utf-8') as f:
                    if not json.load(f).get('verified', True):
                        unverified.append(r['dataset'])
                        continue
            except Exception:
                pass
        if r['dataset'] in seen:
            continue
        recs2.append(r)
        seen.add(r['dataset'])
    # 2) seed 仅补缺(同 dataset 已有 verdict 则跳过);空 dataset 行丢弃
    if seed_tsv and os.path.exists(seed_tsv):
        with open(seed_tsv, encoding='utf-8-sig') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                ds = (row.get('dataset') or '').strip()
                if not ds or ds in seen:
                    continue
                recs2.append(normalize_record(row, kind='seed'))
                seen.add(ds)

    keys = ['dataset', 'n_samples', 'n_genes', 'batch_source', 'confounded_groups',
            'n_confounded_groups', 'batch_separable', 'degenerate',
            'lead_PC_R2_batch', 'lead_PC_R2_group', 'usable', 'cat',
            'recommendation', 'FAIL', 'n_WARN']
    tsv = os.path.join(outroot, 'master_audit_table.tsv')

    def _cell(v):
        # 列表字段(confounded_groups / FAIL)用 ; 连接写进 tsv,避免落出 python repr
        return ';'.join(str(x) for x in v) if isinstance(v, (list, tuple)) else v

    with open(tsv, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in recs2:
            w.writerow({k: _cell(r.get(k, '')) for k in keys})

    # 失真率三档
    cnt = Counter(r['cat'] for r in recs2)
    n = len(recs2)
    summary = {
        'n_datasets': n,
        'usable': cnt.get('usable', 0),
        'caveat': cnt.get('caveat', 0),
        'unusable': cnt.get('unusable', 0),
        'frac_usable': cnt.get('usable', 0) / n if n else 0,
        'frac_caveat': cnt.get('caveat', 0) / n if n else 0,
        'frac_unusable': cnt.get('unusable', 0) / n if n else 0,
        'frac_not_clean': (cnt.get('caveat', 0) + cnt.get('unusable', 0)) / n if n else 0,
        'unverified_group_inference': sorted(set(unverified)),
        'n_unverified': len(set(unverified)),
        'headline': (f'{cnt.get("unusable",0)}/{n} categorically unusable; '
                     f'{cnt.get("caveat",0)}/{n} usable only under documented restrictions; '
                     f'{cnt.get("usable",0)}/{n} unconditionally analyzable'),
    }
    with open(os.path.join(outroot, 'audit_summary.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # 人类可读汇总
    md = [f'# Master audit report — {n} datasets\n']
    md.append(f'**Distortion-rate summary:** {summary["headline"]}\n')
    if unverified:
        md.append('**Unverified group inference** (auto-inferred labels insufficient; '
                  'curated annotation used instead, auto verdict discarded): '
                  + ', '.join(sorted(set(unverified))) + '\n')
    md.append('| dataset | n | batch_source | confounded_groups | cat | usable | R2_batch | R2_group |')
    md.append('|---|---:|---|---|---|---:|---:|---:|')
    for r in sorted(recs2, key=lambda x: (x['cat'] != 'unusable', x['dataset'])):
        cg = ';'.join(r['confounded_groups']) if r['confounded_groups'] else '-'
        md.append(f"| {r['dataset']} | {r['n_samples']} | {r['batch_source']} | {cg} | "
                  f"{r['cat']} | {r['usable']} | {r['lead_PC_R2_batch']} | {r['lead_PC_R2_group']} |")
    md.append('')
    for r in recs2:
        if r['cat'] != 'usable':
            md.append(f"**{r['dataset']}** ({r['cat']}): {r['recommendation']}")
    md.append('')
    with open(os.path.join(outroot, 'master_audit_report.md'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(md))

    print(f'\n汇总 {n} 个数据集 -> {tsv}')
    print(f'失真率: {summary["headline"]}')
    return recs2, summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', nargs='+', help='离线: 含 series_matrix/expression_matrix 的文件或目录(可多个)')
    ap.add_argument('--candidates', default=None, help='在线: 候选池 tsv(逐行下载+审计)')
    ap.add_argument('--verdict-dir', nargs='+', default=None, help='仅汇总: 含 audit_verdict.json 的目录')
    ap.add_argument('--outroot', required=True)
    ap.add_argument('--audit-script', default=os.path.join(HERE, 'batch_audit_generic.py'))
    ap.add_argument('--audit-python', default=sys.executable,
                    help='带 numpy 的 python 解释器(默认当前解释器)')
    ap.add_argument('--seed-tsv', default=os.path.join(HERE, 'P0_1_seed_audits.tsv'))
    ap.add_argument('--org', default=None, choices=['mouse', 'human', 'rat'],
                    help='强制物种(否则按基因符号推断)')
    ap.add_argument('--proxy', default=None)
    ap.add_argument('--skip-done', action='store_true')
    args = ap.parse_args()

    if args.proxy:
        for k in ('http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY'):
            os.environ[k] = args.proxy

    os.makedirs(args.outroot, exist_ok=True)
    verdict_dirs = []

    # 模式 1: 在线下载 + 审计
    if args.candidates:
        with open(args.candidates, encoding='utf-8-sig') as f:
            cands = [r['accession'].strip() for r in csv.DictReader(f, delimiter='\t')
                     if r.get('accession', '').strip()]
        print(f'候选 {len(cands)} 个: {", ".join(cands)}')
        for gse in cands:
            if args.skip_done and os.path.exists(os.path.join(args.outroot, gse, 'audit', 'audit_verdict.json')):
                print(f'. 跳过(已审计) {gse}')
                verdict_dirs.append(os.path.join(args.outroot, gse))
                continue
            mtx = G.download_series_matrix(gse, args.outroot)
            if not mtx:
                continue
            d = process_matrix(mtx, args.outroot, args.audit_script, args.audit_python, args.org, args.skip_done)
            if d:
                verdict_dirs.append(d)

    # 模式 2: 离线审计
    if args.input:
        mats = discover_matrices(args.input)
        print(f'发现 {len(mats)} 个矩阵文件')
        for m in mats:
            d = process_matrix(m, args.outroot, args.audit_script, args.audit_python, args.org, args.skip_done)
            if d:
                verdict_dirs.append(d)

    # 模式 3: 仅汇总
    if args.verdict_dir:
        verdict_dirs += list(args.verdict_dir)

    if not verdict_dirs:
        print('! 没有任何可汇总的审计结果。请检查 --input / --candidates / --verdict-dir。')
        return

    aggregate(verdict_dirs, args.seed_tsv, args.outroot)
    print('\n完成。产物: master_audit_table.tsv / master_audit_report.md / audit_summary.json'
          ' / 各 <GSE>/REPORT.md')


if __name__ == '__main__':
    main()
