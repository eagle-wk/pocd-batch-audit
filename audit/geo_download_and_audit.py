# -*- coding: utf-8 -*-
"""
P0-1 联网补齐 —— 一键下载 GEO 系列矩阵 + 解析样本注释 + 四防线审计 + 汇总

设计目标:
    在**有网络出口**的机器上跑一次,即可把"候选池"里的每个 GSE 下载下来、
    自动解析出样本分组、跑 batch_audit_generic 的四防线审计,最后把所有
    audit_verdict.json 汇总成一张总表(可与本地已审计的 6 个 POCD 数据集合并)。

依赖: Python >= 3.9 + numpy (batch_audit_generic 只依赖 numpy,本脚本只依赖标准库)
   联网: 需要能访问 https://ftp.ncbi.nlm.nih.gov/geo/ (GEO 矩阵)

用法:
    # 1) 仅下载 + 审计候选池(不含已审计的 6 个)
    python geo_download_and_audit.py \
        --candidates P0_1_cross_disease_candidates.tsv \
        --outroot ./p0_1_downloads \
        --seed-dir /path/to/already_audited   # 可选:本地已审计的 6 个 verdict 目录

    # 2) 也可先重扫 GEO 候选池(需另行跑 geo_rescan.py)再喂进来

输出:
    <outroot>/<GSE>/expression_matrix.tsv     基因 x 样本 表达矩阵
    <outroot>/<GSE>/samplesheet.tsv           sample/group[/batch/tissue]
    <outroot>/<GSE>/audit/audit_verdict.json  四防线判定
    <outroot>/master_audit_table.tsv          所有数据集汇总(含 seed)

⚠️ 自动分组映射是启发式,产出后请人工核对 samplesheet.tsv 里的 group 列。
   脚本会把"无法判定"的样本标成 UNKNOWN(审计时该组会被自动排除或标红)。
"""
import argparse, csv, io, json, os, re, subprocess, sys, time, urllib.parse, urllib.request, gzip

GEO_FTP = 'https://ftp.ncbi.nlm.nih.gov/geo/series'

# ------------------------------------------------------------- 分组推断词典
CONTROL_TOKENS = [
    r'\bsham\b', r'\bnaive\b', r'\bnaïve\b', r'\bcontrol\b', r'\bcontrols\b',
    r'\bvehicle\b', r'\bwild[- ]?type\b', r'\bwt\b', r'\bhealthy\b',
    r'\bbaseline\b', r'\buntreated\b', r'\bunexposed\b', r'\bsaline\b',
    r'\bmock\b', r'\bnormal\b', r'\bctr\b', r'\bcon\b',
]
CASE_TOKENS = [
    r'\bsurgery\b', r'\bsurgical\b', r'\banesthes', r'\banaesthes',
    r'\bsevo', r'\bisoflur', r'\bdesflur', r'\bpropofol\b',
    r'\binjur', r'\btbi\b', r'\blesion', r'\bdisease\b', r'\bpatient',
    r'\bad\b', r'\bcase\b', r'\bexposed\b', r'\btreated\b', r'\bmodel\b',
    r'\bpostop', r'\bpostop', r'\bpost-operative', r'\bdelirium\b',
    r'\bcognitive deficit', r'\btrauma\b', r'\baud\b', r'\balzheimer',
]
CONTROL_RE = re.compile('|'.join(CONTROL_TOKENS), re.I)
CASE_RE = re.compile('|'.join(CASE_TOKENS), re.I)

# batch / 技术批次字段
BATCH_FIELD_RE = re.compile(r'(^|[^a-z])(batch|run|lane|lib|library|plate|array|flowcell|date|sequencing|replicate)[^a-z]', re.I)
# tissue 字段
TISSUE_RE = re.compile(r'(hippocamp|CA1|CA3|dentate|cortex|prefrontal|brain|spinal|striatum|hypothalam|cerebell|microgl|neuron|glia)', re.I)


def clean_name(c):
    return c.strip().strip('"').strip()


def _try_download(url, dst, timeout):
    """尝试下载 url 到 dst;成功(>1KB)返回 True,否则 False。"""
    # 1) urllib: 会读取 http_proxy / https_proxy 环境变量
    try:
        print(f'  ↓ 下载 {url}')
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=timeout) as r, open(dst, 'wb') as f:
            f.write(r.read())
        if os.path.getsize(dst) > 1000:
            return True
    except Exception as e:
        print(f'    !! urllib: {e}')
    # 2) curl fallback: Windows 自带 curl 通常能走系统代理,比 urllib 更稳
    try:
        print(f'  ↓ curl fallback')
        proxy = os.environ.get('https_proxy') or os.environ.get('HTTPS_PROXY') or \
                os.environ.get('http_proxy') or os.environ.get('HTTP_PROXY')
        for curl in ['curl', 'curl.exe']:
            try:
                cmd = [curl, '-fsSL', '--max-time', str(timeout)]
                if proxy:
                    cmd += ['-x', proxy]
                cmd += ['-o', dst, url]
                res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout+5)
                if res.returncode == 0 and os.path.exists(dst) and os.path.getsize(dst) > 1000:
                    print(f'    · curl 成功')
                    return True
                else:
                    err = (res.stderr or res.stdout or '').strip()
                    if err:
                        print(f'    !! curl({curl}): {err[:200]}')
            except Exception as e:
                print(f'    !! curl({curl}): {e}')
                continue
    except Exception as e:
        print(f'    !! curl: {e}')
    return False


def download_series_matrix(gse, outdir, timeout=120, retries=3):
    """下载 GEO series matrix,多源 fallback(ftp / download 镜像 / https API / 非压缩)。
    返回本地路径;全部失败返回 None 并打印网络提示。"""
    pre = gse[:-3].lower() if len(gse) > 3 else gse.lower()
    base = f'{gse}.series_matrix.txt'
    cand = [
        f'{GEO_FTP}/{pre}/{gse}/matrix/{base}.gz',
        f'https://download.ncbi.nlm.nih.gov/geo/series/{pre}/{gse}/matrix/{base}.gz',
        f'https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse}&format=file&file={base}.gz',
        f'{GEO_FTP}/{pre}/{gse}/matrix/{base}',
    ]
    os.makedirs(outdir, exist_ok=True)
    for url in cand:
        dst = os.path.join(outdir, url.split('?')[0].rsplit('/', 1)[-1])
        if os.path.exists(dst) and os.path.getsize(dst) > 1000:
            print(f'  · 已存在,跳过下载 {dst}')
            return dst
        for attempt in range(retries):
            if _try_download(url, dst, timeout):
                return dst
            time.sleep(2 + attempt)
    print(f'  ✗ 所有下载源均失败 {gse} —— 请检查本机能否访问 NCBI(国内常被限,'
          f'可换学术网络/代理;或手动下载 {gse}_series_matrix.txt.gz 放到 '
          f'{os.path.abspath(outdir)}/ 后重跑)')
    return None


def parse_series_matrix(path):
    """解析 GEO series matrix.
    返回 (genes, sample_ids, meta_list, matrix_lines)
      genes: 基因 ID 列表(第一列)
      sample_ids: 样本列名(不含 ID_REF)
      meta_list: 每样本 dict(字段->值),来自 !Sample_* 元信息
      matrix_lines: 表达矩阵原始行(已含基因 ID 在首列)
    """
    opener = gzip.open if path.endswith('.gz') else open
    mode = 'rt' if path.endswith('.gz') else 'r'
    meta = {}          # 字段名 -> list(按样本列)
    sample_ids = None
    genes, matrix_lines = [], []
    in_table = False
    with opener(path, mode, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if line.startswith('!') and 'series_matrix_table_begin' not in line:
                # 元信息行: !Sample_xxx = val1 \t val2 ...
                m = re.match(r'^!(\S+?)\s*=\s*(.*)$', line)
                if m:
                    key, vals = m.group(1), m.group(2)
                    meta[key] = [clean_name(v) for v in vals.split('\t')]
                continue
            if 'series_matrix_table_begin' in line:
                in_table = True
                continue
            if 'series_matrix_table_end' in line:
                in_table = False
                continue
            if in_table:
                parts = line.split('\t')
                if sample_ids is None:
                    sample_ids = [clean_name(p) for p in parts[1:]]
                    continue
                genes.append(parts[0])
                matrix_lines.append(line)
    # 对齐 meta 到 sample_ids 列序
    n = len(sample_ids) if sample_ids else 0
    meta_list = []
    for j in range(n):
        d = {}
        for k, vals in meta.items():
            if k.startswith('!Sample_'):
                short = k[len('!Sample_'):]
                d[short] = vals[j] if j < len(vals) else ''
        meta_list.append(d)
    return genes, sample_ids, meta_list, matrix_lines


def text_of(meta):
    """把一个样本的元信息拼成一段可搜索文本"""
    return ' \n '.join(f'{k}={v}' for k, v in meta.items() if v)


def infer_group(meta):
    t = text_of(meta)
    is_control = bool(CONTROL_RE.search(t))
    is_case = bool(CASE_RE.search(t))
    if is_control and not is_case:
        return 'control'
    if is_case and not is_control:
        return 'case'
    if is_control and is_case:
        # 同时命中(罕见): 优先看 title/source_name 里的明确对照词
        return 'control' if CONTROL_RE.search(meta.get('title', '') + ' ' + meta.get('source_name_ch1', '')) else 'case'
    return 'UNKNOWN'


def infer_tissue(meta):
    t = text_of(meta)
    m = TISSUE_RE.search(t)
    return m.group(1) if m else ''


def infer_batch(meta):
    """从 characteristics 里找显式 batch/run/lane 等字段值;找不到返回 ''"""
    for k, v in meta.items():
        if BATCH_FIELD_RE.search(k) and v:
            # 只取像批次编号的值(含数字或字母数字)
            if re.search(r'[0-9]', v) or re.match(r'^[A-Za-z]\d', v):
                return v
        # characteristics 里常写成 "batch: 1" / "Batch: B1"
        if 'characteristics_ch1' in k:
            mm = re.search(r'batch\s*[:=]\s*(\S+)', v, re.I)
            if mm:
                return mm.group(1)
    return ''


def write_dataset(gse, genes, sample_ids, meta_list, matrix_lines, outdir):
    """写出 expression_matrix.tsv + samplesheet.tsv,返回 (matrix_path, sheet_path)"""
    d = os.path.join(outdir, gse)
    os.makedirs(d, exist_ok=True)
    matrix_path = os.path.join(d, 'expression_matrix.tsv')
    with open(matrix_path, 'w', encoding='utf-8') as f:
        f.write('ID_REF\t' + '\t'.join(sample_ids) + '\n')
        for ln in matrix_lines:
            f.write(ln + '\n')

    sheet_path = os.path.join(d, 'samplesheet.tsv')
    with open(sheet_path, 'w', encoding='utf-8') as f:
        f.write('sample\tgroup\tbatch\ttissue\n')
        for sid, meta in zip(sample_ids, meta_list):
            g = infer_group(meta)
            b = infer_batch(meta)
            t = infer_tissue(meta)
            f.write(f'{sid}\t{g}\t{b}\t{t}\n')
    return matrix_path, sheet_path


def detect_org(genes):
    """粗略判断物种: 人符号常全大写(RPL, ACTB),小鼠/大鼠首字母大写其余小写(Actb, Rpl)。"""
    head = [g for g in genes[:300] if g and not g.startswith('ENS')]
    if not head:
        return 'mouse'
    up = sum(1 for g in head if g == g.upper() and g.lower() != g.upper())
    return 'human' if up / len(head) > 0.5 else 'mouse'


def run_audit(gse, matrix_path, sheet_path, outdir, audit_script, org):
    audit_dir = os.path.join(outdir, gse, 'audit')
    os.makedirs(audit_dir, exist_ok=True)
    cmd = [sys.executable, audit_script,
           '--matrix', matrix_path, '--samples', sheet_path,
           '--outdir', audit_dir, '--org', org]
    print(f'  → 审计 {gse} (org={org})')
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True,
                       encoding='utf-8', errors='replace')
        return os.path.join(audit_dir, 'audit_verdict.json')
    except subprocess.CalledProcessError as e:
        print(f'    !! 审计失败 {gse}: {e.stderr[:500]}')
        return None


def summarize_combined(json_dirs, seed_tsv, out_tsv):
    """读若干 audit_verdict.json + 可选 seed tsv,汇总成一张总表。
    seed tsv 的列名与本表 keys 完全一致,可直接并入。"""
    rows = []
    for d in json_dirs:
        vj = os.path.join(d, 'audit_verdict.json')
        if not os.path.exists(vj):
            continue
        with open(vj, encoding='utf-8') as f:
            v = json.load(f)
        rows.append({
            'dataset': os.path.basename(d.rstrip('/')),
            'n_samples': v.get('n_samples', ''),
            'n_genes': v.get('n_genes', ''),
            'batch_source': v.get('batch_source', ''),
            'confounded_groups': ';'.join(v.get('confounded_groups', [])),
            'n_confounded_groups': v.get('n_confounded_groups', ''),
            'batch_separable': v.get('batch_separable_from_group', ''),
            'degenerate': v.get('batch_block_degenerate', ''),
            'lead_PC_R2_batch': v.get('lead_PC', {}).get('R2_batch', ''),
            'lead_PC_R2_group': v.get('lead_PC', {}).get('R2_group', ''),
            'usable': v.get('usable', ''),
            'recommendation': (v.get('recommendation', '') or '')[:160],
            'FAIL': ';'.join(v.get('FAIL', [])),
            'n_WARN': len(v.get('WARN', [])),
        })
    # 合并 seed: 列名与本表 keys 对齐(首行 BOM 用 utf-8-sig 去除)
    if seed_tsv and os.path.exists(seed_tsv):
        with open(seed_tsv, encoding='utf-8-sig') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                rows.append(dict(r))
    keys = ['dataset', 'n_samples', 'n_genes', 'batch_source', 'confounded_groups',
            'n_confounded_groups', 'batch_separable', 'degenerate',
            'lead_PC_R2_batch', 'lead_PC_R2_group', 'usable', 'recommendation', 'FAIL', 'n_WARN']
    with open(out_tsv, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})
    print(f'\n汇总 {len(rows)} 个数据集 → {out_tsv}')
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', required=True,
                    help='候选池 tsv(列: accession[, taxon, n_samples, platform, tissue, notes, cluster])')
    ap.add_argument('--outroot', required=True)
    ap.add_argument('--audit-script',
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          'batch_audit_generic.py'),
                    help='batch_audit_generic.py 路径(默认与本脚本同目录)')
    ap.add_argument('--seed-dir', default=None,
                    help='可选: 本地已审计数据集目录(含若干 audit_verdict.json),合并进总表')
    ap.add_argument('--seed-tsv',
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         'P0_1_seed_audits.tsv'),
                    help='已审计 POCD 数据集判定快照 tsv(默认同目录 P0_1_seed_audits.tsv),合并进总表')
    ap.add_argument('--skip-done', action='store_true', help='跳过已有 audit_verdict.json 的数据集')
    ap.add_argument('--only', default=None, help='只处理指定 accession(逗号分隔),用于补跑')
    ap.add_argument('--proxy', default=None,
                    help='代理地址,例如 http://127.0.0.1:59407;会同时传给 urllib 与 curl')
    args = ap.parse_args()
    if args.proxy:
        for k in ('http_proxy', 'https_proxy', 'HTTP_PROXY', 'HTTPS_PROXY'):
            os.environ[k] = args.proxy
        print(f'已启用代理: {args.proxy}')

    os.makedirs(args.outroot, exist_ok=True)
    # 读候选
    cands = []
    taxon_map = {}
    with open(args.candidates, encoding='utf-8-sig') as f:
        rd = csv.DictReader(f, delimiter='\t')
        for row in rd:
            acc = row.get('accession', '').strip()
            if acc:
                cands.append(acc)
                taxon_map[acc] = (row.get('taxon', '') or '').strip()
    cands = [c for c in cands if c]
    if args.only:
        cands = [c for c in cands if c in set(a.strip() for a in args.only.split(','))]
    print(f'候选 {len(cands)} 个: {", ".join(cands)}')

    verdict_dirs = []
    if args.seed_dir:
        for root, _, files in os.walk(args.seed_dir):
            if 'audit_verdict.json' in files:
                verdict_dirs.append(root)

    for gse in cands:
        d = os.path.join(args.outroot, gse)
        verdict_path = os.path.join(d, 'audit', 'audit_verdict.json')
        if args.skip_done and os.path.exists(verdict_path):
            print(f'· 跳过(已审计) {gse}')
            verdict_dirs.append(d)
            continue
        print(f'\n=== {gse} ===')
        mtx_gz = download_series_matrix(gse, args.outroot)
        if not mtx_gz:
            print(f'  ✗ 下载失败,跳过 {gse}')
            continue
        try:
            genes, sids, meta_list, mlines = parse_series_matrix(mtx_gz)
        except Exception as e:
            print(f'  ✗ 解析失败: {e}')
            continue
        if not sids:
            print(f'  ✗ 未解析到样本列,跳过')
            continue
        if len(mlines) < 500:
            print(f'  ✗ 表达数据行过少({len(mlines)} 行), 可能为 RNA-seq 无 processed matrix, '
                  f'跳过审计 {gse}')
            continue
        matrix_path, sheet_path = write_dataset(gse, genes, sids, meta_list, mlines, args.outroot)
        taxon = taxon_map.get(gse, '')
        org = {'Homo sapiens': 'human', 'Mus musculus': 'mouse',
               'Rattus norvegicus': 'rat'}.get(taxon)
        if not org:
            org = detect_org(genes)
        vp = run_audit(gse, matrix_path, sheet_path, args.outroot, args.audit_script, org)
        if vp:
            verdict_dirs.append(d)
        # 友好提示
        print(f'  样本表: {sheet_path}  (请人工核对 group 列!)')

    out_tsv = os.path.join(args.outroot, 'master_audit_table.tsv')
    summarize_combined(verdict_dirs, args.seed_tsv, out_tsv)
    print('\n完成。产出已落在本机工作区 (p0_1_downloads/master_audit_table.tsv),'
          'WorkBuddy 可直接读取并入汇总表。')


if __name__ == '__main__':
    main()
