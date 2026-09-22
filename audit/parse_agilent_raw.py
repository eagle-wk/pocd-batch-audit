# -*- coding: utf-8 -*-
"""从 GEO 的 RAW tar 里构建表达矩阵 + 样本表(离线, 不联网)。

支持两种单通道芯片原始格式:
  A) Agilent Feature Extraction (表头行以 FEATURES 开头)
     关键列: ControlType(0=真实探针), GeneName / SystematicName, gProcessedSignal
     元数据: FEPARAMS 块含 Scan_Date / Scan_ScannerName → **真实批次**
  B) CodeLink / GenUs (表头行含 Probe_name, 无 FEPARAMS)
     关键列: Probe_name, Normalized_intensity, Probe_type(DISCOVERY=真实探针)

🔴 教训: 列名一律**按名字定位**, 不写死下标。
   GSE59645 有 GeneName, GSE115614 只有 SystematicName(RefSeq 号),
   GSE31357 是 CodeLink 全无 Agilent 列名 —— 写死下标会静默产出 0 基因而不报错。

用法:
  python parse_agilent_raw.py --tar X_RAW.tar --outdir prepared/GSEXXXX \
      [--include 'GSM.*_GE_'] [--group-regex '_GE_(?:Hippo_)?(.+)_24hr_\d+\.txt']
"""
import argparse, gzip, io, os, re, statistics, tarfile
from collections import Counter, defaultdict

GENE_CANDIDATES = ['GeneName', 'SystematicName', 'ProbeName', 'Probe_name']
SIGNAL_CANDIDATES = ['gProcessedSignal', 'Normalized_intensity',
                     'Signal_strength', 'Raw_intensity']

def split_rows(raw):
    with gzip.open(io.BytesIO(raw), 'rt', errors='replace') as f:
        return [ln.rstrip('\n').rstrip('\r') for ln in f]

def parse_meta(lines):
    for i, ln in enumerate(lines[:12]):
        if ln.startswith('FEPARAMS') and i + 1 < len(lines):
            return dict(zip(ln.split('\t'), lines[i + 1].split('\t')))
    return {}

def find_header(lines):
    """返回 (表头行下标, 列名 list) 或 (None, None)"""
    for i, ln in enumerate(lines[:40]):
        cols = ln.split('\t')
        if ln.startswith('FEATURES') or 'Probe_name' in cols or 'ProbeName' in cols:
            return i, cols
    return None, None

def sample_signal(raw):
    """返回 (gene -> [signal...], 实际使用的基因 id 列名)"""
    lines = split_rows(raw)
    hi, cols = find_header(lines)
    if hi is None:
        return {}, None
    idx = {c: j for j, c in enumerate(cols)}
    gcol = next((idx[c] for c in GENE_CANDIDATES if c in idx), None)
    scol = next((idx[c] for c in SIGNAL_CANDIDATES if c in idx), None)
    if gcol is None or scol is None:
        return {}, None
    if 'ControlType' in idx:                       # Agilent
        ci = idx['ControlType']
        keep = lambda p: p[ci] == '0'
    elif 'Probe_type' in idx:                      # CodeLink
        ci = idx['Probe_type']
        keep = lambda p: p[ci].upper() == 'DISCOVERY'
    else:
        keep = lambda p: True
    acc = defaultdict(list)
    for ln in lines[hi + 1:]:
        if not (ln.startswith('DATA') or re.match(r'^\d+\t', ln)):
            continue
        p = ln.split('\t')
        if len(p) <= max(gcol, scol):
            continue
        if not keep(p):
            continue
        g = p[gcol].strip()
        if not g or g.startswith(('GE_BrightCorner', 'DarkCorner')):
            continue
        try:
            acc[g].append(float(p[scol]))
        except ValueError:
            continue
    return acc, cols[gcol]

GROUP_RULES = [('Naive', 'naive'), ('Sham', 'sham'), ('TBI+PMI', 'tbi_pmi'),
               ('TBI+JM6', 'tbi_jm6'), ('TBI+E33', 'tbi_e33'), ('TBI', 'tbi')]

def group_of(fname, gref=None):
    base = os.path.basename(fname)
    if gref:
        m = re.search(gref, base)
        return m.group(1) if m else re.sub(r'[^A-Za-z0-9]+', '_', base)
    m = re.search(r'_GE_(.+?)\.txt', base)
    tag = m.group(1) if m else base
    for pat, lab in GROUP_RULES:
        if tag.startswith(pat):
            return lab
    return re.sub(r'[^A-Za-z0-9]+', '_', tag)

def pick_batch(dates, scanners):
    """优先实测扫描日期, 其次扫描仪; 每层次至少 2 个样本才当作可用批次"""
    for cand, label in ((dates, '实测:扫描日期(Scan_Date)'),
                        (scanners, '实测:扫描仪(Scan_ScannerName)')):
        c = Counter(x for x in cand if x)
        if len(c) >= 2 and min(c.values()) >= 2:
            return list(cand), label
    return [''] * len(dates), '无(无可用实测批次元数据)'

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tar', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--include', default=None, help='只处理成员名匹配该正则的文件')
    ap.add_argument('--group-regex', default=None, help='从文件名提取分组, 取第 1 捕获组')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    inc = re.compile(a.include) if a.include else None
    per_sample, metas, used_col = {}, {}, None
    with tarfile.open(a.tar) as t:
        members = [m for m in t.getmembers()
                   if m.isfile() and (inc is None or inc.search(m.name))]
        print(f'  选中 {len(members)} 个文件')
        for m in members:
            raw = t.extractfile(m).read()          # 只读一次, 元数据与表达共用
            lines = split_rows(raw)
            metas[m.name] = parse_meta(lines)
            acc, gcol = sample_signal(raw)
            used_col = used_col or gcol
            per_sample[m.name] = acc
            print(f'  · {m.name}: {len(acc)} 基因(id 列={gcol}), '
                  f'group={group_of(m.name, a.group_regex)}')
            del raw, lines
    if any(len(v) == 0 for v in per_sample.values()):
        print('  ⚠️ 有文件解析出 0 基因, 请检查表头识别')

    names = list(per_sample)
    common = None
    for n in names:
        s = set(per_sample[n])
        common = s if common is None else (common & s)
    common = sorted(common or [])
    print(f'\n共同基因(id={used_col}): {len(common)}')

    samples = [os.path.splitext(os.path.splitext(n)[0])[0] for n in names]
    groups = [group_of(n, a.group_regex) for n in names]
    dates = [metas[n].get('Scan_Date', '')[:10] for n in names]
    scanners = [metas[n].get('Scan_ScannerName', '') for n in names]
    batch, bsrc = pick_batch(dates, scanners)
    print(f'批次来源: {bsrc}')
    print('  日期分布:', dict(Counter(dates)))
    print('  分组分布:', dict(Counter(groups)))

    with io.open(os.path.join(a.outdir, 'expression_matrix.tsv'), 'w', encoding='utf-8') as f:
        f.write('gene\t' + '\t'.join(samples) + '\n')
        for g in common:
            vals = [f'{statistics.median(per_sample[n][g]):.4f}' for n in names]
            f.write(g + '\t' + '\t'.join(vals) + '\n')
    with io.open(os.path.join(a.outdir, 'samplesheet.tsv'), 'w', encoding='utf-8') as f:
        f.write('sample\tgroup\tbatch\n')
        for s, gg, b in zip(samples, groups, batch):
            f.write(f'{s}\t{gg}\t{b}\n')
    print(f'\n✓ 写出 {len(common)} 基因 × {len(samples)} 样本 → {a.outdir}')

if __name__ == '__main__':
    main()
