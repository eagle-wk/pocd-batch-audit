# -*- coding: utf-8 -*-
"""生成本地已审计 6 个 POCD 数据集的判定快照 P0_1_seed_audits.tsv。
4 个读 audit_verdict.json,2 个(276942/163943)用已确立结论补填。
仅依赖标准库,本机(无网)即可跑。"""
import csv, json, os

BASE = os.path.dirname(os.path.abspath(__file__))
POCD = os.path.abspath(os.path.join(BASE, '..', 'C_m6A_PCD_patterns'))

json_dirs = {
    'GSE330865': os.path.join(POCD, 'GSE330865_data', 'audit'),
    'GSE297195': os.path.join(POCD, 'GSE297195_data', 'audit'),
    'GSE283401': os.path.join(POCD, 'GSE283401_data', 'audit'),
    'GSE304956': os.path.join(POCD, 'GSE304956_data', 'audit'),
}

# 手工补填(无 audit_verdict.json,但基于已确立结论)
manual = {
    'GSE276942': {
        'dataset': 'GSE276942', 'n_samples': 17, 'n_genes': '',
        'batch_source': '实测', 'confounded_groups': 'all(批次≡时间点)',
        'n_confounded_groups': 5, 'batch_separable': False, 'degenerate': True,
        'lead_PC_R2_batch': 0.361, 'lead_PC_R2_group': 0.0,
        'usable': False,
        'recommendation': '不可用:批次与时间点 100% 混杂(PC1 36.1% 被批次占据),不可校正',
        'FAIL': '批次与时间序列 100% 混杂', 'n_WARN': 0,
    },
    'GSE163943': {
        'dataset': 'GSE163943', 'n_samples': 8, 'n_genes': '',
        'batch_source': '无批次标签', 'confounded_groups': '',
        'n_confounded_groups': 0, 'batch_separable': True, 'degenerate': False,
        'lead_PC_R2_batch': 0.0, 'lead_PC_R2_group': 0.0,
        'usable': True,
        'recommendation': '可用(限制条件): 人类外周血 4v4, n=8 已触及天花板,不能作为独立验证队列',
        'FAIL': '', 'n_WARN': 1,
    },
}

keys = ['dataset', 'n_samples', 'n_genes', 'batch_source', 'confounded_groups',
        'n_confounded_groups', 'batch_separable', 'degenerate',
        'lead_PC_R2_batch', 'lead_PC_R2_group', 'usable', 'recommendation', 'FAIL', 'n_WARN']

rows = []
for gse, d in json_dirs.items():
    vj = os.path.join(d, 'audit_verdict.json')
    if not os.path.exists(vj):
        print(f'!! 缺失 {vj}')
        continue
    with open(vj, encoding='utf-8') as f:
        v = json.load(f)
    rows.append({
        'dataset': gse,
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

for gse, m in manual.items():
    rows.append(m)

out = os.path.join(BASE, 'P0_1_seed_audits.tsv')
with open(out, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f'已写出 seed: {len(rows)} 行 → {out}')
for r in rows:
    print(f"  {r['dataset']:12s} n={str(r['n_samples']):>3} usable={r['usable']}")
