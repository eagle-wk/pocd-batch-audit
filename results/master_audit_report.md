# Master audit report — 14 contrast-units

**Distortion-rate summary:** 6/14 categorically unusable; 6/14 usable only under documented restrictions; 2/14 unconditionally analyzable

| unit | n | batch_source | confounded_groups | cat | usable | R2_batch | R2_group |
|---|---:|---|---|---|---:|---:|---:|
| GSE276942 | 17 | design-derived | all(批次≡时间点) | unusable | False | 0.361 | 0.0 |
| GSE304956 | 29 | measured | - | unusable | False | 0.017822633741692684 | 0.031338428073278246 |
| GSE59645 | 16 | measured | naive;sham;tbi;tbi_e33;tbi_jm6;tbi_pmi | unusable | False | 0.9612295827438243 | 0.9959803247849396 |
| GSE64978 | 10 | inferred | - | unusable | False | 0.04353621830703225 | 0.1083889349259558 |
| GSE68207 | 8 | inferred | - | unusable | False | 0.8132851435187268 | 0.022081338905395098 |
| GSE80174_Hippocampus | 9 | inferred | injured | unusable | False | 0.6694424217728403 | 0.4684066416388716 |
| GSE115614 | 14 | measured | Sham;TBI;TBI+Imipramine;TBI+Prozac;TBI+Zoloft | caveat | True | 0.8200205687850002 | 0.9356142269009855 |
| GSE163943 | 8 | none | - | caveat | True | 0.0 | 0.0 |
| GSE283401 | 96 | inferred | old_control_48 hours_hemisphere;old_control_48 hours_hippocampus;old_control_6 hours_hemisphere;old_control_6 hours_hippocampus;old_surgery_48 hours_hemisphere;old_surgery_48 hours_hippocampus;old_surgery_6 hours_hemisphere;young_control_48 hours_hemisphere;young_control_48 hours_hippocampus;young_control_6 hours_hemisphere;young_surgery_48 hours_hemisphere;young_surgery_6 hours_hemisphere;young_surgery_6 hours_hippocampus | caveat | True | 0.049951083715533695 | 0.3279675526849526 |
| GSE297195 | 33 | GSM block | Aged_24h_sham;Aged_24h_surgery;Aged_5w_sham;Aged_5w_surgery;Young_24h_sham;Young_24h_surgery | caveat | True | 0.991335906628195 | 0.991335906628195 |
| GSE31357 | 32 | inferred | 24_hr_Control;24_hr_MT+TBI;24_hr_Saline+TBI;4_hr_CB+TBI;4_hr_Control;4_hr_MT+TBI;4_hr_Saline+TBI | caveat | True | 0.005549300125768486 | 0.4784646552457432 |
| GSE330865 | 30 | GSM block | 17m_sham;17m_surgery;27m_sham;27m_surgery;3m_sham;3m_surgery | caveat | True | 0.25642788985215714 | 0.25642788985215703 |
| GSE80174_Cortex | 10 | inferred | - | usable | True | 0.16062313410798623 | 0.683374231490931 |
| GSE80174_Thalamus | 10 | inferred | - | usable | True | 0.0034580564362242106 | 0.6868895991722899 |

**GSE59645** (unusable): 不可用:见 FAIL 列表
**GSE64978** (unusable): 不可用:见 FAIL 列表
**GSE68207** (unusable): 不可用:见 FAIL 列表
**GSE80174_Hippocampus** (unusable): 不可用:见 FAIL 列表
**GSE31357** (caveat): 慎用:批次与分组不可分离,需人工核查建库记录
**GSE115614** (caveat): 可用,但批次身份待确认:分组与批次标签完全重叠属编号构造,且这些标签在技术指标上不可区分 → 无证据表明存在真实批次效应。建议先用 lib.size + detected genes 作协变量的敏感性分析,并向作者/补充材料核实建库批次记录
**GSE330865** (caveat): 可用,但批次身份待确认:分组与批次标签完全重叠属编号构造,且这些标签在技术指标上不可区分 → 无证据表明存在真实批次效应。建议先用 lib.size + detected genes 作协变量的敏感性分析,并向作者/补充材料核实建库批次记录
**GSE297195** (caveat): usable under restriction: 6/6 groups single GSM block; technically separable (0.52 vs 0.21, perm p<0.001) but block difference fully explained by nesting in timepoint -> not independent batch evidence
**GSE283401** (caveat): 慎用:批次与分组不可分离,需人工核查建库记录
**GSE304956** (unusable): 不可用:见 FAIL 列表
**GSE276942** (unusable): 不可用:批次与时间点 100% 混杂(PC1 36.1% 被批次占据),不可校正
**GSE163943** (caveat): 可用(限制条件): 人类外周血 4v4, n=8 已触及天花板,不能作为独立验证队列
