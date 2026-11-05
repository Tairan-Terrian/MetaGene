# 干实验交付说明：潜在合成致死基因组合

所属课题：**AI主动学习驱动的“代谢-表观遗传”跨维度合成致死聚焦筛选与奥沙利铂耐药逆转**。

本说明用于计算团队向医学院团队移交基因组合。计算团队交付排序候选、预测结果、生物学解释及结果回传接口；医学院团队衔接聚焦筛选和代表组合验证，评价协同效应及奥沙利铂敏感性变化。实测结果回传后用于模型评价和下一轮主动学习。

统一交付入口为[MetaGene GitHub仓库](https://github.com/Tairan-Terrian/MetaGene)，网页阅读无需运行模型。建议依次查看[README](../README.md)、[实验汇报](../results/MetaGene实验汇报.md)、[Top500](../results/candidates/top500.csv)／[Top30](../results/candidates/top30.csv)／[证据卡](../results/candidates/candidate_evidence_cards.md)，再按本说明使用[回传模板](../results/handoff/measured_feedback_template.csv)。完整数据、权重和运行记录按README中的Git LFS步骤获取。

## 本次交付

版本：sl_v2_2026-11-05。[Top500](../results/candidates/top500.csv)中500/500对预测dLFC<0，范围-0.711至-0.310，即模型预测双干预产生低于单基因加性预期的适应性变化。每对提供分数、单基因条件、背景证据、功能来源及生物学解释。

首批讨论使用[Top30](../results/candidates/top30.csv)、[逐对证据卡](../results/candidates/candidate_evidence_cards.md)与[3对代表组合](../results/handoff/first_round_representatives.csv)：UCK2 × KDM3B、PKM × HDAC2、HK2 × KAT6B。

## 字段与选择顺序

| 字段 | 含义 |
|---|---|
| pred_dlfc | 预测的加性残差；负值为协同抑制，越负额外效应越强 |
| synergy_probability | HT29 dLFC<−1强协同事件的校准分数 |
| delivery_rank | 综合单基因条件与协同分数的交付顺序 |
| synergy_rank | 纯分类模型分数顺序 |
| dlfc_pi90_lower / upper | 由校准残差得到的名义90%预测区间 |
| seed_score_sd / seed_dlfc_sd | 多次训练的优化随机性 |

完整说明见[字段字典](../results/candidates/data_dictionary.csv)。Top30在交付顺序上加入每个基因最多3次的多样性条件。

## 单基因分层与生物学讨论

Top500中190对两臂均未触发强单基因耗竭，272对需补充单基因测量，38对作为强单基因背景的协同抑制扩展。独立清单为[单基因初筛兼容](../results/candidates/single_compatible_candidates.csv)、[测量待补充](../results/candidates/single_evidence_missing_candidates.csv)和[强单基因扩展](../results/candidates/strong_single_synergy_extension.csv)。功能来源用于解释两臂各自的催化角色，背景数据用于讨论细胞模型的相关性。

本次交付为模型预测候选；合作方的双干预测量用于确认配对效应，匹配亲本与耐药背景用于评价耐药相关性。

## 容量安排与结果回传

[500配对靶点表](../results/handoff/500_pair_targets.csv)定义基因级对象；[1,000容量草案](../results/handoff/1000_slot_capacity_draft.csv)安排500个首设计、450个次设计和50个控制占位，guide设计由合作方衔接。

[回传模板](../results/handoff/measured_feedback_template.csv)记录配对键、两臂基因、细胞模型、亲本/耐药状态、药物条件、重复、单A/单B/双干预及对照LFC、标尺、质控和数据来源。未中心化LFC按每个重复计算`double − single_A − single_B + control`；已中心化数据的对照项为0。保留生物重复层级，用重复测量评价统计支持。

新数据首先用于评价本版预测，随后纳入下一轮学习。GitHub统一管理代码、配置、报告、图表、候选结果、数据、模型权重和运行记录；大文件使用Git LFS，交接直接使用仓库。
