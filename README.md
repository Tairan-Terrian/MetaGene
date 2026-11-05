# AI主动学习驱动的“代谢-表观遗传”跨维度合成致死聚焦筛选与奥沙利铂耐药逆转

本课题面向结直肠癌奥沙利铂耐药，围绕代谢重编程与表观遗传重塑之间的组合脆弱性，开展AI主动学习引导的跨通路合成致死聚焦筛选，并通过后续功能验证评价耐药逆转。

**MetaGene是本课题的干实验计算流程。** 计算团队负责整合公开数据、构建预测模型、评价主动学习策略，并交付排序后的**Top500代谢酶×组蛋白修饰酶潜在合成致死基因组合**及逐对解释。医学院团队承接聚焦筛选和代表组合验证，评价组合效应与奥沙利铂敏感性变化，并回传结果供模型更新。当前干实验版本：**sl_v2_2026-11-05**。

整体研究路径：**多源数据整合 → AI预测与主动学习 → 基因组合交付 → 聚焦筛选与耐药逆转验证 → 结果回传与模型更新**。

## 干实验主要结果

**Top500预测协同方向：500/500对的预测dLFC<0，范围为-0.711至-0.310，中位数-0.543。** 该结果表示模型预测的双干预效应低于两单基因效应的加性预期。

- 对38个代谢酶与15个直接组蛋白修饰酶形成的550对未测候选排序，输出Top500、Top30、Priority8及3对代表组合。
- 每对Top500提供预测dLFC、排序分数、单基因状态、背景证据、来源链接和一段生物学解释。
- 12轮回顾性主动学习中，上界策略在2,582条测量记录下平均发现14个强协同信号，同预算随机策略约4.67个，是随机策略的约3倍。
- 已有运行记录汇总232次模型拟合及103项完整性检查。

## 阅读与交接入口

医学院直接访问[MetaGene GitHub仓库](https://github.com/Tairan-Terrian/MetaGene)即可阅读交付内容，无需运行模型。建议顺序：**本README → [实验汇报](results/MetaGene实验汇报.md) → [Top500](results/candidates/top500.csv)／[Top30](results/candidates/top30.csv)／[证据卡](results/candidates/candidate_evidence_cards.md) → [交接说明](docs/CANDIDATE_HANDOFF.md)／[回传模板](results/handoff/measured_feedback_template.csv)**。

| 文件 | 内容 |
|---|---|
| [实验汇报](results/MetaGene实验汇报.md) | 总课题目标与干实验阶段的方法、结果、候选及图表 |
| [Top500](results/candidates/top500.csv) | 排序后的潜在合成致死／协同抑制候选及逐对解释 |
| [Top30](results/candidates/top30.csv) / [证据卡](results/candidates/candidate_evidence_cards.md) | 用于首批候选讨论的多样化短名单 |
| [首轮3对代表组合](results/handoff/first_round_representatives.csv) | UCK2 × KDM3B、PKM × HDAC2、HK2 × KAT6B |
| [字段字典](results/candidates/data_dictionary.csv) | 分数、区间、单基因条件及来源字段 |
| [研究设计](docs/RESEARCH_DESIGN.md) / [交接说明](docs/CANDIDATE_HANDOFF.md) | 总课题路径、干实验四模块、团队分工与结果回传 |
| [审计记录](results/integrity_audit.json) | 标签、分割、对称性、排序等检查 |

## 如何读懂协同信号

中心化LFC下，`dLFC = LFC_AB − LFC_A − LFC_B`。**dLFC<0表示相对于加性预期的额外抑制，越负表示协同恶化越强。** 对原始未中心化数据，逐重复加回双对照基线。本项目以这一相对效应筛选候选，无须设定绝对致死水平。

模型的主分类头和分类评价使用**dLFC<−1的强协同事件**，连续回归头输出预测dLFC。Top500的预测均值处于−1与0之间，表示模型预测有协同抑制。`synergy_probability`对应前述强协同事件分数，`pred_dlfc`用于阅读预测协同方向与幅度，二者各有用途。

交付顺序综合单基因条件与协同分数。Top500中190对两臂均未触发强单基因耗竭，272对需补充单基因测量，38对作为强单基因背景下的协同抑制扩展；各类独立清单位于`results/candidates`。

## 研究与验证衔接

当前结果为HT29筛选数据上学习得到的模型预测。核心配对属于外推应用，主测试和校准分别有5个、3个强协同阳性；实际组合效应及HCT116-R/SW480-R中的表现通过合作方回传测量继续评价。详细数据覆盖与模型比较集中列在实验汇报中。

## 复现与存储

本GitHub仓库提供实现代码、配置、报告、图表、候选结果、原始与处理后数据、模型权重及逐次运行记录。网页可直接阅读报告和候选；完整克隆前安装Git和Git LFS：

```bash
git lfs install
git clone https://github.com/Tairan-Terrian/MetaGene.git
cd MetaGene
git lfs pull
```

在仓库根目录运行`python scripts/write_sl_documents.py`，可根据已有精简结果重新生成阅读文档；网页阅读无需运行任何程序。

数据、权重和运行记录随仓库提供，大文件由Git LFS管理；交接直接使用GitHub，不生成交接压缩包。设备自动选择、指定设备以及修改模型后的独立运行目录见[复现说明](docs/REPRODUCIBILITY.md)。安装[项目依赖](requirements.txt)后按该说明运行；参数见[冻结配置](configs/experiment_sl.json)，当前结果位于`results/sl`。
