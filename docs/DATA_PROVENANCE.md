# 数据来源与固定输入

MetaGene将数据分为三层：**公开研究与数据库来源 → 固定的处理后输入 → 从头训练与结果生成**。当前`sl_v2_2026-11-05`的默认复现从第二层开始：仓库提供本版所用的小型输入快照，重新构建标签、划分数据、训练模型并生成评价与候选结果。无需预训练权重、原计算服务器或重新下载大型公共数据库。运行步骤见[复现说明](REPRODUCIBILITY.md)。

## 1. 随仓库提供的九个固定输入

以下文件共同固定基因范围、特征行顺序、网络、标签来源及候选证据。文件大小和SHA256见[输入清单](../data/reference/input_manifest.json)，实验参数另见[冻结配置](../configs/experiment_sl.json)。

| 文件 | 数据角色与用途 |
|---|---|
| [HT29_gene_lfc.csv](../data/pairs_sources/HT29_gene_lfc.csv) | 作者公布的基因级LFC矩阵导出表，保留三个生物重复；用于重建逐重复加性残差与监督标签 |
| [HT29_annotations.csv](../data/pairs_sources/HT29_annotations.csv) | 作者实验设计与基因对应表，提供组合、单基因对照、方向及锚点身份 |
| [symbol_mapping.csv](../data/annotations/symbol_mapping.csv) | 本版HGNC符号与别名映射，用于统一标签、设计和基因目录 |
| [model_genes.csv](../data/processed/model_genes.csv) | 固定898个基因及其行顺序；特征矩阵与网络节点均使用这一顺序 |
| [integrated_features.npz](../data/evidence/integrated_features.npz) | 898个基因的1,831维处理后特征，整合通路、功能标记、依赖、表达、通量与组蛋白峰信息及缺失标记；训练时按各分割拟合标准化参数 |
| [graphs.npz](../data/evidence/graphs.npz) | 映射到本版基因目录的独立网络；多图模型使用Reactome、STRING和BioPlex网络 |
| [external_pairs_model_domain.csv](../data/evidence/external/external_pairs_model_domain.csv) | In4mer与Dede研究中落入本模型基因范围的实测配对，供独立评价使用 |
| [genes_annotation.csv](../data/evidence/genes_annotation.csv) | 基因功能、催化证据及代谢酶／直接组蛋白修饰酶分类，确定候选对象范围并支持解释 |
| [candidate_annotations.csv](../data/reference/candidate_annotations.csv) | 候选配对的固定注释与背景证据，用于衔接新训练所得预测、候选排序和逐对解释 |

`data/sl/ht29_labels.csv`由前两份HT29来源表结合符号映射重新生成；它是本次运行的产物。`integrated_features.npz`和`graphs.npz`是数据表示，模型参数则在每次训练中重新拟合。

## 2. 标签与基因范围的来源

HT29来源为[COLO1组合筛选研究（DOI：10.1038/s41467-025-67256-9）](https://doi.org/10.1038/s41467-025-67256-9)。作者公开的[基因LFC RDS](https://raw.githubusercontent.com/ibarrioh/DualGuide_COLO1/main/Analysis/ENCORE/input/colo1_plm_GEN.rds)和[设计注释RDS](https://raw.githubusercontent.com/ibarrioh/DualGuide_COLO1/main/Analysis/ENCORE/input/IDs_colo1.rds)分别导出为上述两份CSV。`prepare_pairs.py --convert-rds`保留数值矩阵、设计列名和R对象的列优先布局；`prepare_sl_labels.py`随后按重复和guide位置匹配单基因对照，并计算含共同对照校正的dLFC。

898基因目录来自既有公开组合数据的统一整理。Horlbeck数据采用[NAIAD的K562表](https://raw.githubusercontent.com/NeptuneBio/NAIAD/main/tutorials/data/horlbeck_k562_gamma.csv)与[Jurkat表](https://raw.githubusercontent.com/NeptuneBio/NAIAD/main/tutorials/data/horlbeck_jurkat_gamma.csv)，原始研究DOI为[10.1016/j.cell.2018.06.010](https://doi.org/10.1016/j.cell.2018.06.010)。Norman数据整理保留[gamma表](https://raw.githubusercontent.com/NeptuneBio/NAIAD/main/tutorials/data/norman_gamma.csv)和[GI表](https://raw.githubusercontent.com/NeptuneBio/NAIAD/main/tutorials/data/norman_gi.csv)的来源说明，原始研究DOI为[10.1126/science.aax4438](https://doi.org/10.1126/science.aax4438)。本版SL监督任务使用HT29标签；其他来源在基因目录构建与历史数据整理中的角色见[配对来源定义](../data/pairs_sources/source_definitions.json)。

## 3. 特征、网络与候选证据的来源

| 来源 | 本项目中的用途 | 公开入口与本地溯源记录 |
|---|---|---|
| HGNC、Reactome、EpiFactors | 统一基因符号、通路特征和代谢／表观遗传功能目录 | [HGNC完整表](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt)、[Reactome通路文件](https://reactome.org/download/current/ReactomePathways.gmt.zip)、[Reactome人类事件层级](https://reactome.org/ContentService/data/eventsHierarchy/9606)、[EpiFactors 2.1](https://epifactors.autosome.org/public_data/v2.1/EpiGenes_main.csv)；[注释来源清单](../data/annotations/source_manifest.json) |
| UniProt、GO与人工审定 | 审查催化功能，区分代谢酶、直接组蛋白修饰酶及其他调控因子 | [UniProt查询批次清单](../data/evidence/sources/uniprot_download_manifest.json)、[严格注释清单](../data/evidence/annotation_source_manifest.json)；EpiFactors论文DOI：[10.1093/nar/gkac989](https://doi.org/10.1093/nar/gkac989) |
| DepMap 24Q4 Public | 基因依赖及细胞背景特征 | DOI：[10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1)；[GeneEffect文件](https://ndownloader.figshare.com/files/51064667)、[Model文件](https://ndownloader.figshare.com/files/51065297) |
| TCGA-COAD、GSE132465 | 肿瘤／正常组织表达与CRC单细胞表达摘要 | [UCSC Xena表达文件](https://tcga.xenahubs.net/download/TCGA.COAD.sampleMap/HiSeqV2.gz)、[GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)；处理方式见[整合输入说明](../data/evidence/integrated_input_manifest.json) |
| STRING 12.5、BioPlex | 独立网络特征 | [STRING 12.5映射接口](https://version-12-5.string-db.org/api/tsv/get_string_ids)、[BioPlex HCT116网络](https://bioplex.hms.harvard.edu/data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv)、[BioPlex 293T网络](https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv)；BioPlex论文DOI：[10.1016/j.cell.2021.04.011](https://doi.org/10.1016/j.cell.2021.04.011) |
| HCT116/OxR代谢模型 | 作者模型的通量可行范围，经生长归一化和基因映射后形成特征 | DOI：[10.3390/cancers13164130](https://doi.org/10.3390/cancers13164130)；[作者仓库固定提交](https://github.com/HAHerrmann/Hct116_DrugRes/tree/3ba911ee3b931e0909add599beba48bdca9f69a9) |
| GSE138654组蛋白图谱 | H3K4me3、H3K27ac、H3K27me3峰的基因级描述 | [GSE138654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138654)，样本GSM4115600–GSM4115605；与通量一起记录在[通量与表观组清单](../data/evidence/flux_epigenome_manifest.json) |
| HT29单基因辅助证据 | 候选单基因效应分层 | 作者公开的[control-scaled效应](https://raw.githubusercontent.com/ibarrioh/DualGuide_COLO1/main/Analysis/ENCORE/input/colo1_plm_scale.rds)与[SCORE单基因LFC](https://raw.githubusercontent.com/ibarrioh/DualGuide_COLO1/main/Analysis/ENCORE/input/SCORE_median_lfc.rds)；[处理说明](../data/evidence/feature_manifest.json) |
| 耐药表达及单基因药物筛选 | 候选背景证据与功能讨论 | GSE76092、GSE299428、GSE289766、GSE42387、GSE254136；[当前耐药分析来源](../results/disease/source_manifest.json)与[早期数据整理记录](../data/disease/source_manifest.json) |
| In4mer与Dede | 独立配对评价 | 论文DOI：[10.1038/s41467-024-47795-3](https://doi.org/10.1038/s41467-024-47795-3)、[10.1186/s13059-020-02173-2](https://doi.org/10.1186/s13059-020-02173-2)；[外部来源清单](../data/evidence/external/source_manifest.json)记录Figshare文件44435639、44435633及各来源的命中定义 |

DepMap、TCGA、单细胞、通量、组蛋白峰、hg19参考注释及STRING的25份获取记录集中在[公共输入获取清单](../data/evidence/requested_sources/acquisition_manifest.json)，包含文件名、URL、字节数、SHA256和获取时间。该清单用于追溯特征来源；默认训练直接使用九个固定输入。

## 4. 人工审定与版本选择

[代谢功能审定表](../data/evidence/metabolic_adjudication.csv)和[染色质功能审定表](../data/evidence/chromatin_adjudication.csv)随仓库保留。它们记录催化功能、分类理由及引用，是定义研究对象的可审查决策输入。`annotate_strict.py`结合这些表与公开注释生成`genes_annotation.csv`。

本版使用Reactome 98、EpiFactors 2.1、UniProt 2026_03及GO的2026-07-26记录。动态接口会继续更新：2026-11-05审查时，UniProt已返回2026_04，其首批内容与本版快照不同。因此，当前版本训练固定使用已处理的输入；研究者更新公共注释时，应建立新的数据版本并保存新的来源与哈希。

## 5. 研究者可以修改什么

- **比较模型与训练策略**：修改模型代码或实验配置，保留输入快照，重新训练并比较评价结果。具体运行与输出约定见[复现说明](REPRODUCIBILITY.md)。
- **修改候选对象或生物学解释**：先审查两份人工决策表及`genes_annotation.csv`，再更新候选注释。基因范围变化时，同步维护基因目录、特征行顺序和网络节点映射。
- **加入新实测组合**：保留单A、单B、双干预、对照、细胞背景和重复信息，按本版标签定义处理；先评价当前预测，再纳入新版本训练。
- **更新公共数据**：复用下列上游脚本与来源清单，保存下载校验值和处理参数，生成新输入快照后重新训练。

## 6. 上游处理脚本

| 脚本 | 对应处理环节 |
|---|---|
| `prepare_annotations.py` | HGNC符号映射、Reactome通路、EpiFactors目录 |
| `prepare_pairs.py`、`assemble_data.py` | 作者RDS转换、配对整理、基因域和基础通路表示 |
| `annotate_strict.py` | 结合公开注释和人工决策定义代谢酶与组蛋白修饰酶范围 |
| `prepare_features.py` | BioPlex网络与作者单基因辅助证据 |
| `acquire_requested_sources.py` | 获取DepMap、TCGA、单细胞、通量、组蛋白峰和STRING等公共输入 |
| `prepare_flux_epigenome.py`、`prepare_integrated_inputs.py` | 生成多组学基因摘要、整合特征和网络 |
| `prepare_disease.py`、`prepare_resistance.py` | 耐药背景表达与药物筛选证据 |
| `prepare_external.py` | 独立In4mer／Dede组合评价数据 |
| `prepare_sl_labels.py` | 从固定HT29来源表重建本版dLFC标签 |

这些脚本保留了各环节的处理方法。完整上游更新需要按来源准备文件、衔接处理顺序并选择注释版本；当前默认复现入口直接从固定输入开始，从头训练模型并重新生成本版结果。
