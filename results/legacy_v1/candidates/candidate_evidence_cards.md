# 01. PANK1 + HDAC1

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 95；多样化 Top30 排名 1。Priority8：1。

## 计算预测与已测单基因

- 预测加性 GI：**-0.54**；上游 90% GI 预测区间 [-0.734, -0.347]；跨随机种子 SD 0.0299。
- 预测联合效应 -0.0645 = 已测单基因和 0.476 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HDAC1=0.163；PANK1=0.313，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HDAC1=0.0883；PANK1=0.274。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HDAC1（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt Q13547；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/Q13547/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 10220385。
- 作者control-scaled single 0.0883，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.00838；HCT116 -0.151；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 12.4；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.634。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.44/3.28；H3K27ac 1.61/1.42；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.162，FDR 0.105；SW480高耐药对亲本：log2FC -0.153，FDR 0.0474。
- GSE76092 log2FC -0.35，95%均值对比区间 [-0.442, -0.258]，BH FDR 0.000133，BY FDR 0.0014；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 14.9，强度分位数 0.977；耐药株均值 14.6，强度分位数 0.97。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.105 / FDR 0.511；GSE289766 HCT116 log2FC -0.866；GSE42387 HCT116 log2FC -0.212、HT29 log2FC -0.412。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 100 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes](https://reactome.org/content/detail/R-HSA-9615017)：CAMERA Direction=Down，FDR=0.564，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Estrogen-dependent gene expression](https://reactome.org/content/detail/R-HSA-9018519)：CAMERA Direction=Down，FDR=0.564，nGenes=116；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### PANK1（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q8TE04；EC 2.7.1.33；[来源1](https://www.uniprot.org/uniprotkb/Q8TE04/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.274，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.0626；HCT116 -0.0532；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.92；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.128。
- HCT116/OxR作者模型通量：映射反应数 5；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.0894、OxR 0.188。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.52/3.23；H3K27ac 1.79/2.02；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.0719，FDR 0.746；SW480高耐药对亲本：log2FC -0.567，FDR 0.000238。
- GSE76092 log2FC 0.403，95%均值对比区间 [0.252, 0.554]，BH FDR 0.00125，BY FDR 0.0132；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 10.3，强度分位数 0.736；耐药株均值 10.7，强度分位数 0.774。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.229 / FDR 0.335；GSE289766 HCT116 log2FC 0.452；GSE42387 HCT116 log2FC -0.0718、HT29 log2FC 0.191。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC -3.11；有效 guide 数 2；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Metabolism of vitamins and cofactors](https://reactome.org/content/detail/R-HSA-196854)：CAMERA Direction=Down，FDR=0.968，nGenes=195；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of water-soluble vitamins and cofactors](https://reactome.org/content/detail/R-HSA-196849)：CAMERA Direction=Up，FDR=0.975，nGenes=124；ORA上调集合FDR=1、下调集合FDR=1。
- [Vitamin B5 (pantothenate) metabolism](https://reactome.org/content/detail/R-HSA-199220)：CAMERA Direction=Up，FDR=0.99，nGenes=20；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 02. PANK1 + HDAC2

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 109；多样化 Top30 排名 2。Priority8：2。

## 计算预测与已测单基因

- 预测加性 GI：**-0.531**；上游 90% GI 预测区间 [-0.724, -0.338]；跨随机种子 SD 0.0181。
- 预测联合效应 -0.326 = 已测单基因和 0.205 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HDAC2=-0.107；PANK1=0.313，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HDAC2=-0.063；PANK1=0.274。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HDAC2（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt Q92769；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/Q92769/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 9346952。
- 作者control-scaled single -0.063，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0122；HCT116 0.129；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.5；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.637。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.65/3.59；H3K27ac 2.4/2.43；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.292，FDR 0.0183；SW480高耐药对亲本：log2FC -0.298，FDR 0.000211。
- GSE76092 log2FC -0.303，95%均值对比区间 [-0.446, -0.161]，BH FDR 0.00484，BY FDR 0.051；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 10.8，强度分位数 0.781；耐药株均值 10.5，强度分位数 0.752。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.123 / FDR 0.412；GSE289766 HCT116 log2FC -0.874；GSE42387 HCT116 log2FC -0.185、HT29 log2FC 0.067。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 67 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes](https://reactome.org/content/detail/R-HSA-9615017)：CAMERA Direction=Down，FDR=0.564，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Death Receptor Signaling](https://reactome.org/content/detail/R-HSA-73887)：CAMERA Direction=Down，FDR=0.564，nGenes=155；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### PANK1（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q8TE04；EC 2.7.1.33；[来源1](https://www.uniprot.org/uniprotkb/Q8TE04/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.274，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.0626；HCT116 -0.0532；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.92；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.128。
- HCT116/OxR作者模型通量：映射反应数 5；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.0894、OxR 0.188。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.52/3.23；H3K27ac 1.79/2.02；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.0719，FDR 0.746；SW480高耐药对亲本：log2FC -0.567，FDR 0.000238。
- GSE76092 log2FC 0.403，95%均值对比区间 [0.252, 0.554]，BH FDR 0.00125，BY FDR 0.0132；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 10.3，强度分位数 0.736；耐药株均值 10.7，强度分位数 0.774。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.229 / FDR 0.335；GSE289766 HCT116 log2FC 0.452；GSE42387 HCT116 log2FC -0.0718、HT29 log2FC 0.191。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC -3.11；有效 guide 数 2；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Metabolism of vitamins and cofactors](https://reactome.org/content/detail/R-HSA-196854)：CAMERA Direction=Down，FDR=0.968，nGenes=195；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of water-soluble vitamins and cofactors](https://reactome.org/content/detail/R-HSA-196849)：CAMERA Direction=Up，FDR=0.975，nGenes=124；ORA上调集合FDR=1、下调集合FDR=1。
- [Vitamin B5 (pantothenate) metabolism](https://reactome.org/content/detail/R-HSA-199220)：CAMERA Direction=Up，FDR=0.99，nGenes=20；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 03. NAMPT + HDAC2

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 122；多样化 Top30 排名 3。Priority8：3。

## 计算预测与已测单基因

- 预测加性 GI：**-0.525**；上游 90% GI 预测区间 [-0.719, -0.332]；跨随机种子 SD 0.0178。
- 预测联合效应 -0.449 = 已测单基因和 0.0768 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HDAC2=-0.107；NAMPT=0.184，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HDAC2=-0.063；NAMPT=0.19。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HDAC2（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt Q92769；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/Q92769/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 9346952。
- 作者control-scaled single -0.063，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0122；HCT116 0.129；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.5；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.637。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.65/3.59；H3K27ac 2.4/2.43；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.292，FDR 0.0183；SW480高耐药对亲本：log2FC -0.298，FDR 0.000211。
- GSE76092 log2FC -0.303，95%均值对比区间 [-0.446, -0.161]，BH FDR 0.00484，BY FDR 0.051；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 10.8，强度分位数 0.781；耐药株均值 10.5，强度分位数 0.752。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.123 / FDR 0.412；GSE289766 HCT116 log2FC -0.874；GSE42387 HCT116 log2FC -0.185、HT29 log2FC 0.067。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 67 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes](https://reactome.org/content/detail/R-HSA-9615017)：CAMERA Direction=Down，FDR=0.564，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Death Receptor Signaling](https://reactome.org/content/detail/R-HSA-73887)：CAMERA Direction=Down，FDR=0.564，nGenes=155；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### NAMPT（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P43490；EC 2.4.2.12；[来源1](https://www.uniprot.org/uniprotkb/P43490/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.19，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.283；HCT116 -0.316；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.627。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.527、OxR 0.706。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.72/3.53；H3K27ac 2.27/2.58；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.97，FDR 6.24e-14；SW480高耐药对亲本：log2FC -0.276，FDR 3.75e-05。
- GSE76092 log2FC 0.433，95%均值对比区间 [0.245, 0.622]，BH FDR 0.0032，BY FDR 0.0338；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 9.88，强度分位数 0.699；耐药株均值 10.3，强度分位数 0.738。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0696 / FDR 0.707；GSE289766 HCT116 log2FC 1.19；GSE42387 HCT116 log2FC -0.129、HT29 log2FC 0.898。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 9 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Transcriptional Regulation by NPAS4](https://reactome.org/content/detail/R-HSA-9634815)：CAMERA Direction=Down，FDR=0.761，nGenes=34；ORA上调集合FDR=1、下调集合FDR=1。
- [NAD+ and Sirtuins regulate the circadian clock](https://reactome.org/content/detail/R-HSA-9958957)：CAMERA Direction=Down，FDR=0.852，nGenes=43；ORA上调集合FDR=1、下调集合FDR=1。
- [NPAS4 regulates expression of target genes](https://reactome.org/content/detail/R-HSA-9768919)：CAMERA Direction=Down，FDR=0.954，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 04. NAMPT + KAT6A

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 142；多样化 Top30 排名 4。Priority8：4。

## 计算预测与已测单基因

- 预测加性 GI：**-0.516**；上游 90% GI 预测区间 [-0.71, -0.323]；跨随机种子 SD 0.0215。
- 预测联合效应 -0.416 = 已测单基因和 0.101 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：KAT6A=-0.0832；NAMPT=0.184，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：KAT6A=-0.0921；NAMPT=0.19。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### KAT6A（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt Q92794；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q92794/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 11313971。
- 作者control-scaled single -0.0921，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.115；HCT116 0.174；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.215。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.79/3.48；H3K27ac 2.37/2.1；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.161，FDR 0.179；SW480高耐药对亲本：log2FC 0.304，FDR 0.000196。
- GSE76092 log2FC 0.23，95%均值对比区间 [0.0823, 0.378]，BH FDR 0.0244，BY FDR 0.257；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 10.5，强度分位数 0.753；耐药株均值 10.7，强度分位数 0.776。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.0936 / FDR 0.588；GSE289766 HCT116 log2FC 0.891；GSE42387 HCT116 log2FC 0.131、HT29 log2FC 0.178。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Regulation of TP53 Activity](https://reactome.org/content/detail/R-HSA-5633007)：CAMERA Direction=Up，FDR=0.784，nGenes=161；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of TP53 Activity through Acetylation](https://reactome.org/content/detail/R-HSA-6804758)：CAMERA Direction=Down，FDR=0.954，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### NAMPT（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P43490；EC 2.4.2.12；[来源1](https://www.uniprot.org/uniprotkb/P43490/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.19，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.283；HCT116 -0.316；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.627。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.527、OxR 0.706。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.72/3.53；H3K27ac 2.27/2.58；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.97，FDR 6.24e-14；SW480高耐药对亲本：log2FC -0.276，FDR 3.75e-05。
- GSE76092 log2FC 0.433，95%均值对比区间 [0.245, 0.622]，BH FDR 0.0032，BY FDR 0.0338；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 9.88，强度分位数 0.699；耐药株均值 10.3，强度分位数 0.738。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0696 / FDR 0.707；GSE289766 HCT116 log2FC 1.19；GSE42387 HCT116 log2FC -0.129、HT29 log2FC 0.898。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 9 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Transcriptional Regulation by NPAS4](https://reactome.org/content/detail/R-HSA-9634815)：CAMERA Direction=Down，FDR=0.761，nGenes=34；ORA上调集合FDR=1、下调集合FDR=1。
- [NAD+ and Sirtuins regulate the circadian clock](https://reactome.org/content/detail/R-HSA-9958957)：CAMERA Direction=Down，FDR=0.852，nGenes=43；ORA上调集合FDR=1、下调集合FDR=1。
- [NPAS4 regulates expression of target genes](https://reactome.org/content/detail/R-HSA-9768919)：CAMERA Direction=Down，FDR=0.954，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 05. NAMPT + HDAC1

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 146；多样化 Top30 排名 5。Priority8：5。

## 计算预测与已测单基因

- 预测加性 GI：**-0.516**；上游 90% GI 预测区间 [-0.709, -0.322]；跨随机种子 SD 0.0231。
- 预测联合效应 -0.168 = 已测单基因和 0.347 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HDAC1=0.163；NAMPT=0.184，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HDAC1=0.0883；NAMPT=0.19。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HDAC1（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt Q13547；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/Q13547/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 10220385。
- 作者control-scaled single 0.0883，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.00838；HCT116 -0.151；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 12.4；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.634。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.44/3.28；H3K27ac 1.61/1.42；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.162，FDR 0.105；SW480高耐药对亲本：log2FC -0.153，FDR 0.0474。
- GSE76092 log2FC -0.35，95%均值对比区间 [-0.442, -0.258]，BH FDR 0.000133，BY FDR 0.0014；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 14.9，强度分位数 0.977；耐药株均值 14.6，强度分位数 0.97。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.105 / FDR 0.511；GSE289766 HCT116 log2FC -0.866；GSE42387 HCT116 log2FC -0.212、HT29 log2FC -0.412。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 100 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes](https://reactome.org/content/detail/R-HSA-9615017)：CAMERA Direction=Down，FDR=0.564，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Estrogen-dependent gene expression](https://reactome.org/content/detail/R-HSA-9018519)：CAMERA Direction=Down，FDR=0.564，nGenes=116；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### NAMPT（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P43490；EC 2.4.2.12；[来源1](https://www.uniprot.org/uniprotkb/P43490/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.19，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.283；HCT116 -0.316；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.627。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.527、OxR 0.706。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.72/3.53；H3K27ac 2.27/2.58；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.97，FDR 6.24e-14；SW480高耐药对亲本：log2FC -0.276，FDR 3.75e-05。
- GSE76092 log2FC 0.433，95%均值对比区间 [0.245, 0.622]，BH FDR 0.0032，BY FDR 0.0338；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 9.88，强度分位数 0.699；耐药株均值 10.3，强度分位数 0.738。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0696 / FDR 0.707；GSE289766 HCT116 log2FC 1.19；GSE42387 HCT116 log2FC -0.129、HT29 log2FC 0.898。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 9 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Transcriptional Regulation by NPAS4](https://reactome.org/content/detail/R-HSA-9634815)：CAMERA Direction=Down，FDR=0.761，nGenes=34；ORA上调集合FDR=1、下调集合FDR=1。
- [NAD+ and Sirtuins regulate the circadian clock](https://reactome.org/content/detail/R-HSA-9958957)：CAMERA Direction=Down，FDR=0.852，nGenes=43；ORA上调集合FDR=1、下调集合FDR=1。
- [NPAS4 regulates expression of target genes](https://reactome.org/content/detail/R-HSA-9768919)：CAMERA Direction=Down，FDR=0.954，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 06. PANK1 + KAT6A

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 181；多样化 Top30 排名 6。Priority8：6。

## 计算预测与已测单基因

- 预测加性 GI：**-0.5**；上游 90% GI 预测区间 [-0.694, -0.307]；跨随机种子 SD 0.013。
- 预测联合效应 -0.271 = 已测单基因和 0.229 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：KAT6A=-0.0832；PANK1=0.313，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：KAT6A=-0.0921；PANK1=0.274。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### KAT6A（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt Q92794；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q92794/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 11313971。
- 作者control-scaled single -0.0921，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.115；HCT116 0.174；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.215。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.79/3.48；H3K27ac 2.37/2.1；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.161，FDR 0.179；SW480高耐药对亲本：log2FC 0.304，FDR 0.000196。
- GSE76092 log2FC 0.23，95%均值对比区间 [0.0823, 0.378]，BH FDR 0.0244，BY FDR 0.257；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 10.5，强度分位数 0.753；耐药株均值 10.7，强度分位数 0.776。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.0936 / FDR 0.588；GSE289766 HCT116 log2FC 0.891；GSE42387 HCT116 log2FC 0.131、HT29 log2FC 0.178。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Regulation of TP53 Activity](https://reactome.org/content/detail/R-HSA-5633007)：CAMERA Direction=Up，FDR=0.784，nGenes=161；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of TP53 Activity through Acetylation](https://reactome.org/content/detail/R-HSA-6804758)：CAMERA Direction=Down，FDR=0.954，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### PANK1（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q8TE04；EC 2.7.1.33；[来源1](https://www.uniprot.org/uniprotkb/Q8TE04/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.274，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.0626；HCT116 -0.0532；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.92；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.128。
- HCT116/OxR作者模型通量：映射反应数 5；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.0894、OxR 0.188。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.52/3.23；H3K27ac 1.79/2.02；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.0719，FDR 0.746；SW480高耐药对亲本：log2FC -0.567，FDR 0.000238。
- GSE76092 log2FC 0.403，95%均值对比区间 [0.252, 0.554]，BH FDR 0.00125，BY FDR 0.0132；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 10.3，强度分位数 0.736；耐药株均值 10.7，强度分位数 0.774。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.229 / FDR 0.335；GSE289766 HCT116 log2FC 0.452；GSE42387 HCT116 log2FC -0.0718、HT29 log2FC 0.191。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC -3.11；有效 guide 数 2；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Metabolism of vitamins and cofactors](https://reactome.org/content/detail/R-HSA-196854)：CAMERA Direction=Down，FDR=0.968，nGenes=195；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of water-soluble vitamins and cofactors](https://reactome.org/content/detail/R-HSA-196849)：CAMERA Direction=Up，FDR=0.975，nGenes=124；ORA上调集合FDR=1、下调集合FDR=1。
- [Vitamin B5 (pantothenate) metabolism](https://reactome.org/content/detail/R-HSA-199220)：CAMERA Direction=Up，FDR=0.99，nGenes=20；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 07. IMPDH2 + HDAC2

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 221；多样化 Top30 排名 7。Priority8：7。

## 计算预测与已测单基因

- 预测加性 GI：**-0.484**；上游 90% GI 预测区间 [-0.677, -0.291]；跨随机种子 SD 0.0326。
- 预测联合效应 -0.653 = 已测单基因和 -0.169 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HDAC2=-0.107；IMPDH2=-0.0614，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HDAC2=-0.063；IMPDH2=-0.082。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HDAC2（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt Q92769；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/Q92769/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 9346952。
- 作者control-scaled single -0.063，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0122；HCT116 0.129；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.5；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.637。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.65/3.59；H3K27ac 2.4/2.43；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.292，FDR 0.0183；SW480高耐药对亲本：log2FC -0.298，FDR 0.000211。
- GSE76092 log2FC -0.303，95%均值对比区间 [-0.446, -0.161]，BH FDR 0.00484，BY FDR 0.051；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 10.8，强度分位数 0.781；耐药株均值 10.5，强度分位数 0.752。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.123 / FDR 0.412；GSE289766 HCT116 log2FC -0.874；GSE42387 HCT116 log2FC -0.185、HT29 log2FC 0.067。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 67 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes](https://reactome.org/content/detail/R-HSA-9615017)：CAMERA Direction=Down，FDR=0.564，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Death Receptor Signaling](https://reactome.org/content/detail/R-HSA-73887)：CAMERA Direction=Down，FDR=0.564，nGenes=155；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### IMPDH2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P12268；EC 1.1.1.205；[来源1](https://www.uniprot.org/uniprotkb/P12268/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single -0.082，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.471；HCT116 -0.885；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 12.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.687。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.311、OxR 0.481。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.92/3.04；H3K27ac 2.03/1.65；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.149，FDR 0.143；SW480高耐药对亲本：log2FC 0.315，FDR 2.47e-06。
- GSE76092 log2FC 0.187，95%均值对比区间 [0.0639, 0.309]，BH FDR 0.0272，BY FDR 0.287；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 13.5，强度分位数 0.941；耐药株均值 13.7，强度分位数 0.946。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0307 / FDR 0.869；GSE289766 HCT116 log2FC -1.28；GSE42387 HCT116 log2FC -0.173、HT29 log2FC -0.000991。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -0.971；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 0.667。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 8 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Azathioprine ADME](https://reactome.org/content/detail/R-HSA-9748787)：CAMERA Direction=Up，FDR=0.784，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。
- [Drug ADME](https://reactome.org/content/detail/R-HSA-9748784)：CAMERA Direction=Down，FDR=0.784，nGenes=97；ORA上调集合FDR=1、下调集合FDR=1。
- [Purine ribonucleoside monophosphate biosynthesis](https://reactome.org/content/detail/R-HSA-73817)：CAMERA Direction=Down，FDR=0.968，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

共有 2 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Potential therapeutics for SARS](https://reactome.org/content/detail/R-HSA-9679191)：CAMERA Direction=Up，FDR=0.985，nGenes=100；ORA上调集合FDR=1、下调集合FDR=1。
- [SARS-CoV Infections](https://reactome.org/content/detail/R-HSA-9679506)：CAMERA Direction=Down，FDR=0.99，nGenes=415；ORA上调集合FDR=1、下调集合FDR=1。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 08. HK2 + KDM3A

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 408；多样化 Top30 排名 8。Priority8：8。

## 计算预测与已测单基因

- 预测加性 GI：**-0.378**；上游 90% GI 预测区间 [-0.572, -0.185]；跨随机种子 SD 0.0254。
- 预测联合效应 -0.248 = 已测单基因和 0.13 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HK2=0.151；KDM3A=-0.0212，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HK2=0.14；KDM3A=0.00238。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HK2（严格代谢酶）

- 机制：Reviewed EC 2.7.1.1 and direct hexose phosphorylation reactions establish hexokinase catalysis in glycolysis.
- UniProt P52789；EC 2.7.1.1；[来源1](https://www.uniprot.org/uniprotkb/P52789/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.14，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.229；HCT116 -0.427；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.5；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.284。
- HCT116/OxR作者模型通量：映射反应数 7；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.472、OxR 0.497。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.96/3.6；H3K27ac 2.05/1.83；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.221，FDR 0.341；SW480高耐药对亲本：log2FC -0.43，FDR 2.86e-05。
- GSE76092 log2FC 0.177，95%均值对比区间 [0.0453, 0.308]，BH FDR 0.0462，BY FDR 0.486；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 9.23，强度分位数 0.638；耐药株均值 9.41，强度分位数 0.654。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.324 / FDR 0.0307；GSE289766 HCT116 log2FC -0.765；GSE42387 HCT116 log2FC -0.232、HT29 log2FC -0.378。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.01；有效 guide 数 4；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Glycolysis](https://reactome.org/content/detail/R-HSA-70171)：CAMERA Direction=Down，FDR=0.869，nGenes=72；ORA上调集合FDR=1、下调集合FDR=1。
- [Glucose metabolism](https://reactome.org/content/detail/R-HSA-70326)：CAMERA Direction=Down，FDR=0.911，nGenes=83；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of carbohydrates and carbohydrate derivatives](https://reactome.org/content/detail/R-HSA-71387)：CAMERA Direction=Down，FDR=0.968，nGenes=296；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KDM3A（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q9Y4C1；EC 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q9Y4C1/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 16603237。
- 作者control-scaled single 0.00238，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.123；HCT116 -0.228；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10.3；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.186。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.58/3.24；H3K27ac 1.92/1.96；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.0337，FDR 0.864；SW480高耐药对亲本：log2FC 0.482，FDR 4.62e-06。
- GSE76092 log2FC -0.353，95%均值对比区间 [-0.577, -0.13]，BH FDR 0.0227，BY FDR 0.239；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 7.87，强度分位数 0.51；耐药株均值 7.52，强度分位数 0.476。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0893 / FDR 0.738；GSE289766 HCT116 log2FC 1.58；GSE42387 HCT116 log2FC -0.22、HT29 log2FC -0.117。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [NR1H2 and NR1H3-mediated signaling](https://reactome.org/content/detail/R-HSA-9024446)：CAMERA Direction=Down，FDR=0.679，nGenes=46；ORA上调集合FDR=1、下调集合FDR=1。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Signaling by Nuclear Receptors](https://reactome.org/content/detail/R-HSA-9006931)：CAMERA Direction=Down，FDR=0.745，nGenes=260；ORA上调集合FDR=1、下调集合FDR=0.812。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 09. GALK1 + HDAC1

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 69；多样化 Top30 排名 9。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.553**；上游 90% GI 预测区间 [-0.747, -0.36]；跨随机种子 SD 0.0202。
- 预测联合效应 -0.201 = 已测单基因和 0.352 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：GALK1=0.189；HDAC1=0.163，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：GALK1=0.138；HDAC1=0.0883。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### GALK1（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P51570；EC 2.7.1.6；[来源1](https://www.uniprot.org/uniprotkb/P51570/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.138，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.165；HCT116 -0.0566；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.89；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.308。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.132、OxR 0.291。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.94/2.88；H3K27ac 1.57/1.33；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.0624，FDR 0.801；SW480高耐药对亲本：log2FC -0.303，FDR 0.0423。
- GSE76092 log2FC -0.0933，95%均值对比区间 [-0.234, 0.0473]，BH FDR 0.33，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 9.41，强度分位数 0.655；耐药株均值 9.32，强度分位数 0.645。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.103 / FDR 0.609；GSE289766 HCT116 log2FC -0.551；GSE42387 HCT116 log2FC -0.467、HT29 log2FC -0.0579。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -0.33；有效 guide 数 3；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 4 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Metabolism of carbohydrates and carbohydrate derivatives](https://reactome.org/content/detail/R-HSA-71387)：CAMERA Direction=Down，FDR=0.968，nGenes=296；ORA上调集合FDR=1、下调集合FDR=1。
- [Diseases of metabolism](https://reactome.org/content/detail/R-HSA-5668914)：CAMERA Direction=Up，FDR=0.986，nGenes=268；ORA上调集合FDR=1、下调集合FDR=1。
- [Diseases associated with glycosylation precursor biosynthesis](https://reactome.org/content/detail/R-HSA-5609975)：CAMERA Direction=Down，FDR=0.993，nGenes=16；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### HDAC1（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt Q13547；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/Q13547/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 10220385。
- 作者control-scaled single 0.0883，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.00838；HCT116 -0.151；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 12.4；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.634。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.44/3.28；H3K27ac 1.61/1.42；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.162，FDR 0.105；SW480高耐药对亲本：log2FC -0.153，FDR 0.0474。
- GSE76092 log2FC -0.35，95%均值对比区间 [-0.442, -0.258]，BH FDR 0.000133，BY FDR 0.0014；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 14.9，强度分位数 0.977；耐药株均值 14.6，强度分位数 0.97。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.105 / FDR 0.511；GSE289766 HCT116 log2FC -0.866；GSE42387 HCT116 log2FC -0.212、HT29 log2FC -0.412。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 100 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [FOXO-mediated transcription of oxidative stress, metabolic and neuronal genes](https://reactome.org/content/detail/R-HSA-9615017)：CAMERA Direction=Down，FDR=0.564，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Estrogen-dependent gene expression](https://reactome.org/content/detail/R-HSA-9018519)：CAMERA Direction=Down，FDR=0.564，nGenes=116；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 10. UCK2 + KAT6A

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 92；多样化 Top30 排名 10。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.541**；上游 90% GI 预测区间 [-0.735, -0.348]；跨随机种子 SD 0.0321。
- 预测联合效应 -0.396 = 已测单基因和 0.145 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：KAT6A=-0.0832；UCK2=0.228，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：KAT6A=-0.0921；UCK2=0.211。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### KAT6A（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt Q92794；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q92794/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 11313971。
- 作者control-scaled single -0.0921，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.115；HCT116 0.174；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.215。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.79/3.48；H3K27ac 2.37/2.1；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.161，FDR 0.179；SW480高耐药对亲本：log2FC 0.304，FDR 0.000196。
- GSE76092 log2FC 0.23，95%均值对比区间 [0.0823, 0.378]，BH FDR 0.0244，BY FDR 0.257；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 10.5，强度分位数 0.753；耐药株均值 10.7，强度分位数 0.776。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.0936 / FDR 0.588；GSE289766 HCT116 log2FC 0.891；GSE42387 HCT116 log2FC 0.131、HT29 log2FC 0.178。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Regulation of TP53 Activity](https://reactome.org/content/detail/R-HSA-5633007)：CAMERA Direction=Up，FDR=0.784，nGenes=161；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of TP53 Activity through Acetylation](https://reactome.org/content/detail/R-HSA-6804758)：CAMERA Direction=Down，FDR=0.954，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### UCK2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q9BZX2；EC 2.7.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q9BZX2/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.211，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.00677；HCT116 -0.261；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.25；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.274。
- HCT116/OxR作者模型通量：映射反应数 6；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.266、OxR 0.548。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.68/3.83；H3K27ac 2.29/2.46；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.258，FDR 0.044；SW480高耐药对亲本：log2FC -0.239，FDR 0.00416。
- GSE76092 log2FC -0.04，95%均值对比区间 [-0.183, 0.103]，BH FDR 0.707，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 13，强度分位数 0.919；耐药株均值 12.9，强度分位数 0.917。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.0624 / FDR 0.76；GSE289766 HCT116 log2FC -1.88；GSE42387 HCT116 log2FC -0.00537、HT29 log2FC 0.437。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.25；有效 guide 数 2；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。
- [Pyrimidine salvage](https://reactome.org/content/detail/R-HSA-73614)：CAMERA Direction=Up，FDR=0.843，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 11. HK2 + KAT7

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 93；多样化 Top30 排名 11。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.541**；上游 90% GI 预测区间 [-0.735, -0.348]；跨随机种子 SD 0.0159。
- 预测联合效应 -0.261 = 已测单基因和 0.281 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HK2=0.151；KAT7=0.13，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HK2=0.14；KAT7=0.122。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HK2（严格代谢酶）

- 机制：Reviewed EC 2.7.1.1 and direct hexose phosphorylation reactions establish hexokinase catalysis in glycolysis.
- UniProt P52789；EC 2.7.1.1；[来源1](https://www.uniprot.org/uniprotkb/P52789/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.14，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.229；HCT116 -0.427；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.5；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.284。
- HCT116/OxR作者模型通量：映射反应数 7；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.472、OxR 0.497。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.96/3.6；H3K27ac 2.05/1.83；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.221，FDR 0.341；SW480高耐药对亲本：log2FC -0.43，FDR 2.86e-05。
- GSE76092 log2FC 0.177，95%均值对比区间 [0.0453, 0.308]，BH FDR 0.0462，BY FDR 0.486；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 9.23，强度分位数 0.638；耐药株均值 9.41，强度分位数 0.654。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.324 / FDR 0.0307；GSE289766 HCT116 log2FC -0.765；GSE42387 HCT116 log2FC -0.232、HT29 log2FC -0.378。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.01；有效 guide 数 4；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Glycolysis](https://reactome.org/content/detail/R-HSA-70171)：CAMERA Direction=Down，FDR=0.869，nGenes=72；ORA上调集合FDR=1、下调集合FDR=1。
- [Glucose metabolism](https://reactome.org/content/detail/R-HSA-70326)：CAMERA Direction=Down，FDR=0.911，nGenes=83；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of carbohydrates and carbohydrate derivatives](https://reactome.org/content/detail/R-HSA-71387)：CAMERA Direction=Down，FDR=0.968，nGenes=296；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KAT7（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt O95251；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/O95251/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 10438470。
- 作者control-scaled single 0.122，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.232；HCT116 -0.459；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.15。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.82/3.07；H3K27ac 2.06/2.02；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.146，FDR 0.195；SW480高耐药对亲本：log2FC 0.385，FDR 5.81e-08。
- GSE76092 log2FC 0.147，95%均值对比区间 [-0.0427, 0.336]，BH FDR 0.248，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 9.09，强度分位数 0.625；耐药株均值 9.23，强度分位数 0.638。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0764 / FDR 0.676；GSE289766 HCT116 log2FC 0.167；GSE42387 HCT116 log2FC -0.0478、HT29 log2FC 0.0471。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 12. HK2 + KAT6B

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 110；多样化 Top30 排名 12。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.53**；上游 90% GI 预测区间 [-0.723, -0.337]；跨随机种子 SD 0.0238。
- 预测联合效应 -0.382 = 已测单基因和 0.148 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HK2=0.151；KAT6B=-0.00315，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HK2=0.14；KAT6B=0.0533。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HK2（严格代谢酶）

- 机制：Reviewed EC 2.7.1.1 and direct hexose phosphorylation reactions establish hexokinase catalysis in glycolysis.
- UniProt P52789；EC 2.7.1.1；[来源1](https://www.uniprot.org/uniprotkb/P52789/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.14，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.229；HCT116 -0.427；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.5；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.284。
- HCT116/OxR作者模型通量：映射反应数 7；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.472、OxR 0.497。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.96/3.6；H3K27ac 2.05/1.83；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.221，FDR 0.341；SW480高耐药对亲本：log2FC -0.43，FDR 2.86e-05。
- GSE76092 log2FC 0.177，95%均值对比区间 [0.0453, 0.308]，BH FDR 0.0462，BY FDR 0.486；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 9.23，强度分位数 0.638；耐药株均值 9.41，强度分位数 0.654。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.324 / FDR 0.0307；GSE289766 HCT116 log2FC -0.765；GSE42387 HCT116 log2FC -0.232、HT29 log2FC -0.378。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.01；有效 guide 数 4；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Glycolysis](https://reactome.org/content/detail/R-HSA-70171)：CAMERA Direction=Down，FDR=0.869，nGenes=72；ORA上调集合FDR=1、下调集合FDR=1。
- [Glucose metabolism](https://reactome.org/content/detail/R-HSA-70326)：CAMERA Direction=Down，FDR=0.911，nGenes=83；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of carbohydrates and carbohydrate derivatives](https://reactome.org/content/detail/R-HSA-71387)：CAMERA Direction=Down，FDR=0.968，nGenes=296；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KAT6B（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt Q8WYB5；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q8WYB5/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 10497217。
- 作者control-scaled single 0.0533，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.121；HCT116 0.159；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.168。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.52/3.19；H3K27ac 2.13/2.12；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.771，FDR 7.87e-09；SW480高耐药对亲本：log2FC 0.0843，FDR 0.515。
- GSE76092 log2FC 0，95%均值对比区间 [-0.147, 0.147]，BH FDR 1，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 7.46，强度分位数 0.472；耐药株均值 7.46，强度分位数 0.47。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0249 / FDR 0.923；GSE289766 HCT116 log2FC 0.688；GSE42387 HCT116 log2FC 0.138、HT29 log2FC -0.312。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 13. ADK + KAT7

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 175；多样化 Top30 排名 13。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.505**；上游 90% GI 预测区间 [-0.698, -0.311]；跨随机种子 SD 0.0285。
- 预测联合效应 -0.322 = 已测单基因和 0.182 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：ADK=0.0528；KAT7=0.13，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：ADK=0.0762；KAT7=0.122。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### ADK（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P55263；EC 2.7.1.-; 2.7.1.20；[来源1](https://www.uniprot.org/uniprotkb/P55263/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.0762，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0577；HCT116 -0.0966；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.79；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.39。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.554、OxR 0.856。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.34/3.49；H3K27ac 2.23/1.94；H3K27me3 1.12/1.1。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.421，FDR 8.16e-05；SW480高耐药对亲本：log2FC -0.241，FDR 0.0128。
- GSE76092 log2FC 0.28，95%均值对比区间 [0.141, 0.419]，BH FDR 0.00644，BY FDR 0.0678；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 11.4，强度分位数 0.832；耐药株均值 11.7，强度分位数 0.85。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.321 / FDR 0.0135；GSE289766 HCT116 log2FC -0.623；GSE42387 HCT116 log2FC 0.13、HT29 log2FC 0.674。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -1.46；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 5 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Ribavirin ADME](https://reactome.org/content/detail/R-HSA-9755088)：CAMERA Direction=Up，FDR=0.564，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Purine salvage](https://reactome.org/content/detail/R-HSA-74217)：CAMERA Direction=Up，FDR=0.564，nGenes=13；ORA上调集合FDR=1、下调集合FDR=1。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KAT7（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt O95251；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/O95251/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 10438470。
- 作者control-scaled single 0.122，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.232；HCT116 -0.459；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.15。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.82/3.07；H3K27ac 2.06/2.02；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.146，FDR 0.195；SW480高耐药对亲本：log2FC 0.385，FDR 5.81e-08。
- GSE76092 log2FC 0.147，95%均值对比区间 [-0.0427, 0.336]，BH FDR 0.248，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 9.09，强度分位数 0.625；耐药株均值 9.23，强度分位数 0.638。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0764 / FDR 0.676；GSE289766 HCT116 log2FC 0.167；GSE42387 HCT116 log2FC -0.0478、HT29 log2FC 0.0471。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 14. IMPDH2 + KAT7

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 196；多样化 Top30 排名 14。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.494**；上游 90% GI 预测区间 [-0.688, -0.301]；跨随机种子 SD 0.024。
- 预测联合效应 -0.426 = 已测单基因和 0.0682 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：IMPDH2=-0.0614；KAT7=0.13，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：IMPDH2=-0.082；KAT7=0.122。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### IMPDH2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P12268；EC 1.1.1.205；[来源1](https://www.uniprot.org/uniprotkb/P12268/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single -0.082，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.471；HCT116 -0.885；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 12.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.687。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.311、OxR 0.481。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.92/3.04；H3K27ac 2.03/1.65；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.149，FDR 0.143；SW480高耐药对亲本：log2FC 0.315，FDR 2.47e-06。
- GSE76092 log2FC 0.187，95%均值对比区间 [0.0639, 0.309]，BH FDR 0.0272，BY FDR 0.287；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 13.5，强度分位数 0.941；耐药株均值 13.7，强度分位数 0.946。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0307 / FDR 0.869；GSE289766 HCT116 log2FC -1.28；GSE42387 HCT116 log2FC -0.173、HT29 log2FC -0.000991。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -0.971；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 0.667。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 8 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Azathioprine ADME](https://reactome.org/content/detail/R-HSA-9748787)：CAMERA Direction=Up，FDR=0.784，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。
- [Drug ADME](https://reactome.org/content/detail/R-HSA-9748784)：CAMERA Direction=Down，FDR=0.784，nGenes=97；ORA上调集合FDR=1、下调集合FDR=1。
- [Purine ribonucleoside monophosphate biosynthesis](https://reactome.org/content/detail/R-HSA-73817)：CAMERA Direction=Down，FDR=0.968，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KAT7（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt O95251；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/O95251/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 10438470。
- 作者control-scaled single 0.122，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.232；HCT116 -0.459；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.15。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.82/3.07；H3K27ac 2.06/2.02；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.146，FDR 0.195；SW480高耐药对亲本：log2FC 0.385，FDR 5.81e-08。
- GSE76092 log2FC 0.147，95%均值对比区间 [-0.0427, 0.336]，BH FDR 0.248，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 9.09，强度分位数 0.625；耐药株均值 9.23，强度分位数 0.638。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0764 / FDR 0.676；GSE289766 HCT116 log2FC 0.167；GSE42387 HCT116 log2FC -0.0478、HT29 log2FC 0.0471。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 15. IMPDH2 + KAT6B

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 224；多样化 Top30 排名 15。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.483**；上游 90% GI 预测区间 [-0.676, -0.289]；跨随机种子 SD 0.0262。
- 预测联合效应 -0.547 = 已测单基因和 -0.0646 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：IMPDH2=-0.0614；KAT6B=-0.00315，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：IMPDH2=-0.082；KAT6B=0.0533。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### IMPDH2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P12268；EC 1.1.1.205；[来源1](https://www.uniprot.org/uniprotkb/P12268/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single -0.082，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.471；HCT116 -0.885；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 12.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.687。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.311、OxR 0.481。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.92/3.04；H3K27ac 2.03/1.65；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.149，FDR 0.143；SW480高耐药对亲本：log2FC 0.315，FDR 2.47e-06。
- GSE76092 log2FC 0.187，95%均值对比区间 [0.0639, 0.309]，BH FDR 0.0272，BY FDR 0.287；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 13.5，强度分位数 0.941；耐药株均值 13.7，强度分位数 0.946。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0307 / FDR 0.869；GSE289766 HCT116 log2FC -1.28；GSE42387 HCT116 log2FC -0.173、HT29 log2FC -0.000991。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -0.971；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 0.667。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 8 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Azathioprine ADME](https://reactome.org/content/detail/R-HSA-9748787)：CAMERA Direction=Up，FDR=0.784，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。
- [Drug ADME](https://reactome.org/content/detail/R-HSA-9748784)：CAMERA Direction=Down，FDR=0.784，nGenes=97；ORA上调集合FDR=1、下调集合FDR=1。
- [Purine ribonucleoside monophosphate biosynthesis](https://reactome.org/content/detail/R-HSA-73817)：CAMERA Direction=Down，FDR=0.968，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KAT6B（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt Q8WYB5；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q8WYB5/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 10497217。
- 作者control-scaled single 0.0533，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.121；HCT116 0.159；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.168。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.52/3.19；H3K27ac 2.13/2.12；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.771，FDR 7.87e-09；SW480高耐药对亲本：log2FC 0.0843，FDR 0.515。
- GSE76092 log2FC 0，95%均值对比区间 [-0.147, 0.147]，BH FDR 1，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 7.46，强度分位数 0.472；耐药株均值 7.46，强度分位数 0.47。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0249 / FDR 0.923；GSE289766 HCT116 log2FC 0.688；GSE42387 HCT116 log2FC 0.138、HT29 log2FC -0.312。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 16. ADK + KDM3B

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 364；多样化 Top30 排名 16。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.411**；上游 90% GI 预测区间 [-0.605, -0.218]；跨随机种子 SD 0.0161。
- 预测联合效应 -0.3 = 已测单基因和 0.111 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：ADK=0.0528；KDM3B=0.058，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：ADK=0.0762；KDM3B=0.0541。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### ADK（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P55263；EC 2.7.1.-; 2.7.1.20；[来源1](https://www.uniprot.org/uniprotkb/P55263/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.0762，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0577；HCT116 -0.0966；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.79；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.39。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.554、OxR 0.856。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.34/3.49；H3K27ac 2.23/1.94；H3K27me3 1.12/1.1。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.421，FDR 8.16e-05；SW480高耐药对亲本：log2FC -0.241，FDR 0.0128。
- GSE76092 log2FC 0.28，95%均值对比区间 [0.141, 0.419]，BH FDR 0.00644，BY FDR 0.0678；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 11.4，强度分位数 0.832；耐药株均值 11.7，强度分位数 0.85。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.321 / FDR 0.0135；GSE289766 HCT116 log2FC -0.623；GSE42387 HCT116 log2FC 0.13、HT29 log2FC 0.674。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -1.46；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 5 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Ribavirin ADME](https://reactome.org/content/detail/R-HSA-9755088)：CAMERA Direction=Up，FDR=0.564，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Purine salvage](https://reactome.org/content/detail/R-HSA-74217)：CAMERA Direction=Up，FDR=0.564，nGenes=13；ORA上调集合FDR=1、下调集合FDR=1。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KDM3B（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q7LBC6；EC 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q7LBC6/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 16603237。
- 作者control-scaled single 0.0541，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.0286；HCT116 -0.0249；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.199。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.92/3.55；H3K27ac 2.8/2.6；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.107，FDR 0.334；SW480高耐药对亲本：log2FC 0.0938，FDR 0.215。
- GSE76092 log2FC -0.0467，95%均值对比区间 [-0.207, 0.114]，BH FDR 0.696，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 11.3，强度分位数 0.825；耐药株均值 11.3，强度分位数 0.822。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.00687 / FDR 0.974；GSE289766 HCT116 log2FC 0.0336；GSE42387 HCT116 log2FC -0.172、HT29 log2FC -0.518。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 17. ADK + PHF8

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 387；多样化 Top30 排名 17。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.397**；上游 90% GI 预测区间 [-0.59, -0.203]；跨随机种子 SD 0.0188。
- 预测联合效应 -0.451 = 已测单基因和 -0.0548 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：ADK=0.0528；PHF8=-0.108，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：ADK=0.0762；PHF8=-0.104。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### ADK（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P55263；EC 2.7.1.-; 2.7.1.20；[来源1](https://www.uniprot.org/uniprotkb/P55263/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.0762，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0577；HCT116 -0.0966；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.79；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.39。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.554、OxR 0.856。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.34/3.49；H3K27ac 2.23/1.94；H3K27me3 1.12/1.1。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.421，FDR 8.16e-05；SW480高耐药对亲本：log2FC -0.241，FDR 0.0128。
- GSE76092 log2FC 0.28，95%均值对比区间 [0.141, 0.419]，BH FDR 0.00644，BY FDR 0.0678；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 11.4，强度分位数 0.832；耐药株均值 11.7，强度分位数 0.85。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.321 / FDR 0.0135；GSE289766 HCT116 log2FC -0.623；GSE42387 HCT116 log2FC 0.13、HT29 log2FC 0.674。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -1.46；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 5 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Ribavirin ADME](https://reactome.org/content/detail/R-HSA-9755088)：CAMERA Direction=Up，FDR=0.564，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Purine salvage](https://reactome.org/content/detail/R-HSA-74217)：CAMERA Direction=Up，FDR=0.564，nGenes=13；ORA上调集合FDR=1、下调集合FDR=1。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### PHF8（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q9UPP1；EC 1.14.11.-; 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q9UPP1/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 21423274。
- 作者control-scaled single -0.104，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.318；HCT116 0.16；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.99；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.145。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.8/2.87；H3K27ac 1.74/1.97；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.289，FDR 0.0436；SW480高耐药对亲本：log2FC -0.65，FDR 1.27e-14。
- GSE76092 log2FC -0.103，95%均值对比区间 [-0.26, 0.0531]，BH FDR 0.332，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 8.57，强度分位数 0.573；耐药株均值 8.46，强度分位数 0.565。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.156 / FDR 0.351；GSE289766 HCT116 log2FC 0.576；GSE42387 HCT116 log2FC 0.375、HT29 log2FC -0.267。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Mitotic Prophase](https://reactome.org/content/detail/R-HSA-68875)：CAMERA Direction=Down，FDR=0.777，nGenes=108；ORA上调集合FDR=1、下调集合FDR=1。
- [M Phase](https://reactome.org/content/detail/R-HSA-68886)：CAMERA Direction=Up，FDR=0.801，nGenes=366；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 18. UCK2 + KAT6B

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 42；多样化 Top30 排名 18。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.572**；上游 90% GI 预测区间 [-0.765, -0.378]；跨随机种子 SD 0.0223。
- 预测联合效应 -0.346 = 已测单基因和 0.225 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：KAT6B=-0.00315；UCK2=0.228，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：KAT6B=0.0533；UCK2=0.211。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### KAT6B（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone acetyltransferase activity and defined HAT complexes establish a direct histone modifier; chromatin catalysis is a primary biochemical mechanism.
- UniProt Q8WYB5；EC 2.3.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q8WYB5/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 10497217。
- 作者control-scaled single 0.0533，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.121；HCT116 0.159；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.168。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.52/3.19；H3K27ac 2.13/2.12；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.771，FDR 7.87e-09；SW480高耐药对亲本：log2FC 0.0843，FDR 0.515。
- GSE76092 log2FC 0，95%均值对比区间 [-0.147, 0.147]，BH FDR 1，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 7.46，强度分位数 0.472；耐药株均值 7.46，强度分位数 0.47。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0249 / FDR 0.923；GSE289766 HCT116 log2FC 0.688；GSE42387 HCT116 log2FC 0.138、HT29 log2FC -0.312。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HATs acetylate histones](https://reactome.org/content/detail/R-HSA-3214847)：CAMERA Direction=Down，FDR=0.959，nGenes=106；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### UCK2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q9BZX2；EC 2.7.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q9BZX2/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.211，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.00677；HCT116 -0.261；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.25；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.274。
- HCT116/OxR作者模型通量：映射反应数 6；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.266、OxR 0.548。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.68/3.83；H3K27ac 2.29/2.46；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.258，FDR 0.044；SW480高耐药对亲本：log2FC -0.239，FDR 0.00416。
- GSE76092 log2FC -0.04，95%均值对比区间 [-0.183, 0.103]，BH FDR 0.707，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 13，强度分位数 0.919；耐药株均值 12.9，强度分位数 0.917。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.0624 / FDR 0.76；GSE289766 HCT116 log2FC -1.88；GSE42387 HCT116 log2FC -0.00537、HT29 log2FC 0.437。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.25；有效 guide 数 2；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。
- [Pyrimidine salvage](https://reactome.org/content/detail/R-HSA-73614)：CAMERA Direction=Up，FDR=0.843，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 19. UCK2 + KDM3B

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 331；多样化 Top30 排名 19。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.432**；上游 90% GI 预测区间 [-0.625, -0.238]；跨随机种子 SD 0.031。
- 预测联合效应 -0.145 = 已测单基因和 0.286 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：KDM3B=0.058；UCK2=0.228，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：KDM3B=0.0541；UCK2=0.211。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### KDM3B（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q7LBC6；EC 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q7LBC6/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 16603237。
- 作者control-scaled single 0.0541，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.0286；HCT116 -0.0249；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.199。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.92/3.55；H3K27ac 2.8/2.6；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.107，FDR 0.334；SW480高耐药对亲本：log2FC 0.0938，FDR 0.215。
- GSE76092 log2FC -0.0467，95%均值对比区间 [-0.207, 0.114]，BH FDR 0.696，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 11.3，强度分位数 0.825；耐药株均值 11.3，强度分位数 0.822。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.00687 / FDR 0.974；GSE289766 HCT116 log2FC 0.0336；GSE42387 HCT116 log2FC -0.172、HT29 log2FC -0.518。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### UCK2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q9BZX2；EC 2.7.1.48；[来源1](https://www.uniprot.org/uniprotkb/Q9BZX2/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.211，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.00677；HCT116 -0.261；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.25；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.274。
- HCT116/OxR作者模型通量：映射反应数 6；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.266、OxR 0.548。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.68/3.83；H3K27ac 2.29/2.46；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.258，FDR 0.044；SW480高耐药对亲本：log2FC -0.239，FDR 0.00416。
- GSE76092 log2FC -0.04，95%均值对比区间 [-0.183, 0.103]，BH FDR 0.707，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 13，强度分位数 0.919；耐药株均值 12.9，强度分位数 0.917。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.0624 / FDR 0.76；GSE289766 HCT116 log2FC -1.88；GSE42387 HCT116 log2FC -0.00537、HT29 log2FC 0.437。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.25；有效 guide 数 2；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。
- [Pyrimidine salvage](https://reactome.org/content/detail/R-HSA-73614)：CAMERA Direction=Up，FDR=0.843，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 20. DGUOK + KDM3B

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 333；多样化 Top30 排名 20。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.431**；上游 90% GI 预测区间 [-0.624, -0.238]；跨随机种子 SD 0.0128。
- 预测联合效应 -0.117 = 已测单基因和 0.314 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：DGUOK=0.256；KDM3B=0.058，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：DGUOK=0.21；KDM3B=0.0541。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### DGUOK（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q16854；EC 2.7.1.113; 2.7.1.76；[来源1](https://www.uniprot.org/uniprotkb/Q16854/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.21，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0552；HCT116 -0.202；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.74；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.668。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.267、OxR 0.549。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.48/2.79；H3K27ac 1.88/1.75；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.147，FDR 0.193；SW480高耐药对亲本：log2FC -0.271，FDR 0.00105。
- GSE76092 log2FC 0.0433，95%均值对比区间 [-0.0518, 0.139]，BH FDR 0.518，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 12.5，强度分位数 0.896；耐药株均值 12.5，强度分位数 0.897。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0711 / FDR 0.716；GSE289766 HCT116 log2FC -0.517；GSE42387 HCT116 log2FC -0.22、HT29 log2FC 0.132。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.27；有效 guide 数 4；guide 耗竭比例 0.75；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Purine salvage](https://reactome.org/content/detail/R-HSA-74217)：CAMERA Direction=Up，FDR=0.564，nGenes=13；ORA上调集合FDR=1、下调集合FDR=1。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KDM3B（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q7LBC6；EC 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q7LBC6/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 16603237。
- 作者control-scaled single 0.0541，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 0.0286；HCT116 -0.0249；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.2；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.199。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.92/3.55；H3K27ac 2.8/2.6；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.107，FDR 0.334；SW480高耐药对亲本：log2FC 0.0938，FDR 0.215。
- GSE76092 log2FC -0.0467，95%均值对比区间 [-0.207, 0.114]，BH FDR 0.696，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 11.3，强度分位数 0.825；耐药株均值 11.3，强度分位数 0.822。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.00687 / FDR 0.974；GSE289766 HCT116 log2FC 0.0336；GSE42387 HCT116 log2FC -0.172、HT29 log2FC -0.518。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：mixed_or_small；覆盖 3 项，方向支持比例 0。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin modifying enzymes](https://reactome.org/content/detail/R-HSA-3247509)：CAMERA Direction=Down，FDR=0.968，nGenes=236；ORA上调集合FDR=1、下调集合FDR=1。
- [Chromatin organization](https://reactome.org/content/detail/R-HSA-4839726)：CAMERA Direction=Down，FDR=0.968，nGenes=338；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 21. DGUOK + PHF8

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 358；多样化 Top30 排名 21。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.418**；上游 90% GI 预测区间 [-0.611, -0.224]；跨随机种子 SD 0.0253。
- 预测联合效应 -0.269 = 已测单基因和 0.148 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：DGUOK=0.256；PHF8=-0.108，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：DGUOK=0.21；PHF8=-0.104。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### DGUOK（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q16854；EC 2.7.1.113; 2.7.1.76；[来源1](https://www.uniprot.org/uniprotkb/Q16854/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.21，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0552；HCT116 -0.202；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.74；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.668。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.267、OxR 0.549。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.48/2.79；H3K27ac 1.88/1.75；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.147，FDR 0.193；SW480高耐药对亲本：log2FC -0.271，FDR 0.00105。
- GSE76092 log2FC 0.0433，95%均值对比区间 [-0.0518, 0.139]，BH FDR 0.518，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 12.5，强度分位数 0.896；耐药株均值 12.5，强度分位数 0.897。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0711 / FDR 0.716；GSE289766 HCT116 log2FC -0.517；GSE42387 HCT116 log2FC -0.22、HT29 log2FC 0.132。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.27；有效 guide 数 4；guide 耗竭比例 0.75；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Purine salvage](https://reactome.org/content/detail/R-HSA-74217)：CAMERA Direction=Up，FDR=0.564，nGenes=13；ORA上调集合FDR=1、下调集合FDR=1。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### PHF8（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q9UPP1；EC 1.14.11.-; 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q9UPP1/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 21423274。
- 作者control-scaled single -0.104，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.318；HCT116 0.16；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.99；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.145。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.8/2.87；H3K27ac 1.74/1.97；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.289，FDR 0.0436；SW480高耐药对亲本：log2FC -0.65，FDR 1.27e-14。
- GSE76092 log2FC -0.103，95%均值对比区间 [-0.26, 0.0531]，BH FDR 0.332，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 8.57，强度分位数 0.573；耐药株均值 8.46，强度分位数 0.565。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.156 / FDR 0.351；GSE289766 HCT116 log2FC 0.576；GSE42387 HCT116 log2FC 0.375、HT29 log2FC -0.267。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Mitotic Prophase](https://reactome.org/content/detail/R-HSA-68875)：CAMERA Direction=Down，FDR=0.777，nGenes=108；ORA上调集合FDR=1、下调集合FDR=1。
- [M Phase](https://reactome.org/content/detail/R-HSA-68886)：CAMERA Direction=Up，FDR=0.801，nGenes=366；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 22. NME2 + PHF8

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 418；多样化 Top30 排名 22。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.366**；上游 90% GI 预测区间 [-0.559, -0.172]；跨随机种子 SD 0.0181。
- 预测联合效应 -0.325 = 已测单基因和 0.0402 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：NME2=0.148；PHF8=-0.108，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：NME2=0.157；PHF8=-0.104。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### NME2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P22392；EC 2.7.13.3; 2.7.4.6；[来源1](https://www.uniprot.org/uniprotkb/P22392/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.157，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.361；HCT116 -0.437；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 13.3；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.252。
- HCT116/OxR作者模型通量：映射反应数 3e+01；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.633、OxR 1.07。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.86/3.69；H3K27ac 2.64/1.82；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.106，FDR 0.756；SW480高耐药对亲本：log2FC -0.594，FDR 0.299。
- GSE76092 log2FC -0.0133，95%均值对比区间 [-0.086, 0.0593]，BH FDR 0.818，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 17.6，强度分位数 0.998；耐药株均值 17.6，强度分位数 0.998。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.285 / FDR 0.494；GSE289766 HCT116 log2FC -1.76；GSE42387 HCT116 log2FC -0.322、HT29 log2FC 0.226。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -1.13；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 0.667。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Ribavirin ADME](https://reactome.org/content/detail/R-HSA-9755088)：CAMERA Direction=Up，FDR=0.564，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Azathioprine ADME](https://reactome.org/content/detail/R-HSA-9748787)：CAMERA Direction=Up，FDR=0.784，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。
- [Drug ADME](https://reactome.org/content/detail/R-HSA-9748784)：CAMERA Direction=Down，FDR=0.784，nGenes=97；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### PHF8（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q9UPP1；EC 1.14.11.-; 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q9UPP1/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 21423274。
- 作者control-scaled single -0.104，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.318；HCT116 0.16；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.99；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.145。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.8/2.87；H3K27ac 1.74/1.97；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.289，FDR 0.0436；SW480高耐药对亲本：log2FC -0.65，FDR 1.27e-14。
- GSE76092 log2FC -0.103，95%均值对比区间 [-0.26, 0.0531]，BH FDR 0.332，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 8.57，强度分位数 0.573；耐药株均值 8.46，强度分位数 0.565。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.156 / FDR 0.351；GSE289766 HCT116 log2FC 0.576；GSE42387 HCT116 log2FC 0.375、HT29 log2FC -0.267。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Mitotic Prophase](https://reactome.org/content/detail/R-HSA-68875)：CAMERA Direction=Down，FDR=0.777，nGenes=108；ORA上调集合FDR=1、下调集合FDR=1。
- [M Phase](https://reactome.org/content/detail/R-HSA-68886)：CAMERA Direction=Up，FDR=0.801，nGenes=366；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 23. GALK1 + KDM3A

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 374；多样化 Top30 排名 23。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.407**；上游 90% GI 预测区间 [-0.6, -0.214]；跨随机种子 SD 0.0794。
- 预测联合效应 -0.239 = 已测单基因和 0.168 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：GALK1=0.189；KDM3A=-0.0212，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：GALK1=0.138；KDM3A=0.00238。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 upper_quartile_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### GALK1（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P51570；EC 2.7.1.6；[来源1](https://www.uniprot.org/uniprotkb/P51570/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.138，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.165；HCT116 -0.0566；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.89；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.308。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.132、OxR 0.291。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.94/2.88；H3K27ac 1.57/1.33；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.0624，FDR 0.801；SW480高耐药对亲本：log2FC -0.303，FDR 0.0423。
- GSE76092 log2FC -0.0933，95%均值对比区间 [-0.234, 0.0473]，BH FDR 0.33，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 9.41，强度分位数 0.655；耐药株均值 9.32，强度分位数 0.645。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.103 / FDR 0.609；GSE289766 HCT116 log2FC -0.551；GSE42387 HCT116 log2FC -0.467、HT29 log2FC -0.0579。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -0.33；有效 guide 数 3；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 4 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Metabolism of carbohydrates and carbohydrate derivatives](https://reactome.org/content/detail/R-HSA-71387)：CAMERA Direction=Down，FDR=0.968，nGenes=296；ORA上调集合FDR=1、下调集合FDR=1。
- [Diseases of metabolism](https://reactome.org/content/detail/R-HSA-5668914)：CAMERA Direction=Up，FDR=0.986，nGenes=268；ORA上调集合FDR=1、下调集合FDR=1。
- [Diseases associated with glycosylation precursor biosynthesis](https://reactome.org/content/detail/R-HSA-5609975)：CAMERA Direction=Down，FDR=0.993，nGenes=16；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KDM3A（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q9Y4C1；EC 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q9Y4C1/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 16603237。
- 作者control-scaled single 0.00238，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.123；HCT116 -0.228；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10.3；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.186。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.58/3.24；H3K27ac 1.92/1.96；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.0337，FDR 0.864；SW480高耐药对亲本：log2FC 0.482，FDR 4.62e-06。
- GSE76092 log2FC -0.353，95%均值对比区间 [-0.577, -0.13]，BH FDR 0.0227，BY FDR 0.239；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 7.87，强度分位数 0.51；耐药株均值 7.52，强度分位数 0.476。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0893 / FDR 0.738；GSE289766 HCT116 log2FC 1.58；GSE42387 HCT116 log2FC -0.22、HT29 log2FC -0.117。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [NR1H2 and NR1H3-mediated signaling](https://reactome.org/content/detail/R-HSA-9024446)：CAMERA Direction=Down，FDR=0.679，nGenes=46；ORA上调集合FDR=1、下调集合FDR=1。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Signaling by Nuclear Receptors](https://reactome.org/content/detail/R-HSA-9006931)：CAMERA Direction=Down，FDR=0.745，nGenes=260；ORA上调集合FDR=1、下调集合FDR=0.812。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 24. DGUOK + KDM3A

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 459；多样化 Top30 排名 24。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.322**；上游 90% GI 预测区间 [-0.516, -0.129]；跨随机种子 SD 0.0492。
- 预测联合效应 -0.0876 = 已测单基因和 0.235 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：DGUOK=0.256；KDM3A=-0.0212，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：DGUOK=0.21；KDM3A=0.00238。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 upper_quartile_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### DGUOK（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q16854；EC 2.7.1.113; 2.7.1.76；[来源1](https://www.uniprot.org/uniprotkb/Q16854/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.21，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0552；HCT116 -0.202；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 9.74；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.668。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.267、OxR 0.549。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.48/2.79；H3K27ac 1.88/1.75；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.147，FDR 0.193；SW480高耐药对亲本：log2FC -0.271，FDR 0.00105。
- GSE76092 log2FC 0.0433，95%均值对比区间 [-0.0518, 0.139]，BH FDR 0.518，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 12.5，强度分位数 0.896；耐药株均值 12.5，强度分位数 0.897。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0711 / FDR 0.716；GSE289766 HCT116 log2FC -0.517；GSE42387 HCT116 log2FC -0.22、HT29 log2FC 0.132。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -2.27；有效 guide 数 4；guide 耗竭比例 0.75；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 3 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Purine salvage](https://reactome.org/content/detail/R-HSA-74217)：CAMERA Direction=Up，FDR=0.564，nGenes=13；ORA上调集合FDR=1、下调集合FDR=1。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### KDM3A（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed direct histone lysine demethylation establishes a core chromatin-modifying enzyme. The detailed substrate scope remains as stated in the reviewed record.
- UniProt Q9Y4C1；EC 1.14.11.65；[来源1](https://www.uniprot.org/uniprotkb/Q9Y4C1/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 16603237。
- 作者control-scaled single 0.00238，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.123；HCT116 -0.228；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10.3；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.186。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.58/3.24；H3K27ac 1.92/1.96；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.0337，FDR 0.864；SW480高耐药对亲本：log2FC 0.482，FDR 4.62e-06。
- GSE76092 log2FC -0.353，95%均值对比区间 [-0.577, -0.13]，BH FDR 0.0227，BY FDR 0.239；标签 resistant_down_BH_lt_005。
- 亲本相对芯片强度均值 7.87，强度分位数 0.51；耐药株均值 7.52，强度分位数 0.476。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.0893 / FDR 0.738；GSE289766 HCT116 log2FC 1.58；GSE42387 HCT116 log2FC -0.22、HT29 log2FC -0.117。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [NR1H2 and NR1H3-mediated signaling](https://reactome.org/content/detail/R-HSA-9024446)：CAMERA Direction=Down，FDR=0.679，nGenes=46；ORA上调集合FDR=1、下调集合FDR=1。
- [HDMs demethylate histones](https://reactome.org/content/detail/R-HSA-3214842)：CAMERA Direction=Down，FDR=0.682，nGenes=26；ORA上调集合FDR=1、下调集合FDR=1。
- [Signaling by Nuclear Receptors](https://reactome.org/content/detail/R-HSA-9006931)：CAMERA Direction=Down，FDR=0.745，nGenes=260；ORA上调集合FDR=1、下调集合FDR=0.812。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 25. GALK1 + TAF1

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 59；多样化 Top30 排名 25。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.56**；上游 90% GI 预测区间 [-0.754, -0.367]；跨随机种子 SD 0.01。
- 预测联合效应 -0.862 = 已测单基因和 -0.301 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：GALK1=0.189；TAF1=-0.49，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：GALK1=0.138；TAF1=-0.454。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### GALK1（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P51570；EC 2.7.1.6；[来源1](https://www.uniprot.org/uniprotkb/P51570/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.138，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.165；HCT116 -0.0566；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.89；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.308。
- HCT116/OxR作者模型通量：映射反应数 2；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.132、OxR 0.291。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.94/2.88；H3K27ac 1.57/1.33；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.0624，FDR 0.801；SW480高耐药对亲本：log2FC -0.303，FDR 0.0423。
- GSE76092 log2FC -0.0933，95%均值对比区间 [-0.234, 0.0473]，BH FDR 0.33，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 9.41，强度分位数 0.655；耐药株均值 9.32，强度分位数 0.645。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.103 / FDR 0.609；GSE289766 HCT116 log2FC -0.551；GSE42387 HCT116 log2FC -0.467、HT29 log2FC -0.0579。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -0.33；有效 guide 数 3；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 4 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Metabolism of carbohydrates and carbohydrate derivatives](https://reactome.org/content/detail/R-HSA-71387)：CAMERA Direction=Down，FDR=0.968，nGenes=296；ORA上调集合FDR=1、下调集合FDR=1。
- [Diseases of metabolism](https://reactome.org/content/detail/R-HSA-5668914)：CAMERA Direction=Up，FDR=0.986，nGenes=268；ORA上调集合FDR=1、下调集合FDR=1。
- [Diseases associated with glycosylation precursor biosynthesis](https://reactome.org/content/detail/R-HSA-5609975)：CAMERA Direction=Down，FDR=0.993，nGenes=16；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### TAF1（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed EC 2.3.1.48 and H3/H4 acetyltransferase activity support direct histone modification in addition to its TFIID scaffold and kinase functions; inclusion is not based on generic transcription-factor status.
- UniProt P21675；EC 2.3.1.48; 2.7.11.1；[来源1](https://www.uniprot.org/uniprotkb/P21675/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 11295558。
- 作者control-scaled single -0.454，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.6；HCT116 -0.961；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.154。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.85/2.76；H3K27ac 1.84/1.74；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.161，FDR 0.261；SW480高耐药对亲本：log2FC 0.0653，FDR 0.668。
- GSE76092 log2FC -0.0667，95%均值对比区间 [-0.201, 0.0678]，BH FDR 0.477，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 8.63，强度分位数 0.578；耐药株均值 8.56，强度分位数 0.574。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.166 / FDR 0.238；GSE289766 HCT116 log2FC 0.766；GSE42387 HCT116 log2FC -0.0791、HT29 log2FC -0.18。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 17 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Regulation of TP53 Activity](https://reactome.org/content/detail/R-HSA-5633007)：CAMERA Direction=Up，FDR=0.784，nGenes=161；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of TP53 Activity through Phosphorylation](https://reactome.org/content/detail/R-HSA-6804756)：CAMERA Direction=Up，FDR=0.801，nGenes=93；ORA上调集合FDR=1、下调集合FDR=1。
- [HIV Transcription Initiation](https://reactome.org/content/detail/R-HSA-167161)：CAMERA Direction=Up，FDR=0.854，nGenes=47；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 26. PNP + TAF1

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 397；多样化 Top30 排名 26。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.392**；上游 90% GI 预测区间 [-0.585, -0.198]；跨随机种子 SD 0.0241。
- 预测联合效应 -0.799 = 已测单基因和 -0.407 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：PNP=0.0833；TAF1=-0.49，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：PNP=0.0659；TAF1=-0.454。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### PNP（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P00491；EC 2.4.2.1；[来源1](https://www.uniprot.org/uniprotkb/P00491/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.0659，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.0216；HCT116 -0.145；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 11.1；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.427。
- HCT116/OxR作者模型通量：映射反应数 1e+01；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.234、OxR 0.415。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.04/3.12；H3K27ac 1.61/1.65；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.0387，FDR 0.792；SW480高耐药对亲本：log2FC -0.142，FDR 0.198。
- GSE76092 log2FC -0.03，95%均值对比区间 [-0.165, 0.105]，BH FDR 0.777，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 14.7，强度分位数 0.974；耐药株均值 14.7，强度分位数 0.973。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.22 / FDR 0.14；GSE289766 HCT116 log2FC -1.52；GSE42387 HCT116 log2FC -0.0209、HT29 log2FC 0.43。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 9 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Ribavirin ADME](https://reactome.org/content/detail/R-HSA-9755088)：CAMERA Direction=Up，FDR=0.564，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Purine salvage](https://reactome.org/content/detail/R-HSA-74217)：CAMERA Direction=Up，FDR=0.564，nGenes=13；ORA上调集合FDR=1、下调集合FDR=1。
- [Nucleotide salvage](https://reactome.org/content/detail/R-HSA-8956321)：CAMERA Direction=Up，FDR=0.564，nGenes=23；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### TAF1（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed EC 2.3.1.48 and H3/H4 acetyltransferase activity support direct histone modification in addition to its TFIID scaffold and kinase functions; inclusion is not based on generic transcription-factor status.
- UniProt P21675；EC 2.3.1.48; 2.7.11.1；[来源1](https://www.uniprot.org/uniprotkb/P21675/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 11295558。
- 作者control-scaled single -0.454，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.6；HCT116 -0.961；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.154。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.85/2.76；H3K27ac 1.84/1.74；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.161，FDR 0.261；SW480高耐药对亲本：log2FC 0.0653，FDR 0.668。
- GSE76092 log2FC -0.0667，95%均值对比区间 [-0.201, 0.0678]，BH FDR 0.477，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 8.63，强度分位数 0.578；耐药株均值 8.56，强度分位数 0.574。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.166 / FDR 0.238；GSE289766 HCT116 log2FC 0.766；GSE42387 HCT116 log2FC -0.0791、HT29 log2FC -0.18。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 17 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Regulation of TP53 Activity](https://reactome.org/content/detail/R-HSA-5633007)：CAMERA Direction=Up，FDR=0.784，nGenes=161；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of TP53 Activity through Phosphorylation](https://reactome.org/content/detail/R-HSA-6804756)：CAMERA Direction=Up，FDR=0.801，nGenes=93；ORA上调集合FDR=1、下调集合FDR=1。
- [HIV Transcription Initiation](https://reactome.org/content/detail/R-HSA-167161)：CAMERA Direction=Up，FDR=0.854，nGenes=47；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 27. DHODH + TAF1

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 485；多样化 Top30 排名 27。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.291**；上游 90% GI 预测区间 [-0.484, -0.0976]；跨随机种子 SD 0.0127。
- 预测联合效应 -0.538 = 已测单基因和 -0.247 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：DHODH=0.244；TAF1=-0.49，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：DHODH=0.258；TAF1=-0.454。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 lower_half_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### DHODH（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q02127；EC 1.3.5.2；[来源1](https://www.uniprot.org/uniprotkb/Q02127/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.258，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.437；HCT116 -0.616；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 7.78；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.113。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.33、OxR 0.436。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.25/3.05；H3K27ac 1.29/1.48；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.144，FDR 0.469；SW480高耐药对亲本：log2FC 0.0102，FDR 0.953。
- GSE76092 log2FC 0.127，95%均值对比区间 [0.0252, 0.228]，BH FDR 0.0621，BY FDR 0.654；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 11.8，强度分位数 0.858；耐药株均值 11.9，强度分位数 0.865。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.31 / FDR 0.0609；GSE289766 HCT116 log2FC -1.1；GSE42387 HCT116 log2FC 0.227、HT29 log2FC 0.181。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -1.61；有效 guide 数 2；guide 耗竭比例 1；处理重复耗竭比例 1。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 2 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Nucleotide biosynthesis](https://reactome.org/content/detail/R-HSA-8956320)：CAMERA Direction=Down，FDR=0.968，nGenes=14；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### TAF1（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed EC 2.3.1.48 and H3/H4 acetyltransferase activity support direct histone modification in addition to its TFIID scaffold and kinase functions; inclusion is not based on generic transcription-factor status.
- UniProt P21675；EC 2.3.1.48; 2.7.11.1；[来源1](https://www.uniprot.org/uniprotkb/P21675/entry)。
- EpiFactors 功能：Histone modification write；功能文献 PMID 11295558。
- 作者control-scaled single -0.454，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.6；HCT116 -0.961；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.154。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.85/2.76；H3K27ac 1.84/1.74；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.161，FDR 0.261；SW480高耐药对亲本：log2FC 0.0653，FDR 0.668。
- GSE76092 log2FC -0.0667，95%均值对比区间 [-0.201, 0.0678]，BH FDR 0.477，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 8.63，强度分位数 0.578；耐药株均值 8.56，强度分位数 0.574。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC 0.166 / FDR 0.238；GSE289766 HCT116 log2FC 0.766；GSE42387 HCT116 log2FC -0.0791、HT29 log2FC -0.18。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：up_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 17 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Regulation of TP53 Activity](https://reactome.org/content/detail/R-HSA-5633007)：CAMERA Direction=Up，FDR=0.784，nGenes=161；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of TP53 Activity through Phosphorylation](https://reactome.org/content/detail/R-HSA-6804756)：CAMERA Direction=Up，FDR=0.801，nGenes=93；ORA上调集合FDR=1、下调集合FDR=1。
- [HIV Transcription Initiation](https://reactome.org/content/detail/R-HSA-167161)：CAMERA Direction=Up，FDR=0.854，nGenes=47；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 28. NME2 + HDAC3

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 488；多样化 Top30 排名 28。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.288**；上游 90% GI 预测区间 [-0.481, -0.0946]；跨随机种子 SD 0.0493。
- 预测联合效应 -0.526 = 已测单基因和 -0.238 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HDAC3=-0.386；NME2=0.148，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HDAC3=-0.376；NME2=0.157。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 upper_quartile_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HDAC3（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt O15379；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/O15379/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 10655483。
- 作者control-scaled single -0.376，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.952；HCT116 -0.987；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10.1；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.325。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.6/3.09；H3K27ac 1.62/1.63；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.207，FDR 0.0488；SW480高耐药对亲本：log2FC -0.0255，FDR 0.787。
- GSE76092 log2FC -0.09，95%均值对比区间 [-0.191, 0.0111]，BH FDR 0.18，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 13.2，强度分位数 0.931；耐药株均值 13.2，强度分位数 0.928。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.131 / FDR 0.381；GSE289766 HCT116 log2FC -0.761；GSE42387 HCT116 log2FC -0.386、HT29 log2FC 0.356。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 62 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [RUNX2 regulates osteoblast differentiation](https://reactome.org/content/detail/R-HSA-8940973)：CAMERA Direction=Up，FDR=0.564，nGenes=24；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Transcriptional activation of mitochondrial biogenesis](https://reactome.org/content/detail/R-HSA-2151201)：CAMERA Direction=Down，FDR=0.564，nGenes=54；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### NME2（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt P22392；EC 2.7.13.3; 2.7.4.6；[来源1](https://www.uniprot.org/uniprotkb/P22392/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.157，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.361；HCT116 -0.437；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 13.3；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.252。
- HCT116/OxR作者模型通量：映射反应数 3e+01；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.633、OxR 1.07。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.86/3.69；H3K27ac 2.64/1.82；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.106，FDR 0.756；SW480高耐药对亲本：log2FC -0.594，FDR 0.299。
- GSE76092 log2FC -0.0133，95%均值对比区间 [-0.086, 0.0593]，BH FDR 0.818，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 17.6，强度分位数 0.998；耐药株均值 17.6，强度分位数 0.998。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.285 / FDR 0.494；GSE289766 HCT116 log2FC -1.76；GSE42387 HCT116 log2FC -0.322、HT29 log2FC 0.226。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -1.13；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 0.667。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 6 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Ribavirin ADME](https://reactome.org/content/detail/R-HSA-9755088)：CAMERA Direction=Up，FDR=0.564，nGenes=11；ORA上调集合FDR=1、下调集合FDR=1。
- [Azathioprine ADME](https://reactome.org/content/detail/R-HSA-9748787)：CAMERA Direction=Up，FDR=0.784，nGenes=22；ORA上调集合FDR=1、下调集合FDR=1。
- [Drug ADME](https://reactome.org/content/detail/R-HSA-9748784)：CAMERA Direction=Down，FDR=0.784，nGenes=97；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 29. MVK + HDAC3

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 362；多样化 Top30 排名 29。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.413**；上游 90% GI 预测区间 [-0.606, -0.219]；跨随机种子 SD 0.0527。
- 预测联合效应 -0.736 = 已测单基因和 -0.323 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：HDAC3=-0.386；MVK=0.0628，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：HDAC3=-0.376；MVK=0.0237。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 upper_quartile_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### HDAC3（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt O15379；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/O15379/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 10655483。
- 作者control-scaled single -0.376，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.952；HCT116 -0.987；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10.1；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.325。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.6/3.09；H3K27ac 1.62/1.63；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.207，FDR 0.0488；SW480高耐药对亲本：log2FC -0.0255，FDR 0.787。
- GSE76092 log2FC -0.09，95%均值对比区间 [-0.191, 0.0111]，BH FDR 0.18，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 13.2，强度分位数 0.931；耐药株均值 13.2，强度分位数 0.928。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.131 / FDR 0.381；GSE289766 HCT116 log2FC -0.761；GSE42387 HCT116 log2FC -0.386、HT29 log2FC 0.356。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 62 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [RUNX2 regulates osteoblast differentiation](https://reactome.org/content/detail/R-HSA-8940973)：CAMERA Direction=Up，FDR=0.564，nGenes=24；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Transcriptional activation of mitochondrial biogenesis](https://reactome.org/content/detail/R-HSA-2151201)：CAMERA Direction=Down，FDR=0.564，nGenes=54；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### MVK（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q03426；EC 2.7.1.36；[来源1](https://www.uniprot.org/uniprotkb/Q03426/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.0237，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.697；HCT116 -1.29；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 8.83；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.161。
- HCT116/OxR作者模型通量：映射反应数 1；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.125、OxR 0.174。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 2.81/2.95；H3K27ac 1.88/1.4；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC -0.153，FDR 0.483；SW480高耐药对亲本：log2FC -0.0799，FDR 0.518。
- GSE76092 log2FC 0.227，95%均值对比区间 [0.0963, 0.357]，BH FDR 0.0142，BY FDR 0.15；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 13.7，强度分位数 0.948；耐药株均值 13.9，强度分位数 0.954。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.228 / FDR 0.261；GSE289766 HCT116 log2FC -0.791；GSE42387 HCT116 log2FC -0.00206、HT29 log2FC 0.552。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC -0.483；有效 guide 数 3；guide 耗竭比例 0.667；处理重复耗竭比例 0.667。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 5 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Cholesterol biosynthesis](https://reactome.org/content/detail/R-HSA-191273)：CAMERA Direction=Up，FDR=0.718，nGenes=27；ORA上调集合FDR=1、下调集合FDR=1。
- [Lanosterol biosynthesis](https://reactome.org/content/detail/R-HSA-9969896)：CAMERA Direction=Up，FDR=0.927，nGenes=14；ORA上调集合FDR=1、下调集合FDR=1。
- [Metabolism of steroids](https://reactome.org/content/detail/R-HSA-8957322)：CAMERA Direction=Down，FDR=0.927，nGenes=154；ORA上调集合FDR=1、下调集合FDR=0.118。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。


---

# 30. AK6 + HDAC3

**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。
交接证据层：hypothesis_design_extrapolation_and_unvalidated_resistant_context；需要交接探索：True。双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。

GI 主长名单排名 435；多样化 Top30 排名 30。Priority8：本轮未选入。

## 计算预测与已测单基因

- 预测加性 GI：**-0.341**；上游 90% GI 预测区间 [-0.534, -0.148]；跨随机种子 SD 0.0743。
- 预测联合效应 -0.654 = 已测单基因和 -0.313 + 预测 GI；这是同一预测的代数重构。
- 原始单基因 LFC：AK6=0.0724；HDAC3=-0.386，用于GI和joint的源标尺计算。
- 作者control-scaled单基因：AK6=0.028；HDAC3=-0.376。负控0、正控−1；强耗竭阈值 ≤ -0.5。覆盖2/2臂，已观测达到阈值0臂，缺测0臂。
- 单基因层：no_strong_control_scaled_single_at_cutoff。两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺。
- 不确定性层：interval_entirely_negative；seed SD 位于核心长名单的 upper_quartile_of_core_candidate_seed_spread。GI区间/single数学分层：A_negative_GI_interval_below_zero_no_strong_single；此A–U编码不是生物学证据等级。
- 模型 pooled_multigraph，5 个拟合；设计范围 both_nonanchor_extrapolation。
- 区间来源：Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping。适用范围：Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage。
- 种子差异范围：Training randomness conditional on the same final split; not a calibrated interval。

GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源15,864个已测组合中，99.01%的proxy小于0，均值-0.464、中位数-0.473。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。
原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。

## 两臂功能与耐药表达证据

GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。

### AK6（严格代谢酶）

- 机制：Reviewed EC and catalytic reaction directly establish intermediary metabolite, cofactor, energy or redox catalysis; protein-bound metabolic cofactors are allowed. Not inferred from a pathway membership or catalytic GO alone.
- UniProt Q9Y3D8；EC 2.7.4.3；[来源1](https://www.uniprot.org/uniprotkb/Q9Y3D8/entry)。
- EpiFactors 功能：未提供；功能文献 PMID 未提供。
- 作者control-scaled single 0.028，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.654；HCT116 -1.04；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 未覆盖/未提供；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.429。
- HCT116/OxR作者模型通量：映射反应数 3；生长率归一化后log1p绝对通量界的基因关联均值，亲本 0.921、OxR 1.38。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.11/3.15；H3K27ac 1.95/2.15；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.0842，FDR 0.717；SW480高耐药对亲本：log2FC -0.469，FDR 0.000164。
- GSE76092 log2FC 0.507，95%均值对比区间 [0.358, 0.655]，BH FDR 0.000272，BY FDR 0.00287；标签 resistant_up_BH_lt_005。
- 亲本相对芯片强度均值 12.2，强度分位数 0.881；耐药株均值 12.7，强度分位数 0.906。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.504 / FDR 0.000817；GSE289766 HCT116 log2FC -0.312；GSE42387 HCT116 log2FC 未覆盖/未提供、HT29 log2FC 未覆盖/未提供。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 2 项，方向支持比例 0.5。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 2 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [Metabolism of nucleotides](https://reactome.org/content/detail/R-HSA-15869)：CAMERA Direction=Up，FDR=0.968，nGenes=98；ORA上调集合FDR=1、下调集合FDR=1。
- [Interconversion of nucleotide di- and triphosphates](https://reactome.org/content/detail/R-HSA-499943)：CAMERA Direction=Down，FDR=0.974，nGenes=30；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

### HDAC3（直接组蛋白修饰酶）

- 机制：Reviewed direct histone-modifying catalytic activity; admitted to the primary research scope. Reviewed histone deacetylation reactions and HDAC chromatin complexes establish a direct histone modifier, despite additional non-histone substrates.
- UniProt O15379；EC 3.5.1.-; 3.5.1.98；[来源1](https://www.uniprot.org/uniprotkb/O15379/entry)。
- EpiFactors 功能：Histone modification erase；功能文献 PMID 10655483。
- 作者control-scaled single -0.376，基于12个gene-control guide组合；独立Project SCORE HT29 LFC 未覆盖/未提供，源文件相对rank 未覆盖/未提供。SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。
- DepMap亲本细胞gene effect：HT29 -0.952；HCT116 -0.987；SW480 未覆盖/未提供。该单基因标尺不能替代耐药株或双基因结果。
- TCGA-COAD肿瘤中位表达 10.1；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 0.325。
- HCT116/OxR作者模型通量：映射反应数 未覆盖/未提供；生长率归一化后log1p绝对通量界的基因关联均值，亲本 未覆盖/未提供、OxR 未覆盖/未提供。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。
- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：H3K4me3 3.6/3.09；H3K27ac 1.62/1.63；H3K27me3 0/0。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。
- GSE299428 HCT116高耐药对亲本：log2FC 0.207，FDR 0.0488；SW480高耐药对亲本：log2FC -0.0255，FDR 0.787。
- GSE76092 log2FC -0.09，95%均值对比区间 [-0.191, 0.0111]，BH FDR 0.18，BY FDR 1；标签 measured_without_BH_lt_005_association。
- 亲本相对芯片强度均值 13.2，强度分位数 0.931；耐药株均值 13.2，强度分位数 0.928。相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。
- 其他表达：GSE299428 HCT116 MD log2FC -0.131 / FDR 0.381；GSE289766 HCT116 log2FC -0.761；GSE42387 HCT116 log2FC -0.386、HT29 log2FC 0.356。后两项来源只提供此处使用的效应量。
- HCT116 跨研究方向：down_in_resistance；覆盖 3 项，方向支持比例 0.333。
- GSE254136 单基因奥沙利铂筛选 log2FC 未覆盖/未提供；有效 guide 数 未覆盖/未提供；guide 耗竭比例 未覆盖/未提供；处理重复耗竭比例 未覆盖/未提供。负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。

**主耐药对比的通路成员证据**

共有 62 条已检验成员通路，其中 0 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。
- [RUNX2 regulates osteoblast differentiation](https://reactome.org/content/detail/R-HSA-8940973)：CAMERA Direction=Up，FDR=0.564，nGenes=24；ORA上调集合FDR=1、下调集合FDR=1。
- [Regulation of MECP2 expression and activity](https://reactome.org/content/detail/R-HSA-9022692)：CAMERA Direction=Down，FDR=0.564，nGenes=31；ORA上调集合FDR=1、下调集合FDR=1。
- [Transcriptional activation of mitochondrial biogenesis](https://reactome.org/content/detail/R-HSA-2151201)：CAMERA Direction=Down，FDR=0.564，nGenes=54；ORA上调集合FDR=1、下调集合FDR=1。

成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。

表达来源：[GSE76092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092)、[GSE299428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428)、[GSE289766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289766)、[GSE42387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387)；药筛来源：[GSE254136](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136)。

主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。

## 两臂共享的主对比通路

本次主对比的已检验 Reactome 集合中没有成员关联记录。

共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。

## 独立研究中的同一组合

当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。

## 合作方待验证的问题与反馈字段

- 操作性信息：需合作方确认。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。
- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。
- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。
- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。
- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。

