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
