# 完整获取与复现

本仓库同时提供已完成的参考结果和生成这些结果的实现。原始数据快照、处理后输入、训练记录、模型权重、候选清单和报告均随仓库保存；大文件使用Git LFS。

## 获取项目

```bash
git lfs install
git clone https://github.com/Tairan-Terrian/MetaGene.git
cd MetaGene
git lfs pull
```

医学院可直接阅读[实验报告](../results/MetaGene实验汇报.md)、[Top500](../results/candidates/top500.csv)和[交接说明](CANDIDATE_HANDOFF.md)。研究者完整克隆后，可使用已有模型继续分析，也可在独立目录重新训练。

## 安装环境与核对输入

本版计算环境为Python 3.10.20、PyTorch 2.6.0+cu124、CUDA 12.4和两张RTX 3090。其余依赖固定在[复现依赖](../requirements-reproduce.txt)。新建环境后运行：

```bash
python -m pip install -r requirements-reproduce.txt
python scripts/reproduce.py --check-only
```

输入检查对照[九个固定输入的SHA256清单](../data/reference/input_manifest.json)。这些文件包含HT29标签源表、基因顺序、特征、网络、独立测试数据和候选静态注释，合计约16 MB；它们不包含预训练参数。数据角色、公开来源和人工审定记录见[数据溯源说明](DATA_PROVENANCE.md)。

## 从头训练当前版本

```bash
python scripts/reproduce.py --output-dir reproductions/sl_v2_rerun --devices cuda:0 cuda:1
```

只有一张GPU时使用`--devices cuda:0`；CPU使用`--devices cpu`。省略设备参数会自动选择最多两张CUDA卡，无CUDA时使用CPU。

该命令把代码、固定输入和配置复制到新目录，重新构建标签，然后依次执行基准训练、模型选择、最终集成训练、主动学习、候选交付、外部评价、完整性检查、图表与报告生成。参考配置对应66次基准拟合、10次最终拟合和156次主动学习拟合。

新目录从空训练缓存开始，不复制仓库中的模型权重或既有运行记录。产物位于：

- `reproductions/sl_v2_rerun/results/MetaGene实验汇报.md`
- `reproductions/sl_v2_rerun/results/candidates/`
- `reproductions/sl_v2_rerun/results/sl/`
- `reproductions/sl_v2_rerun/reproduction.log`

输入、分割、种子和依赖版本固定；运行结果保留实际设备及环境信息。不同GPU、CPU或数值库可能产生浮点差异，比较时同时检查候选顺序、指标和完整性审计。

中断后，在相同代码、输入和配置下恢复：

```bash
python scripts/reproduce.py --output-dir reproductions/sl_v2_rerun --devices cuda:0 cuda:1 --resume
```

## 修改模型或配置

修改`metagene/sl.py`及相关模型代码后，指定新的输出目录运行即可。配置实验可以复制并编辑[实验配置](../configs/experiment_sl.json)：

```bash
python scripts/reproduce.py --config configs/my_experiment.json --output-dir reproductions/my_experiment --devices cuda:0
```

若有意修改固定输入，增加`--allow-modified-inputs`，新运行会记录修改后的校验值。代码、输入或配置变化后，旧拟合缓存不会被自动复用。仓库中的参考报告和候选结果保持原样。

仅重新整理参考阅读文档可运行`python scripts/write_sl_documents.py`。完整原始数据与上游处理脚本也保留在仓库中，更新数据库版本时按数据溯源说明追踪输入变化。
