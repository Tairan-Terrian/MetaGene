"""Write current reading documents from compact results, without refitting models."""
from pathlib import Path
import csv
import json
import statistics

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/sl'
PROJECT_TITLE = 'AI主动学习驱动的“代谢-表观遗传”跨维度合成致死聚焦筛选与奥沙利铂耐药逆转'
REPOSITORY_URL = 'https://github.com/Tairan-Terrian/MetaGene'


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def collect_fit_counts(out, cfg, selection):
    """Count completed fits, using compact evaluation tables for reading copies."""
    expected = [(p, m, s) for p in cfg['protocols'] for m in cfg['models']
                for s in (cfg['seeds'] if p == 'pair_holdout' else cfg['secondary_seeds'])]
    final_jobs = [(m, s) for m in dict.fromkeys(
        [selection['selected_gnn'], selection['best_unrestricted_model']]) for s in cfg['seeds']]
    ac = cfg['active_learning']
    run_root = out / 'runs'
    has_runs = run_root.is_dir() and any(run_root.glob('*/run.json'))
    if has_runs:
        benchmark = sum((run_root / f'{p}_{m}_{s}' / 'run.json').exists() for p, m, s in expected)
        final = sum((run_root / f'final_{m}_{s}' / 'run.json').exists() for m, s in final_jobs)
        active = sum(len(rows(run_root / f'active_{strategy}_{seed}' / 'curves.csv'))
                     for strategy in ac['strategies'] for seed in ac['seeds']
                     if (run_root / f'active_{strategy}_{seed}' / 'curves.csv').exists())
        source = 'completed run.json files and recorded active-learning curve rows'
    else:
        metrics_path = out / 'validation_metrics.csv'
        completed = {(r['protocol'], r['model'], int(r['seed'])) for r in rows(metrics_path)} if metrics_path.exists() else set()
        benchmark = len(set(expected) & completed)
        # A compact delivery retains the ensemble result, not every fit directory.
        final = len(final_jobs) if (out / 'delivery_review.json').exists() else 0
        curve_path = out / 'active_learning_curves.csv'
        active = sum(r['strategy'] in ac['strategies'] and int(r['seed']) in ac['seeds']
                     and int(r['round']) <= ac['rounds'] for r in rows(curve_path)) if curve_path.exists() else 0
        source = 'validation rows, completed delivery ensemble configuration, and active-learning curve rows'
    counts = {'benchmark': benchmark, 'final': final, 'active_learning': active}
    return counts, source


def display_model(name):
    return {'multigraph': '多图GNN', 'string': 'STRING图GNN', 'xgb': 'XGBoost',
            'linear': '线性模型', 'graphless': '无图模型', 'regression_only': '仅回归模型'}.get(name, name)


def compact_number(value):
    return f'{float(value):.2f}'.rstrip('0').rstrip('.')


def table(headers, values):
    return ('| ' + ' | '.join(headers) + ' |\n|' + '|'.join(['---']*len(headers)) + '|\n'
            + '\n'.join('| ' + ' | '.join(str(v).replace('|',' / ') for v in row) + ' |' for row in values))


def render_candidate_cards(records):
    def number(value):
        try:
            x=float(value)
            return f'{x:.3f}' if x==x else '待补充'
        except (TypeError,ValueError):return '待补充'
    text=['# 干实验交付：Top30逐对证据卡\n',
          f'所属课题：**{PROJECT_TITLE}**。\n',
          f'交付入口：[MetaGene GitHub仓库]({REPOSITORY_URL})。可直接在网页阅读候选及解释，无需运行模型。\n',
          '以下为HT29筛选背景下的模型预测候选。预测dLFC<0表示双干预低于单基因加性预期；强协同事件分数对应dLFC<−1。合作方回传测量用于评价实际配对效应。\n']
    for row in records:
        rationale=row['biological_rationale'].split('当前单基因审查为')[0]
        coverage={'not_covered_unknown':'现有独立组合表未覆盖'}.get(row['external_pair_state'],row['external_pair_state'])
        text += [f"## {row['top30_rank']}. {row['metabolic_gene']} × {row['chromatin_gene']}\n",
                 f"Top500第{row['top500_rank']}；{row['single_effect_status']}。\n",rationale+'\n',
                 f"- 预测dLFC：{number(row['pred_dlfc'])}；名义90%区间 [{number(row['dlfc_pi90_lower'])}, {number(row['dlfc_pi90_upper'])}]；强协同事件分数：{float(row['synergy_probability']):.5f}。",
                 f"- 作者control-scaled单基因：{row['gene_a']}={number(row['gene_a_colo1_single_control_scaled'])}；{row['gene_b']}={number(row['gene_b_colo1_single_control_scaled'])}。",
                 f"- 耐药表达关联背景数：{int(row['target_resistance_association_count'])}；外部组合数据：{coverage}。",
                 f"- 功能来源：[代谢酶]({row['metabolic_function_source']})；[组蛋白修饰酶]({row['histone_function_source']})。\n"]
    return '\n'.join(text)+'\n'


def main():
    candidates = ROOT/'results/candidates'
    top = rows(candidates/'top500.csv')
    short = rows(candidates/'top30.csv')
    representatives = rows(ROOT/'results/handoff/first_round_representatives.csv')
    cfg = json.loads((ROOT/'configs/experiment_sl.json').read_text())
    review = json.loads((OUT/'delivery_review.json').read_text())
    labels = json.loads((ROOT/'data/sl/label_definition.json').read_text())
    metrics = rows(OUT/'ensemble_metrics.csv')
    ac = cfg['active_learning']
    active = [r for r in rows(OUT/'active_learning_summary.csv')
              if int(r['round']) == ac['rounds'] and r['strategy'] in ac['strategies']]
    external = rows(OUT/'external_metrics.csv')
    ci = next(r for r in json.loads((OUT/'bootstrap_intervals.json').read_text()) if r['protocol']=='pair_holdout')
    audit = json.loads((OUT/'integrity_audit.json').read_text())
    delta = [float(r['pred_dlfc']) for r in top]
    summary = {
        'version':cfg['version'], 'n_candidates':len(top),
        'predicted_dlfc_below_zero':sum(x<0 for x in delta),
        'predicted_dlfc_min':min(delta), 'predicted_dlfc_max':max(delta),
        'predicted_dlfc_median':statistics.median(delta),
        'predicted_dlfc_below_minus1':sum(x< -1 for x in delta),
        'direction_definition':'Negative predicted dLFC means predicted double-perturbation fitness below the additive expectation.',
        'training_definition':'The frozen main classifier and reported classification benchmarks retain dLFC < -1.',
        'result_type':'Model predictions; pair effects are to be measured by the collaborating laboratory.',
        'prediction_values_and_rankings_changed':False,
    }
    (OUT/'prediction_direction_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    selection = json.loads((OUT/'model_selection.json').read_text())
    selected = selection['selected_gnn']
    model_name = display_model(selected)
    best_model_name = display_model(selection['best_unrestricted_model'])
    primary = next(r for r in metrics if r['protocol']=='pair_holdout' and r['model']==selected)
    calibration = review['calibration']
    fit_counts, fit_count_source = collect_fit_counts(OUT, cfg, selection)
    total_fits = sum(fit_counts.values())
    ann = rows(ROOT/'data/evidence/genes_annotation.csv')
    n_metabolic = sum(int(float(r['is_strict_metabolic_enzyme'])) for r in ann)
    n_histone = sum(int(float(r['is_direct_histone_enzyme'])) for r in ann)
    candidate_space = n_metabolic * n_histone
    candidate_universe = int(review['candidate_universe'])
    previously_assayed = candidate_space - candidate_universe
    tiers = {tier: sum(int(float(r['single_effect_tier'])) == tier for r in top) for tier in range(3)}
    background_count = sum(float(r['target_resistance_association_count']) > 0 for r in top)
    external_count = sum(float(r['external_measurement_count']) > 0 for r in top)
    n_top, n_short, n_representatives = len(top), len(short), len(representatives)
    n_priority = len(rows(candidates/'priority8.csv'))
    feature_meta = json.loads((ROOT/'data/evidence/integrated_input_manifest.json').read_text())
    coverage = rows(OUT/'core_measured_label_coverage.csv')
    one_replicate = sum(int(float(r['n_replicates'])) == 1 for r in coverage)
    no_replicate = sum(int(float(r['n_replicates'])) == 0 for r in coverage)
    fractions = cfg['fractions']
    split_fractions = '/'.join(f'{x:.0%}' for x in [1-sum(fractions.values()), fractions['validation'], fractions['calibration'], fractions['test']])
    negative_count = sum(x < 0 for x in delta)
    direction_text = '预测dLFC全部小于0' if negative_count == n_top else f'{negative_count}/{n_top}对预测dLFC小于0'
    direction_interpretation = (f'这意味着当前模型对这{n_top}对组合均预测出负向协同效应'
                                if negative_count == n_top else f'这意味着当前模型对其中{negative_count}对组合预测出负向协同效应')
    range_text = ('Top500的预测均值处于−1与0之间，表示模型预测有协同抑制。'
                  if all(-1 < x < 0 for x in delta) else
                  f'Top500中{negative_count}对预测均值小于0，其中{sum(x < -1 for x in delta)}对小于−1。')
    upper = next((r for r in active if r['strategy']=='upper_confidence'), None)
    chosen_active = upper or max(active, key=lambda r:float(r['strong_hits_mean']))
    strategy_name = {'upper_confidence':'上界策略','probability':'事件分数优先策略',
                     'uncertainty':'不确定性优先策略','random':'随机策略'}.get(chosen_active['strategy'],chosen_active['strategy'])
    random = next((r for r in active if r['strategy']=='random'), None)
    budget = float(chosen_active['n_measured'])
    active_hits = float(chosen_active['strong_hits_mean'])
    active_brief = f'{strategy_name}在{budget:,.0f}条测量记录下平均发现{compact_number(active_hits)}个强协同信号'
    active_comparison = ''
    if random is not None and random is not chosen_active:
        random_hits = float(random['strong_hits_mean'])
        same_budget = float(random['n_measured']) == budget
        active_brief += (f'，同预算随机策略约{compact_number(random_hits)}个' if same_budget else
                         f'；随机策略在{float(random["n_measured"]):,.0f}条记录下约{compact_number(random_hits)}个')
        if same_budget and random_hits > 0:
            active_comparison = f'，是随机策略的约{compact_number(active_hits/random_hits)}倍'
    active_result = f'{active_brief}{active_comparison}。'
    best_test = max((r for r in metrics if r['protocol']=='pair_holdout'),key=lambda r:float(r['auprc']))
    ranking_comparison = f'{display_model(best_test["model"])}在该固定测试集上的AUPRC最高。'
    regression = next((r for r in metrics if r['protocol']=='pair_holdout' and r['model']=='regression_only'),None)
    if regression is not None and float(primary['auprc']) > float(regression['auprc']) and float(primary['hits_at_500']) > float(regression['hits_at_500']):
        ranking_comparison += f'{model_name}相较仅回归模型具有更高的AUPRC和Top500回收数。'
    slots = rows(ROOT/'results/handoff/1000_slot_capacity_draft.csv')
    slot_counts = {kind:sum(r['slot_type']==kind for r in slots) for kind in {r['slot_type'] for r in slots}}
    external_coverage_text = ('Top500目前在外部组合表中未被覆盖' if not external_count else f'Top500目前有{external_count}对被外部组合表覆盖')
    pair_table = table(['模型','AUPRC','AUROC','Top50强协同命中','Top500强协同命中','dLFC RMSE'],
        [[r['model'],f"{float(r['auprc']):.4f}",f"{float(r['auroc']):.4f}",r['hits_at_50'],r['hits_at_500'],f"{float(r['dlfc_rmse']):.4f}"]
         for r in metrics if r['protocol']=='pair_holdout'])
    generalization = table(['留出方式','测试配对','强协同阳性','AUPRC','AUROC'],
        [[r['protocol'],r['n'],r['positives'],f"{float(r['auprc']):.4f}",f"{float(r['auroc']):.4f}"] for r in metrics if r['model']==selected])
    active_table = table(['采集策略','累计测量数','dLFC<−1命中均值 ± 种子SD','dLFC<−0.75命中均值'],
        [[r['strategy'],int(float(r['n_measured'])),f"{float(r['strong_hits_mean']):.2f} ± {float(r['strong_hits_sd']):.2f}",f"{float(r['tail_hits_mean']):.2f}"] for r in active])
    short_table = table(['Top30序号','代谢酶','组蛋白修饰酶','预测dLFC','强协同事件分数'],
        [[r['top30_rank'],r['metabolic_gene'],r['chromatin_gene'],f"{float(r['pred_dlfc']):.3f}",f"{float(r['synergy_probability']):.4f}"] for r in short[:10]])
    external_table = table(['独立数据集','未在HT29出现的配对数','源定义阳性','AUPRC','AUROC'],
        [[r['dataset'],r['n'],r['positives'],f"{float(r['auprc']):.4f}",f"{float(r['auroc']):.4f}"] for r in external if r['scope']=='pair_not_in_ht29'])
    representative_names = '、'.join(f"{r['metabolic_gene']} × {r['chromatin_gene']}" for r in representatives)
    reference_representatives = {('UCK2','KDM3B'),('PKM','HDAC2'),('HK2','KAT6B')}
    representative_rationale = 'UCK2参与嘧啶核苷补救，KDM3B参与组蛋白去甲基化，组合假说连接核苷酸稳态与染色质适应；PKM与HK2参与糖酵解，分别与HDAC2的去乙酰化和KAT6B的乙酰化调控配对，形成代谢供给与染色质调节的双重干预方向。' if {(r['metabolic_gene'],r['chromatin_gene']) for r in representatives} == reference_representatives else '；'.join(r.get('biological_rationale','').split('当前单基因审查为')[0] for r in representatives)
    version = cfg['version']
    report = f'''# {PROJECT_TITLE}

**干实验阶段报告｜AI模型、主动学习与基因组合交付**

计算流程名称：MetaGene。实验版本：{version}。

## 1. 课题目标与本次干实验成果

总课题围绕结直肠癌奥沙利铂耐药，研究代谢重编程与表观遗传重塑之间的组合脆弱性。通过AI主动学习优选跨通路合成致死组合，衔接聚焦筛选与代表组合验证，评价双基因干预能否恢复奥沙利铂敏感性。

整体研究路径为：**多源数据整合 → AI预测与主动学习 → Top500基因组合交付 → 医学院聚焦筛选与耐药逆转验证 → 实测结果回传与模型更新**。计算团队负责模型、候选排序、解释和更新；医学院团队负责后续湿实验及功能验证。

本次干实验已完成“**同时干预代谢酶与组蛋白修饰酶，预测超出单基因加性预期的协同抑制，形成潜在合成致死候选**”的计算流程，交付Top500排序、Top30短名单、逐对生物学解释及结果回传接口。下文汇报该计算阶段的方法与结果。

**Top500的{direction_text}，范围为{min(delta):.3f}至{max(delta):.3f}，中位数为{statistics.median(delta):.3f}。** {direction_interpretation}：双干预后的适应性低于两单基因效应的加性预期。这里评价的是相对预期的协同恶化，无须设定绝对致死水平。

| 计算产出 | 结果 |
|---|---|
| 研究对象 | {n_metabolic}个代谢酶 × {n_histone}个直接组蛋白修饰酶 |
| 候选空间 | {candidate_space}种配对，排除{previously_assayed}种既有设计后进入排序{candidate_universe}对 |
| Top500预测方向 | {negative_count}/{n_top}对预测dLFC<0 |
| 交接层级 | Top500（{n_top}对）→ Top30（{n_short}对）→ Priority8（{n_priority}对）→ {n_representatives}对首轮代表组合 |
| 首轮代表组合 | {representative_names} |
| 主动学习 | {ac['rounds']}轮；{active_brief} |
| 计算规模 | {fit_counts['benchmark']}次基准拟合、{fit_counts['final']}次最终拟合、{fit_counts['active_learning']}次主动学习拟合 |
| 完整性检查 | {audit['count']}项通过 |

阅读入口：[Top500](candidates/top500.csv)、[Top30](candidates/top30.csv)、[逐对证据卡](candidates/candidate_evidence_cards.md)、[首轮代表组合](handoff/first_round_representatives.csv)、[字段字典](candidates/data_dictionary.csv)。每对Top500均附催化功能、组合机制假说、细胞背景信息及来源链接。

## 2. 协同效应的定义

当LFC按共同对照中心化时：

`dLFC = LFC_AB − (LFC_A + LFC_B)`

**dLFC<0表示负向协同，数值越负，双干预相对于加性预期的额外抑制越强。** 单基因效应在加性期望中显式扣除，因而分析对象是组合带来的额外效应。

本研究保留三个不同用途的量：

| 量或阈值 | 用途 |
|---|---|
| 预测dLFC<0 | 判断模型预测的协同方向 |
| dLFC<−1 | 本版主分类头与分类评价使用的强协同事件 |
| dLFC<−0.75 | 辅助分类头与阈值敏感性分析 |
| synergy_probability | 校准后的HT29 dLFC<−1事件分数，用于候选优先级 |

{range_text}−1用于区分本次训练定义的强协同事件；候选的协同方向按0判断。训练标签、原有评价指标和候选数值均沿用已完成的实验配置。

## 3. 数据处理与模型

### 3.1 重建加性残差

HT29原始LFC以质粒为参照。本次按生物重复和guide位置匹配单基因对照，并加回共同基线：

`dLFC_r = LFC_AB,r − median(LFC_A+control,r) − median(LFC_control+B,r) + median(LFC_control+control,r)`

双非必需基因对照的三个重复中位数分别为{', '.join(f'{v:.4f}' for v in labels['baseline_nonessential_by_replicate'].values())}。先在每个重复内计算残差、平均互换方向，再对有效生物重复求均值；监督标签要求至少2个有效重复。缺测保留为空。

最终监督集合包含{labels['model_pairs']:,}对组合，其中{labels['positive_model_pairs']}对达到dLFC<−1，{labels['replicate_significant_model_pairs']}对同时具有重复层面单侧t检验BH q<0.1支持。该统计注释基于聚合基因LFC的重复差异。标签定义见[来源与处理说明](../data/sl/label_definition.json)，不同基线和阈值的计数见[敏感性分析](sl/label_threshold_sensitivity.csv)。

### 3.2 多源表示与对称模型

模型整合{feature_meta['genes']}个基因的{feature_meta['features']:,}维特征：Reactome功能注释、DepMap单敲依赖、TCGA-COAD表达、CRC单细胞表达、HCT116/OxR公开模型通量范围和组蛋白峰描述。STRING、Reactome、BioPlex提供外部网络。基因特征标准化仅使用训练配对涉及的基因。

配对表示由两臂嵌入的和、绝对差及逐元素积构成，交换基因顺序保持预测一致。{model_name}联合学习强协同分类、辅助分类和连续dLFC，损失为加权BCE + {cfg['regression_weight']:g} Smooth-L1 + {cfg['ranking_weight']:g}配对排序损失。隐藏维{cfg['hidden_dim']}、dropout {cfg['dropout']:g}、学习率{cfg['learning_rate']:g}、weight decay {cfg['weight_decay']:g}、batch {cfg['batch_size']}，最多{cfg['max_epochs']}轮，{cfg['patience']}轮早停。

### 3.3 分割与模型选择

主协议按{split_fractions}划分训练、验证、校准和测试集；固定分割种子{cfg['partition_seed']}。主测试有{int(primary['n']):,}对组合、{primary['positives']}个强协同阳性；校准集有{calibration['n']:,}对、{calibration['positives']}个强协同阳性。未见配对采用{len(cfg['seeds'])}个训练种子，未见基因和未见锚点采用{len(cfg['secondary_seeds'])}个种子。基因留出要求每个测试配对至少一臂未在训练中出现。

按主验证集跨种子平均AUPRC选择{model_name}；全部模型中验证最优为{best_model_name}。测试评价固定后，对选中模型以原训练+测试集合重新拟合，每个模型{len(cfg['seeds'])}个种子，保留验证和校准集合。候选分数由{len(cfg['seeds'])}个模型集成并进行Platt校准；dLFC名义90%预测区间来自校准集绝对残差，种子SD表示优化随机性。具体配置见[冻结配置](../configs/experiment_sl.json)。

![标签与重复支持](sl/figures/01_labels.png)

## 4. 模型评价

以下指标针对**dLFC<−1的强协同识别任务**，用于评价模型在已测配对中发现强信号的能力。候选Top500的dLFC<0方向判断与该分类评价分别报告。

{pair_table}

{model_name}的AUPRC为{float(primary['auprc']):.4f}，AUROC为{float(primary['auroc']):.4f}；测试集中前50对找回{primary['hits_at_50']}个强协同阳性，前500对找回{primary['hits_at_500']}个（测试集共{primary['positives']}个）。{ranking_comparison}逐次结果见[测试指标](sl/test_metrics.csv)，跨种子均值与SD见[汇总](sl/benchmark_summary.csv)。

{generalization}

主配对留出的基因簇bootstrap AUPRC 95%区间为[{ci['lower_95']:.4f}, {ci['upper_95']:.4f}]。对应重采样次数和子集结果分别保存在[区间文件](sl/bootstrap_intervals.json)及[子集评价](sl/subgroup_metrics.csv)。

![排序能力](sl/figures/02_ranking.png)

![模型比较](sl/figures/03_benchmarks.png)

独立研究采用各来源的命中定义。下表列出未在HT29原始设计中出现的配对；全部可用配对的结果一并保存在[外部评价](sl/external_metrics.csv)。

{external_table}

## 5. 主动学习

在固定训练池的已有真实标签上进行{ac['rounds']}轮回顾性查询模拟。每个种子下，各策略从相同的{ac['initial_fraction']:.0%}初始集合出发，每轮增加{ac['query_fraction']:.0%}测量预算；本次策略为{'、'.join(ac['strategies'])}。每轮重新拟合，最多{min(cfg['max_epochs'],60)}轮、10轮早停，MC-dropout {ac['mc_passes']}次；验证、校准和测试记录独立保留。

{active_table}

{active_result}该结果用于评价主动学习安排后续组合查询的能力。图中阴影为{len(ac['seeds'])}个种子的SD；完整曲线见[主动学习结果](sl/active_learning_curves.csv)。

![主动学习](sl/figures/04_active_learning.png)

## 6. 候选排序与生物学解释

候选先按单基因效应条件分层，再按强协同事件分数降序、预测dLFC升序和配对键排序；另保留纯模型分数顺序`synergy_rank`。Top30限制每个基因最多出现3次，Priority8最多2次，首轮3对互不共享基因。

| 单基因分层 | Top500数量 | 交接用途 |
|---|---:|---|
| 两臂有测量且均未触发强耗竭阈值 | {tiers[0]} | 优先讨论潜在合成致死配对 |
| 至少一臂单基因测量待补充 | {tiers[1]} | 结合单基因条件开展组合判断 |
| 至少一臂存在强单基因效应 | {tiers[2]} | 协同抑制／合成病弱扩展 |

单基因阈值采用作者control-scaled标尺≤{cfg['candidate_policy']['strong_single_cutoff']:g}。对应文件为[单基因初筛兼容候选](candidates/single_compatible_candidates.csv)、[单基因测量待补充候选](candidates/single_evidence_missing_candidates.csv)和[强单基因背景扩展](candidates/strong_single_synergy_extension.csv)。

{short_table}

首轮代表组合为 **{representative_names}**。{representative_rationale}逐对说明将已知催化功能、预测效应和机制假说分别写明，并链接功能来源。

Top500中{background_count}对具有HCT116/SW480耐药表达关联。DepMap、肿瘤/单细胞表达、模型通量与组蛋白峰提供对应背景信息。耐药关联及单基因数据作为候选讨论依据，配对效应通过合作方的结果回传继续评价。

![候选分层](sl/figures/05_candidates.png)

## 7. 计算结果与课题衔接

本版交付的是**模型预测的潜在合成致死／协同抑制基因组合**，用于总课题的聚焦筛选阶段。训练来自HT29基础生长筛选；{len(coverage)}种既有核心设计中，{one_replicate}种只有1个有效匹配重复，{no_replicate}种缺少有效残差，核心配对因此属于模型外推应用（[覆盖表](sl/core_measured_label_coverage.csv)）。主测试和校准分别有{primary['positives']}个、{calibration['positives']}个强协同阳性，{external_coverage_text}。医学院后续在耐药及匹配亲本模型中评价组合效应与奥沙利铂敏感性变化，并将实测结果回传计算团队，推进主动学习闭环。

## 8. 交付与复现

医学院交接入口为[MetaGene GitHub仓库]({REPOSITORY_URL})。建议依次阅读仓库README、本报告、Top500／Top30及逐对证据卡，再查看交接说明与结果回传模板。网页阅读无需运行模型。

交付包含Top500、Top30、Priority8、首轮代表组合、逐对解释、[500配对靶点表](handoff/500_pair_targets.csv)、[1,000容量草案](handoff/1000_slot_capacity_draft.csv)及[真实结果回传模板](handoff/measured_feedback_template.csv)。容量草案安排{slot_counts.get('primary_candidate_capacity',0)}个首设计、{slot_counts.get('secondary_design_capacity',0)}个次设计和{slot_counts.get('control_or_single_capacity_reservation',0)}个控制占位；guide由合作方衔接。

回传内容包括细胞模型、亲本/耐药状态、药物条件、重复、单A/单B/双干预及对照LFC、标尺、质控和来源。后续先评价当前预测在新测量上的表现，再将新标签纳入下一轮学习。

{audit['count']}项完整性检查覆盖标签代数、分割隔离、配对交换对称性、候选唯一性、类别、排序、多样性及冻结输入，详见[审计记录](integrity_audit.json)。GitHub提供代码、配置、报告、数据、模型权重和运行结果，大文件使用Git LFS，获取方式见[仓库README]({REPOSITORY_URL}#复现与存储)。以下为计算阶段顺序；设备选择与独立运行目录请按[复现说明]({REPOSITORY_URL}/blob/main/docs/REPRODUCIBILITY.md)执行：

```bash
python scripts/prepare_sl_labels.py
python -u scripts/run_sl_experiments.py --phase all
python scripts/deliver_sl.py
python scripts/analyze_sl.py
python scripts/audit_sl.py
python scripts/make_sl_report.py
python scripts/finalize_sl_release.py
```

GitHub网页可直接阅读结果，也可运行`python scripts/write_sl_documents.py`重新生成阅读文档。报告中的拟合次数汇总已有运行记录；修改模型或配置后，按复现说明在独立运行目录训练并生成新报告。交接直接使用GitHub，不生成交接压缩包。

## 9. 主要来源

- [HT29组合筛选原始研究](https://www.nature.com/articles/s41467-025-67256-9)：原始LFC、实验背景及源分析方法。
- [In4mer原始研究](https://www.nature.com/articles/s41467-024-47795-3)：加性期望、dLFC与源研究评分。
- [组合筛选评分方法研究（2024预印本）](https://pmc.ncbi.nlm.nih.gov/articles/PMC11370406/)：空模型与归一化方法。
- [EpiFactors](https://epifactors.autosome.org/)及候选表中的UniProt链接：组蛋白修饰及催化功能注释。
'''
    readme = f'''# {PROJECT_TITLE}

本课题面向结直肠癌奥沙利铂耐药，围绕代谢重编程与表观遗传重塑之间的组合脆弱性，开展AI主动学习引导的跨通路合成致死聚焦筛选，并通过后续功能验证评价耐药逆转。

**MetaGene是本课题的干实验计算流程。** 计算团队负责整合公开数据、构建预测模型、评价主动学习策略，并交付排序后的**Top500代谢酶×组蛋白修饰酶潜在合成致死基因组合**及逐对解释。医学院团队承接聚焦筛选和代表组合验证，评价组合效应与奥沙利铂敏感性变化，并回传结果供模型更新。当前干实验版本：**{version}**。

整体研究路径：**多源数据整合 → AI预测与主动学习 → 基因组合交付 → 聚焦筛选与耐药逆转验证 → 结果回传与模型更新**。

## 干实验主要结果

**Top500预测协同方向：{negative_count}/{n_top}对的预测dLFC<0，范围为{min(delta):.3f}至{max(delta):.3f}，中位数{statistics.median(delta):.3f}。** 该结果表示模型预测的双干预效应低于两单基因效应的加性预期。

- 对{n_metabolic}个代谢酶与{n_histone}个直接组蛋白修饰酶形成的{candidate_universe}对未测候选排序，输出Top500、Top30、Priority8及{n_representatives}对代表组合。
- 每对Top500提供预测dLFC、排序分数、单基因状态、背景证据、来源链接和一段生物学解释。
- {ac['rounds']}轮回顾性主动学习中，{active_result}
- 已有运行记录汇总{total_fits}次模型拟合及{audit['count']}项完整性检查。

## 阅读与交接入口

医学院直接访问[MetaGene GitHub仓库]({REPOSITORY_URL})即可阅读交付内容，无需运行模型。建议顺序：**本README → [实验汇报](results/MetaGene实验汇报.md) → [Top500](results/candidates/top500.csv)／[Top30](results/candidates/top30.csv)／[证据卡](results/candidates/candidate_evidence_cards.md) → [交接说明](docs/CANDIDATE_HANDOFF.md)／[回传模板](results/handoff/measured_feedback_template.csv)**。

| 文件 | 内容 |
|---|---|
| [实验汇报](results/MetaGene实验汇报.md) | 总课题目标与干实验阶段的方法、结果、候选及图表 |
| [Top500](results/candidates/top500.csv) | 排序后的潜在合成致死／协同抑制候选及逐对解释 |
| [Top30](results/candidates/top30.csv) / [证据卡](results/candidates/candidate_evidence_cards.md) | 用于首批候选讨论的多样化短名单 |
| [首轮3对代表组合](results/handoff/first_round_representatives.csv) | {representative_names} |
| [字段字典](results/candidates/data_dictionary.csv) | 分数、区间、单基因条件及来源字段 |
| [研究设计](docs/RESEARCH_DESIGN.md) / [交接说明](docs/CANDIDATE_HANDOFF.md) | 总课题路径、干实验四模块、团队分工与结果回传 |
| [审计记录](results/integrity_audit.json) | 标签、分割、对称性、排序等检查 |

## 如何读懂协同信号

中心化LFC下，`dLFC = LFC_AB − LFC_A − LFC_B`。**dLFC<0表示相对于加性预期的额外抑制，越负表示协同恶化越强。** 对原始未中心化数据，逐重复加回双对照基线。本项目以这一相对效应筛选候选，无须设定绝对致死水平。

模型的主分类头和分类评价使用**dLFC<−1的强协同事件**，连续回归头输出预测dLFC。{range_text}`synergy_probability`对应前述强协同事件分数，`pred_dlfc`用于阅读预测协同方向与幅度，二者各有用途。

交付顺序综合单基因条件与协同分数。Top500中{tiers[0]}对两臂均未触发强单基因耗竭，{tiers[1]}对需补充单基因测量，{tiers[2]}对作为强单基因背景下的协同抑制扩展；各类独立清单位于`results/candidates`。

## 研究与验证衔接

当前结果为HT29筛选数据上学习得到的模型预测。核心配对属于外推应用，主测试和校准分别有{primary['positives']}个、{calibration['positives']}个强协同阳性；实际组合效应及HCT116-R/SW480-R中的表现通过合作方回传测量继续评价。详细数据覆盖与模型比较集中列在实验汇报中。

## 复现与存储

本GitHub仓库提供实现代码、配置、报告、图表、候选结果、原始与处理后数据、模型权重及逐次运行记录。网页可直接阅读报告和候选；完整克隆前安装Git和Git LFS：

```bash
git lfs install
git clone {REPOSITORY_URL}.git
cd MetaGene
git lfs pull
```

在仓库根目录运行`python scripts/write_sl_documents.py`，可根据已有精简结果重新生成阅读文档；网页阅读无需运行任何程序。

数据、权重和运行记录随仓库提供，大文件由Git LFS管理；交接直接使用GitHub，不生成交接压缩包。设备自动选择、指定设备以及修改模型后的独立运行目录见[复现说明](docs/REPRODUCIBILITY.md)。安装[项目依赖](requirements.txt)后按该说明运行；参数见[冻结配置](configs/experiment_sl.json)，当前结果位于`results/sl`。
'''
    design = f'''# {PROJECT_TITLE}

**研究设计｜总课题路径与干实验分工**

## 科学问题

围绕结直肠癌奥沙利铂耐药，研究代谢重编程与表观遗传重塑之间的组合脆弱性：能否通过AI主动学习优选代谢酶与组蛋白修饰酶配对，以双基因干预产生超出加性预期的协同抑制，并在后续验证中恢复奥沙利铂敏感性？

总课题串联**AI候选发现、聚焦筛选、耐药逆转验证和反馈学习**。MetaGene承载其中的干实验计算流程，当前阶段交付可供医学院衔接的排序基因组合。

实现代码、配置、报告、图表、数据、权重和运行结果统一在[GitHub仓库]({REPOSITORY_URL})交付，大文件使用Git LFS。阅读顺序与完整克隆步骤见[README](../README.md)；医学院可直接在网页查看报告、候选与交接说明。

## 整体研究路径与分工

| 阶段 | 负责方 | 工作与交付 |
|---|---|---|
| AI预测与候选发现 | 计算团队 | 多源数据整合、GNN建模、主动学习策略评价；交付Top500基因组合、短名单、逐对解释及回传接口 |
| 聚焦筛选与功能验证 | 医学院团队 | 承接基因组合、衔接实验设计，在耐药及匹配亲本模型中评价协同效应与奥沙利铂敏感性变化 |
| 反馈学习与下一轮优选 | 双方衔接，计算团队更新模型 | 汇总实测结果、评价当前预测、更新模型并提出下一轮候选 |

当前已有主动学习评价为{ac['rounds']}轮回顾性模拟；后续以医学院回传数据推进前瞻性迭代。

## 干实验预测目标

中心化LFC下，`dLFC = LFC_AB − LFC_A − LFC_B`。dLFC<0表示负向协同，关注组合相对于预期的额外抑制。本版训练以dLFC<−1为强协同分类事件，以连续dLFC辅助学习；候选的协同方向按预测dLFC是否小于0解读。

## 干实验四模块

| 模块 | 工作内容 | 产出 |
|---|---|---|
| 对象定义与多源表示 | UniProt催化功能、EpiFactors组蛋白修饰范围；DepMap、TCGA、CRC单细胞、通量模型、组蛋白峰及独立网络 | 基因目录、特征、网络与来源 |
| 协同信号预测 | 按重复及guide位置重建加性残差；交换对称GNN联合分类、回归与排序学习 | 预测dLFC、事件分数、留出与独立来源评价 |
| 耐药背景关联 | HT29/HTOXAR3及HCT116/SW480表达、药物筛选、单基因依赖和多组学背景 | 逐基因和逐配对的背景证据 |
| 主动学习与交接 | {ac['rounds']}轮回顾性查询、同预算随机对照、候选分层与多样化选择 | Top500→Top30→Priority8→3对代表组合、解释卡及回传接口 |

## 模型与计算实验设计

监督标签来自HT29 knockout，至少2个有效匹配重复。特征标准化在训练配对涉及的基因上拟合，外部网络保持固定。配对表示采用两臂嵌入的和、绝对差及逐元素积；每个候选使用同一预测器，实测单基因效应用于候选分层。

固定分割种子{cfg['partition_seed']}，主协议按{split_fractions}分配训练、验证、校准和测试。模型按验证AUPRC选择，测试集报告所有基线；另评价未见基因和未见锚点。测试报告固定后，用原训练+测试重新拟合最终模型，保留验证和校准集合。

强协同事件分数由{len(cfg['seeds'])}个模型集成并校准；dLFC区间基于校准残差，种子SD记录优化随机性。训练、分割和评价定义保存在[冻结配置](../configs/experiment_sl.json)。

## 干实验候选交付

{n_metabolic}个代谢酶与{n_histone}个直接组蛋白修饰酶组成{candidate_space}种配对，排除{previously_assayed}种既有设计后排序{candidate_universe}对。当前Top500的{direction_text}。候选先按单基因条件分层，再按协同事件分数、预测dLFC及配对键排序；Top30限制每个基因最多3次。

交付分别提供催化功能事实、组合机制假说、预测效应和耐药背景证据。{tiers[0]}对为单基因初筛兼容候选，{tiers[1]}对需补充单基因测量，{tiers[2]}对归入强单基因背景协同抑制扩展。

本版是HT29基础生长背景的计算预测，核心配对属于外推应用。合作方回传同背景的单A、单B、双干预及对照测量后，先检验当前预测，再更新下一轮模型。实际结果与覆盖信息见[实验汇报](../results/MetaGene实验汇报.md)。
'''
    handoff = f'''# 干实验交付说明：潜在合成致死基因组合

所属课题：**{PROJECT_TITLE}**。

本说明用于计算团队向医学院团队移交基因组合。计算团队交付排序候选、预测结果、生物学解释及结果回传接口；医学院团队衔接聚焦筛选和代表组合验证，评价协同效应及奥沙利铂敏感性变化。实测结果回传后用于模型评价和下一轮主动学习。

统一交付入口为[MetaGene GitHub仓库]({REPOSITORY_URL})，网页阅读无需运行模型。建议依次查看[README](../README.md)、[实验汇报](../results/MetaGene实验汇报.md)、[Top500](../results/candidates/top500.csv)／[Top30](../results/candidates/top30.csv)／[证据卡](../results/candidates/candidate_evidence_cards.md)，再按本说明使用[回传模板](../results/handoff/measured_feedback_template.csv)。完整数据、权重和运行记录按README中的Git LFS步骤获取。

## 本次交付

版本：{version}。[Top500](../results/candidates/top500.csv)中{negative_count}/{n_top}对预测dLFC<0，范围{min(delta):.3f}至{max(delta):.3f}，即模型预测双干预产生低于单基因加性预期的适应性变化。每对提供分数、单基因条件、背景证据、功能来源及生物学解释。

首批讨论使用[Top30](../results/candidates/top30.csv)、[逐对证据卡](../results/candidates/candidate_evidence_cards.md)与[3对代表组合](../results/handoff/first_round_representatives.csv)：{representative_names}。

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

Top500中{tiers[0]}对两臂均未触发强单基因耗竭，{tiers[1]}对需补充单基因测量，{tiers[2]}对作为强单基因背景的协同抑制扩展。独立清单为[单基因初筛兼容](../results/candidates/single_compatible_candidates.csv)、[测量待补充](../results/candidates/single_evidence_missing_candidates.csv)和[强单基因扩展](../results/candidates/strong_single_synergy_extension.csv)。功能来源用于解释两臂各自的催化角色，背景数据用于讨论细胞模型的相关性。

本次交付为模型预测候选；合作方的双干预测量用于确认配对效应，匹配亲本与耐药背景用于评价耐药相关性。

## 容量安排与结果回传

[500配对靶点表](../results/handoff/500_pair_targets.csv)定义基因级对象；[1,000容量草案](../results/handoff/1000_slot_capacity_draft.csv)安排{slot_counts.get('primary_candidate_capacity',0)}个首设计、{slot_counts.get('secondary_design_capacity',0)}个次设计和{slot_counts.get('control_or_single_capacity_reservation',0)}个控制占位，guide设计由合作方衔接。

[回传模板](../results/handoff/measured_feedback_template.csv)记录配对键、两臂基因、细胞模型、亲本/耐药状态、药物条件、重复、单A/单B/双干预及对照LFC、标尺、质控和数据来源。未中心化LFC按每个重复计算`double − single_A − single_B + control`；已中心化数据的对照项为0。保留生物重复层级，用重复测量评价统计支持。

新数据首先用于评价本版预测，随后纳入下一轮学习。GitHub统一管理代码、配置、报告、图表、候选结果、数据、模型权重和运行记录；大文件使用Git LFS，交接直接使用仓库。
'''
    (ROOT/'README.md').write_text(readme)
    (ROOT/'results/MetaGene实验汇报.md').write_text(report)
    (ROOT/'docs/RESEARCH_DESIGN.md').write_text(design)
    (ROOT/'docs/CANDIDATE_HANDOFF.md').write_text(handoff)
    cards=render_candidate_cards(short)
    (candidates/'candidate_evidence_cards.md').write_text(cards)
    if (OUT/'candidates').is_dir():
        (OUT/'candidates/candidate_evidence_cards.md').write_text(cards)
        mirrored=report.replace('](sl/','](').replace('](../data/','](../../data/').replace('](../configs/','](../../configs/')
        mirrored=mirrored.replace('](handoff/','](../handoff/')
        (OUT/'report.md').write_text(mirrored)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
