#!/usr/bin/env python3
"""Assemble a medical-research handoff from current candidate predictions.

This stage reads predictions, annotations, disease evidence and independent
pair measurements. It does not train, select or evaluate a predictive model.
Run after finalize_revision.py --phase predict and the GSE76092 limma fit:
    python scripts/rank_candidates.py
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = 'resistant_vs_parental_untreated'
SECONDARY = ['parental_oxa_response', 'resistant_oxa_response', 'state_by_oxa_interaction']
CLASS_ORDER = ['core_metabolism_histone', 'epigenetic_extension', 'signaling_extension']
ANNOTATION_COLUMNS = [
    'is_strict_metabolic_enzyme', 'is_extended_metabolic_signaling',
    'is_chromatin_regulator', 'is_core_chromatin_factor', 'is_direct_histone_enzyme', 'histone_scope_reason', 'metabolic_class',
    'metabolic_subclass', 'metabolic_reason', 'epigenetic_class',
    'chromatin_evidence_tier', 'chromatin_mechanism_class', 'chromatin_tier_reason',
    'chromatin_annotation_caveat', 'chromatin_tier_source',
    'reviewed_accessions', 'ec_numbers', 'reviewed_source_urls',
    'epifactors_functions', 'epifactors_function_pmids',
]
PRIMARY_COLUMNS = [
    'log2fc', 'ci95_low', 'ci95_high', 'average_log2_intensity', 'moderated_t',
    'pvalue', 'padj_bh', 'padj_by', 'significant_bh_005', 'significant_by_005',
    'high_effect_bh', 'HT29_untreated_mean_log2_intensity',
    'HT29_untreated_sd_log2_intensity', 'HT29_untreated_relative_intensity_percentile',
    'HTOXAR3_untreated_mean_log2_intensity', 'HTOXAR3_untreated_sd_log2_intensity',
    'HTOXAR3_untreated_relative_intensity_percentile',
]
SECONDARY_COLUMNS = ['log2fc', 'ci95_low', 'ci95_high', 'padj_bh', 'padj_by']
PRIORITY_SORT = [
    'scaled_single_unknown_count', 'strong_single_count', 'raw_joint_unfavorable',
    'DepMap_strong_single_count', 'priority_band_order', 'seed_sd_band_order',
    'target_resistance_association_count', 'gse76092_high_effect_bh_gene_count',
    'gse76092_bh005_gene_count', 'pred_gi', 'pair_key',
]
PRIORITY_ASCENDING = [True, True, True, True, True, True, False, False, False, True, True]
GEO = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='


def diversified(frame, n, cap):
    """Take the first n rows in an explicit order with a cap across both arms."""
    counts = Counter()
    selected = []
    for index, row in frame.iterrows():
        if counts[row.gene_a] < cap and counts[row.gene_b] < cap:
            selected.append(index)
            counts.update([row.gene_a, row.gene_b])
        if len(selected) == n:
            break
    return frame.loc[selected].copy()


def representative_selection(frame, n):
    """Earliest feasible rank tuple with no repeated gene (small Priority8)."""
    for indices in combinations(frame.index, n):
        candidate = frame.loc[list(indices)]
        if pd.concat([candidate.gene_a, candidate.gene_b]).nunique() == 2*n:
            return candidate.copy()
    raise ValueError(f'Priority shortlist cannot supply {n} disjoint representative pairs')


def classify_pair(a, b, annotation):
    """Check both orientations and prefer the narrowest eligible mechanism."""
    orientations = [(a, b), (b, a)]
    for category, metabolic_flag, chromatin_flag in [
        ('core_metabolism_histone', 'is_strict_metabolic_enzyme', 'is_direct_histone_enzyme'),
        ('epigenetic_extension', 'is_strict_metabolic_enzyme', 'is_chromatin_regulator'),
        ('signaling_extension', 'is_extended_metabolic_signaling', 'is_chromatin_regulator'),
    ]:
        eligible = [(m, e) for m, e in orientations
                    if annotation.at[m, metabolic_flag] == 1 and annotation.at[e, chromatin_flag] == 1]
        if eligible:
            m, e = eligible[0]
            return category, m, e, ';'.join(f'{x}|{y}' for x, y in eligible)
    return 'outside_declared_candidate_scope', '', '', ''


def attach_gene_table(frame, gene_table):
    columns = {}
    for arm in ['a', 'b']:
        for column in gene_table.columns:
            columns[f'gene_{arm}_{column}'] = frame[f'gene_{arm}'].map(gene_table[column])
    return pd.concat([frame, pd.DataFrame(columns, index=frame.index)], axis=1)


def attach_pathway_evidence(frame, disease_dir):
    """Retain every tested member association for the primary contrast."""
    camera = pd.read_csv(disease_dir / 'reactome_camera.csv')
    camera = camera[camera.contrast.eq(PRIMARY)].rename(columns={
        'NGenes': 'camera_n_genes', 'Direction': 'camera_direction',
        'PValue': 'camera_pvalue', 'FDR': 'camera_fdr', 'method': 'camera_method',
    })
    ora = pd.read_csv(disease_dir / 'reactome_ora.csv')
    ora = ora[ora.contrast.eq(PRIMARY)]
    path_stats = camera.copy()
    ora_columns = ['pathway_id', 'n_pathway', 'n_selected', 'n_universe', 'n_overlap',
                   'pvalue', 'padj_bh', 'overlap_genes', 'selection_rule']
    for direction in ['up', 'down']:
        block = ora.loc[ora.direction.eq(direction), ora_columns]
        block = block.rename(columns={col: f'ora_{direction}_{col}' for col in ora_columns if col != 'pathway_id'})
        path_stats = path_stats.merge(block, on='pathway_id', how='left')
    membership = pd.read_csv(disease_dir / 'reactome_tested_membership.csv')
    member = membership[membership.gene.isin(set(frame.gene_a) | set(frame.gene_b))]
    gene_paths = member[['gene', 'pathway_id']].merge(path_stats, on='pathway_id')
    for direction in ['up', 'down']:
        gene_paths[f'gene_in_ora_{direction}_selected_set'] = [
            gene in overlap.split(';') for gene, overlap in
            zip(gene_paths.gene, gene_paths[f'ora_{direction}_overlap_genes'].fillna(''))]
    gene_paths = gene_paths.sort_values(['gene', 'camera_fdr', 'camera_n_genes', 'pathway_id'])
    pair_arms = frame[['pair_key', 'gene_a', 'gene_b']].melt(
        id_vars='pair_key', value_vars=['gene_a', 'gene_b'], var_name='arm', value_name='gene')
    pair_paths = pair_arms.merge(gene_paths, on='gene')
    pair_paths = pair_paths.merge(frame[['pair_key', 'pair_class']], on='pair_key')
    shared_ids = pair_paths.groupby(['pair_key', 'pathway_id']).size().reset_index(name='member_arms')
    shared = shared_ids[shared_ids.member_arms.eq(2)].merge(path_stats, on='pathway_id')
    shared = shared.merge(frame[['pair_key', 'gene_a', 'gene_b', 'pair_class']], on='pair_key')
    shared = shared.sort_values(['pair_key', 'camera_fdr', 'camera_n_genes', 'pathway_id'])
    genes = gene_paths.groupby('gene')
    best = gene_paths.drop_duplicates('gene').set_index('gene')
    summary = pd.DataFrame({'primary_camera_pathway_count': genes.size()})
    summary['primary_camera_fdr005_pathway_count'] = gene_paths.camera_fdr.lt(0.05).groupby(gene_paths.gene).sum()
    summary['primary_camera_min_fdr'] = best.camera_fdr
    summary['primary_camera_top_pathway_id'] = best.pathway_id
    summary['primary_camera_top_pathway_name'] = best.pathway_name
    summary['primary_camera_top_pathway_direction'] = best.camera_direction
    for direction in ['up', 'down']:
        summary[f'primary_ora_{direction}_fdr005_pathway_count'] = gene_paths[f'ora_{direction}_padj_bh'].lt(0.05).groupby(gene_paths.gene).sum()
    frame = attach_gene_table(frame, summary)
    for arm in ['a', 'b']:
        for column in summary.columns[summary.columns.str.endswith('_count')]:
            key = f'gene_{arm}_{column}'
            frame[key] = frame[key].fillna(0).astype(int)
    frame['shared_tested_pathway_count'] = frame.pair_key.map(shared.groupby('pair_key').size()).fillna(0).astype(int)
    frame['shared_camera_fdr005_pathway_count'] = frame.pair_key.map(
        shared.camera_fdr.lt(0.05).groupby(shared.pair_key).sum()).fillna(0).astype(int)
    frame['shared_camera_min_fdr'] = frame.pair_key.map(shared.groupby('pair_key').camera_fdr.min())
    return frame, gene_paths, pair_paths, shared


def pathway_card_lines(pathways):
    if not len(pathways):
        return ['本次主对比的已检验 Reactome 集合中没有成员关联记录。']
    significant = int(pathways.camera_fdr.lt(0.05).sum())
    lines = [f'共有 {len(pathways)} 条已检验成员通路，其中 {significant} 条 CAMERA FDR<0.05；以下按 CAMERA FDR 最小取3条，包含未达阈值结果。']
    for path in pathways.head(3).itertuples(index=False):
        lines.append(
            f'- [{path.pathway_name}](https://reactome.org/content/detail/{path.pathway_id})：'
            f'CAMERA Direction={path.camera_direction}，FDR={fmt(path.camera_fdr)}，nGenes={int(path.camera_n_genes)}；'
            f'ORA上调集合FDR={fmt(path.ora_up_padj_bh)}、下调集合FDR={fmt(path.ora_down_padj_bh)}。')
    return lines


def single_strata(frame, cutoff):
    out = frame.copy()
    out['strong_single_cutoff'] = cutoff
    covered = []
    for arm in ['a', 'b']:
        scaled = out[f'gene_{arm}_colo1_single_control_scaled']
        covered.append(scaled.notna())
        out[f'gene_{arm}_strong_single'] = scaled.le(cutoff).astype('boolean').where(scaled.notna())
    out['scaled_single_covered_gene_count'] = covered[0].astype(int) + covered[1].astype(int)
    out['scaled_single_unknown_count'] = 2 - out.scaled_single_covered_gene_count
    out['known_strong_single_count'] = out[['gene_a_strong_single', 'gene_b_strong_single']].fillna(False).sum(axis=1).astype(int)
    complete = out.scaled_single_unknown_count.eq(0)
    out['strong_single_count'] = out.known_strong_single_count.where(complete).astype('Int64')
    strong = out.known_strong_single_count.gt(0) & complete
    no_strong = out.known_strong_single_count.eq(0) & complete
    out['single_effect_stratum'] = np.select([
        ~complete, out.strong_single_count.fillna(-1).astype(int).eq(1), out.strong_single_count.fillna(-1).astype(int).eq(2),
    ], ['unknown_control_scaled_single', 'one_strong_control_scaled_single',
        'two_strong_control_scaled_singles'], default='no_strong_control_scaled_single_at_cutoff')
    negative = out.pred_gi.lt(0)
    upper_negative = out.pi90_high.lt(0)
    out['priority_band_order'] = np.select([
        no_strong & negative & upper_negative,
        no_strong & negative & ~upper_negative,
        strong & negative & upper_negative,
        strong & negative & ~upper_negative,
        no_strong & ~negative,
        strong & ~negative,
        ~complete & negative,
    ], [0, 1, 2, 3, 4, 5, 6], default=7)
    out['priority_band'] = out.priority_band_order.map({
        0: 'A_negative_GI_interval_below_zero_no_strong_single',
        1: 'B_exploratory_negative_GI_interval_includes_zero_no_strong_single',
        2: 'C_negative_GI_interval_below_zero_strong_single',
        3: 'D_exploratory_negative_GI_interval_includes_zero_strong_single',
        4: 'E_exploratory_nonnegative_GI_no_strong_single',
        5: 'F_exploratory_nonnegative_GI_strong_single',
        6: 'U1_exploratory_negative_GI_unknown_control_scaled_single',
        7: 'U2_exploratory_nonnegative_GI_unknown_control_scaled_single',
    })
    out['interval_exploration_group'] = ~out.pi90_high.lt(0)
    out['single_effect_interpretation'] = np.select([
        ~complete, strong,
    ], ['至少一臂作者control-scaled单基因缺测；强single状态未知，需合作方确认',
        '作者control-scaled单基因达到较强耗竭阈值；组合判别需确认干预强度和动态范围'],
        default='两臂作者control-scaled值均已测且高于当前强耗竭阈值；此标签使用负控0、正控−1标尺')
    return out


def priority_selection(top, cutoff, n, cap):
    ranked = single_strata(top, cutoff).sort_values(PRIORITY_SORT, ascending=PRIORITY_ASCENDING)
    picked = diversified(ranked, n, cap)
    picked['priority8_rank'] = np.arange(1, len(picked) + 1)
    return picked


def fmt(value, digits=3):
    return '未覆盖/未提供' if pd.isna(value) else f'{float(value):.{digits}g}'


def text_value(value):
    return '未提供' if pd.isna(value) or str(value) == '' else str(value)


def source_links(value):
    return '；'.join(f'[来源{i+1}]({url})' for i, url in enumerate(str(value).split(';')) if url)


def evidence_card(row, external, prediction_metadata, cutoff, gene_path_groups, shared_paths, empty_paths, source_proxy_summary):
    """One independent, self-contained evidence card per shortlisted pair."""
    pair_label = f"{row['metabolic_gene']} + {row['chromatin_gene']}"
    lines = [
        f"# {int(row['top30_rank']):02d}. {pair_label}",
        '',
        '**交接状态：待实验验证的组合候选。** 当前 HT29 组合筛选没有该组合的测量；两臂具有同来源单基因测量。',
        f"交接证据层：{row['handoff_evidence_tier']}；需要交接探索：{row['handoff_requires_exploration']}。{row['handoff_evidence_reason']}",
        '',
        f"GI 主长名单排名 {int(row['rank_gi_within_class'])}；多样化 Top30 排名 {int(row['top30_rank'])}。"
        f"Priority8：{int(row['priority8_rank']) if pd.notna(row['priority8_rank']) else '本轮未选入'}。",
        '',
        '## 计算预测与已测单基因',
        '',
        f"- 预测加性 GI：**{fmt(row['pred_gi'])}**；上游 90% GI 预测区间 [{fmt(row['pi90_low'])}, {fmt(row['pi90_high'])}]；跨随机种子 SD {fmt(row['between_seed_sd'])}。",
        f"- 预测联合效应 {fmt(row['pred_joint'])} = 已测单基因和 {fmt(row['single_effect_sum'])} + 预测 GI；这是同一预测的代数重构。",
        f"- 原始单基因 LFC：{row['gene_a']}={fmt(row['single_a'])}；{row['gene_b']}={fmt(row['single_b'])}，用于GI和joint的源标尺计算。",
        f"- 作者control-scaled单基因：{row['gene_a']}={fmt(row['gene_a_colo1_single_control_scaled'])}；{row['gene_b']}={fmt(row['gene_b_colo1_single_control_scaled'])}。"
        f"负控0、正控−1；强耗竭阈值 ≤ {cutoff:g}。覆盖{int(row['scaled_single_covered_gene_count'])}/2臂，已观测达到阈值{int(row['known_strong_single_count'])}臂，缺测{int(row['scaled_single_unknown_count'])}臂。",
        f"- 单基因层：{row['single_effect_stratum']}。{row['single_effect_interpretation']}。",
        f"- 不确定性层：{row['interval_sign_class']}；seed SD 位于核心长名单的 {row['seed_sd_band']}。GI区间/single数学分层：{row['priority_band']}；此A–U编码不是生物学证据等级。",
        f"- 模型 {row['model']}，{int(row['n_models'])} 个拟合；设计范围 {row['design_scope']}。",
        f"- 区间来源：{prediction_metadata['interval_definition']}。适用范围：{prediction_metadata['interval_scope']}。",
        f"- 种子差异范围：{prediction_metadata['seed_sd_scope']}。",
        '',
        f"GI proxy是联合生长效应相对于两个单基因效应之和的偏离，其零点不是生物学零效应或杀伤阈值。该HT29来源{source_proxy_summary['n_measured_pairs']:,}个已测组合中，{source_proxy_summary['fraction_below_zero']:.2%}的proxy小于0，均值{source_proxy_summary['mean']:.3f}、中位数{source_proxy_summary['median']:.3f}。因此区间上界<0仅表示相对这个proxy零点的数学分层，不构成确认或高置信的生物学合成致死结论。",
        '原始joint是相对plasmid参考的文库丰度LFC，受到文库组成归一化影响；它的0也不是无杀伤/绝对存活阈值。负GI可以同时对应正joint。Top500按GI排序；Top30与首轮代表组合另行审查raw joint、单基因效应及耐药关联。未建立生物学null校准。',
        '',
        '## 两臂功能与耐药表达证据',
        '',
        'GSE76092 主比较为未处理 HTOXAR3 耐药株减去未处理 HT29 亲本，独立 3 对 3。正向表示耐药株相对上调，负向表示相对下调；两种方向均描述关联。',
        '',
    ]
    for arm in ['a', 'b']:
        gene = row[f'gene_{arm}']
        role = '严格代谢酶' if gene == row['metabolic_gene'] else '直接组蛋白修饰酶'
        reason = row[f'gene_{arm}_metabolic_reason'] if role == '严格代谢酶' else row[f'gene_{arm}_histone_scope_reason'] + ' ' + str(row[f'gene_{arm}_chromatin_tier_reason'])
        prefix = f'gene_{arm}_GSE76092_'
        lines.extend([
            f"### {gene}（{role}）", '',
            f"- 机制：{text_value(reason)}",
            f"- UniProt {text_value(row[f'gene_{arm}_reviewed_accessions'])}；EC {text_value(row[f'gene_{arm}_ec_numbers'])}；"
            f"{source_links(row[f'gene_{arm}_reviewed_source_urls'])}。",
            f"- EpiFactors 功能：{text_value(row[f'gene_{arm}_epifactors_functions'])}；功能文献 PMID {text_value(row[f'gene_{arm}_epifactors_function_pmids'])}。",
            f"- 作者control-scaled single {fmt(row[f'gene_{arm}_colo1_single_control_scaled'])}，基于{fmt(row[f'gene_{arm}_n_guide_pairs_in_single_estimate'], 2)}个gene-control guide组合；"
            f"独立Project SCORE HT29 LFC {fmt(row[f'gene_{arm}_score_ht29_lfc'])}，源文件相对rank {fmt(row[f'gene_{arm}_score_ht29_relative_rank'])}。"
            'SCORE原始效应与相对rank作为独立单基因证据展示，未转换为正式essentiality概率或肿瘤选择性。',
            f"- DepMap亲本细胞gene effect：HT29 {fmt(row[f'gene_{arm}_DepMap_HT29_gene_effect'])}；HCT116 {fmt(row[f'gene_{arm}_DepMap_HCT116_gene_effect'])}；SW480 {fmt(row[f'gene_{arm}_DepMap_SW480_gene_effect'])}。该单基因标尺不能替代耐药株或双基因结果。",
            f"- TCGA-COAD肿瘤中位表达 {fmt(row[f'gene_{arm}_TCGA_tumor_median_log2_RSEM'])}；GSE132465肿瘤组织上皮细胞的患者等权平均检测率 {fmt(row[f'gene_{arm}_sc_tumor_epithelial_mean_detection'])}。",
            f"- HCT116/OxR作者模型通量：映射反应数 {fmt(row[f'gene_{arm}_flux_HCT116_mapped_reaction_count'], 0)}；生长率归一化后log1p绝对通量界的基因关联均值，亲本 {fmt(row[f'gene_{arm}_flux_HCT116_Control_mean_log1p_abs_bound'])}、OxR {fmt(row[f'gene_{arm}_flux_HCT116_OxR_mean_log1p_abs_bound'])}。这是反应范围的描述性汇总，不是实测细胞内通量、置信区间或该基因敲除效应。",
            '- GSE138654组蛋白图谱：各修饰的启动子峰最大富集log1p，按亲本/OxR顺序：' + '；'.join(
                f"{mark} {fmt(row[f'gene_{arm}_ChIP_{mark}_HCT116_log1p_max_fold_enrichment'])}/{fmt(row[f'gene_{arm}_ChIP_{mark}_HCT116OxR_log1p_max_fold_enrichment'])}"
                for mark in ['H3K4me3','H3K27ac','H3K27me3']) + '。没有重复差异结合检验；0表示未在已映射启动子窗口检出作者调用峰，不证明不存在修饰。不同研究的OxR株不假定相同。',
            f"- GSE299428 HCT116高耐药对亲本：log2FC {fmt(row[f'gene_{arm}_GSE299428_HCT116_HD_vs_PAR_log2fc'])}，FDR {fmt(row[f'gene_{arm}_GSE299428_HCT116_HD_vs_PAR_padj'])}；SW480高耐药对亲本：log2FC {fmt(row[f'gene_{arm}_GSE299428_SW480_HD_vs_PAR_log2fc'])}，FDR {fmt(row[f'gene_{arm}_GSE299428_SW480_HD_vs_PAR_padj'])}。",
            f"- GSE76092 log2FC {fmt(row[prefix+'log2fc'])}，95%均值对比区间 [{fmt(row[prefix+'ci95_low'])}, {fmt(row[prefix+'ci95_high'])}]，"
            f"BH FDR {fmt(row[prefix+'padj_bh'])}，BY FDR {fmt(row[prefix+'padj_by'])}；标签 {row[f'gene_{arm}_resistance_expression_state']}。",
            f"- 亲本相对芯片强度均值 {fmt(row[prefix+'HT29_untreated_mean_log2_intensity'])}，强度分位数 {fmt(row[prefix+'HT29_untreated_relative_intensity_percentile'])}；"
            f"耐药株均值 {fmt(row[prefix+'HTOXAR3_untreated_mean_log2_intensity'])}，强度分位数 {fmt(row[prefix+'HTOXAR3_untreated_relative_intensity_percentile'])}。"
            '相对强度用于比较，绝对表达量、蛋白丰度与可检测性需合作方确认。',
            f"- 其他表达：GSE299428 HCT116 MD log2FC {fmt(row[f'gene_{arm}_GSE299428_HCT116_MD_vs_PAR_log2fc'])} / FDR {fmt(row[f'gene_{arm}_GSE299428_HCT116_MD_vs_PAR_padj'])}；"
            f"GSE289766 HCT116 log2FC {fmt(row[f'gene_{arm}_GSE289766_HCT116_OX_vs_WT_log2fc'])}；"
            f"GSE42387 HCT116 log2FC {fmt(row[f'gene_{arm}_GSE42387_HCT116_OX_vs_PAR_log2fc'])}、HT29 log2FC {fmt(row[f'gene_{arm}_GSE42387_HT29_OX_vs_PAR_log2fc'])}。"
            '后两项来源只提供此处使用的效应量。',
            f"- HCT116 跨研究方向：{text_value(row[f'gene_{arm}_hct116_direction'])}；覆盖 {fmt(row[f'gene_{arm}_hct116_studies_measured'], 2)} 项，"
            f"方向支持比例 {fmt(row[f'gene_{arm}_hct116_directional_support_fraction'])}。",
            f"- GSE254136 单基因奥沙利铂筛选 log2FC {fmt(row[f'gene_{arm}_GSE254136_screen_log2fc'])}；"
            f"有效 guide 数 {fmt(row[f'gene_{arm}_GSE254136_n_usable_guides'], 2)}；guide 耗竭比例 {fmt(row[f'gene_{arm}_GSE254136_guide_depletion_fraction'])}；"
            f"处理重复耗竭比例 {fmt(row[f'gene_{arm}_GSE254136_treated_replicate_depletion_fraction'])}。"
            '负向表示药物条件下 guide 相对耗竭；单个 PBS 对照限制了效应推断，此项是单基因药物关联线索。',
            '',
        ])
        lines.extend(['**主耐药对比的通路成员证据**', ''])
        lines.extend(pathway_card_lines(gene_path_groups.get(gene, empty_paths)))
        lines.extend(['', '成员关系说明该基因属于该集合；CAMERA/ORA描述集合级耐药关联。集合方向不指定每个成员的方向，也不建立该基因的因果作用或双扰动协同。', ''])
    lines.extend([
        f"表达来源：[GSE76092]({GEO}GSE76092)、[GSE299428]({GEO}GSE299428)、[GSE289766]({GEO}GSE289766)、[GSE42387]({GEO}GSE42387)；"
        f"药筛来源：[GSE254136]({GEO}GSE254136)。",
        '',
        '主表另保留两个细胞状态各自的急性用药反应和交互项。HT29 使用 10 μM、HTOXAR3 使用 30 μM 奥沙利铂；交互项具有剂量混杂，不能解释为同剂量耐药逆转。',
        '',
        '## 两臂共享的主对比通路', '',
    ])
    lines.extend(pathway_card_lines(shared_paths))
    lines.extend([
        '', '共享通路表示在所检验集合中的共同成员身份。其集合级统计量与每个基因或该 pair 的扰动效应分别解释。完整成员关联和未达阈值结果保留在通路证据CSV。', '',
        '## 独立研究中的同一组合', '',
    ])
    if len(external):
        for item in external.to_dict('records'):
            result = '达到来源效应量 hit 阈值' if item['study_hit'] else '实测且未达到来源 hit 阈值'
            lines.append(
                f"- {item['study']} / {item['cell_line']}：{result}；作者 GI {fmt(item['author_gi'])}，"
                f"作者联合效应 {fmt(item['joint_effect'])}；阈值 {item['hit_definition']}；"
                f"{item['perturbation_type']}；[原研究](https://doi.org/{item['source_doi']})。"
                f"在任一项目训练来源出现同 pair：{item['pair_seen_in_any_training_study']}。")
        lines.extend(['', '外部 hit 保留来源阈值；实测未过阈值表示该实验的非 hit。跨细胞系结果单独呈现。当前外部库采用基础生长条件，未提供奥沙利铂下的双扰动或获得性耐药验证。'])
    else:
        lines.append('当前独立外部组合库未覆盖这一 pair；其外部组合状态为未知。')
    lines.extend([
        '', '## 合作方待验证的问题与反馈字段', '',
        f"- 操作性信息：{row['experimental_operability_status']}。需反馈可用细胞模型及亲本/耐药对应关系、目标转录本和蛋白表达、可用干预工具、干预效率及选择性证据。",
        '- 组合判别问题：两臂单基因效应与测量动态范围是否允许识别额外组合效应；可提供哪些单基因/组合对照和独立证据来检验相对加性预期的偏离。',
        '- 耐药适用问题：所观察的效应属于基础生长抑制、双扰动交互还是药物敏化；亲本/耐药模型、药物暴露和读数是否具有可比性。',
        '- 已测证据：本来源两臂单基因效应；所列耐药表达和单基因药筛；若上节有记录则包括独立来源同 pair 测量。',
        '- 项目推测：当前 HT29 中的双扰动 GI、组合机制、耐药背景适用性和药物响应变化。最终交接结论为待验证候选。',
        '',
    ])
    return '\n'.join(lines)


def write_handoff(out, top500, representatives):
    """Gene-pair and result-schema handoff, with unfilled design fields explicit."""
    handoff = out.parent / 'handoff'
    handoff.mkdir(exist_ok=True)
    pairs = top500[['top500_rank', 'pair_key', 'metabolic_gene', 'chromatin_gene', 'pred_gi',
                    'observed_HT29_single_count', 'handoff_evidence_tier']].copy()
    pairs.insert(0, 'candidate_id', [f'MG{i:04d}' for i in range(1, len(pairs)+1)])
    pairs['target_resistant_models'] = 'HCT116-R;SW480-R'
    pairs['matched_parental_models'] = 'HCT116;SW480'
    pairs['partner_context_confirmation'] = 'pending: exact acquired-resistance model identity'
    pairs['pair_validation_status'] = 'unvalidated_hypothesis'
    pairs.to_csv(handoff / 'mini_pool_gene_pairs.csv', index=False)
    # These are capacity reservations, not reagent designs or a synthesis order.
    slots = []
    for i, row in pairs.iterrows():
        slots.append({'slot_id': f'S{len(slots)+1:04d}', 'slot_type': 'primary_candidate',
                      'candidate_id': row.candidate_id, 'pair_key': row.pair_key,
                      'guide_a_id': '', 'guide_b_id': '', 'design_status': 'awaiting_partner_guide_assignment'})
    for row in pairs.head(450).itertuples():
        slots.append({'slot_id': f'S{len(slots)+1:04d}', 'slot_type': 'secondary_design_reservation',
                      'candidate_id': row.candidate_id, 'pair_key': row.pair_key,
                      'guide_a_id': '', 'guide_b_id': '', 'design_status': 'awaiting_partner_guide_assignment'})
    for _ in range(50):
        slots.append({'slot_id': f'S{len(slots)+1:04d}', 'slot_type': 'control_reservation',
                      'candidate_id': '', 'pair_key': '', 'guide_a_id': '', 'guide_b_id': '',
                      'design_status': 'awaiting_partner_control_definition'})
    pd.DataFrame(slots).to_csv(handoff / 'library_capacity_plan.csv', index=False)
    fields = ['record_id', 'candidate_id', 'gene_a', 'gene_b', 'cell_model', 'parental_or_resistant',
              'drug_condition_label', 'replicate_id', 'single_a_effect', 'single_b_effect', 'joint_effect',
              'effect_scale', 'control_normalization', 'quality_pass', 'source_file', 'measured_at']
    pd.DataFrame(columns=fields).to_csv(handoff / 'measured_feedback_template.csv', index=False)
    representatives[['representative_rank', 'pair_key', 'metabolic_gene', 'chromatin_gene', 'pred_gi',
                     'pred_joint', 'single_a', 'single_b', 'nomination_status']].to_csv(handoff / 'first_round_representatives.csv', index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    args = parser.parse_args()
    root = args.root
    out = root / 'results/candidates'
    out.mkdir(parents=True, exist_ok=True)
    config = json.loads((root / 'configs/experiment.json').read_text())
    policy = config['candidate_policy']
    cutoff = policy['strong_single_cutoff']
    cap = policy['maximum_pairs_per_gene']
    top_n = policy['shortlist_size']
    priority_n = policy['priority_size']
    prediction_metadata = json.loads((out / 'prediction_metadata.json').read_text())
    frame = pd.read_csv(out / 'predicted_pairs.csv')
    annotation = pd.read_csv(root / 'data/evidence/genes_annotation.csv').set_index('gene')
    disease_dir = root / 'results/disease'
    primary = pd.read_csv(disease_dir / f'{PRIMARY}_gene_effects.csv').set_index('gene')
    effects = pd.read_csv(disease_dir / 'GSE76092_gene_effects.csv')
    old_evidence = pd.read_csv(root / 'results/drug_evidence/disease_evidence.csv', low_memory=False).set_index('gene')
    external = pd.read_csv(root / 'data/evidence/external/external_pairs.csv')
    observed = pd.read_csv(root / 'data/processed/pairs.csv')
    observed_ht29_proxy = observed.loc[observed.dataset.eq('ht29_colo1'), 'additive_gi_proxy'].dropna()
    source_proxy_summary = {
        'source': 'data/processed/pairs.csv / ht29_colo1 / additive_gi_proxy',
        'n_measured_pairs': int(len(observed_ht29_proxy)),
        'fraction_below_zero': float(observed_ht29_proxy.lt(0).mean()),
        'mean': float(observed_ht29_proxy.mean()),
        'median': float(observed_ht29_proxy.median()),
        'use': 'Descriptive source distribution only; no recentering, ranking change or biological null calibration.',
    }
    singles = pd.read_csv(root / 'data/processed/singles.csv')
    singles = singles[singles.dataset.eq('ht29_colo1')].set_index('gene')
    single_evidence = pd.read_csv(root / 'data/evidence/single_gene_evidence.csv').set_index('gene')

    roles = [classify_pair(a, b, annotation) for a, b in frame[['gene_a', 'gene_b']].to_numpy()]
    frame[['pair_class', 'metabolic_gene', 'chromatin_gene', 'eligible_orientations']] = roles
    measured = observed.loc[observed.dataset.eq('ht29_colo1'), ['gene_a', 'gene_b']]
    measured_keys = {'|'.join(sorted(pair)) for pair in measured.to_numpy()}
    training_keys = {'|'.join(sorted(pair)) for pair in observed[['gene_a', 'gene_b']].to_numpy()}
    frame['measured_current_ht29_pair'] = frame.pair_key.isin(measured_keys)
    frame['pair_seen_in_any_training_study'] = frame.pair_key.isin(training_keys)
    mapping = pd.read_csv(root / 'data/annotations/symbol_mapping.csv', keep_default_na=False).set_index('source_symbol').gene
    source_design = pd.read_csv(root / 'data/pairs_sources/HT29_annotations.csv')
    raw_anchors = source_design.loc[source_design.MyNote.eq('AnchorSingletons'), 'Gene2']
    anchors = set(raw_anchors.map(mapping).fillna(raw_anchors))
    frame['anchor_gene_count'] = frame.gene_a.isin(anchors).astype(int) + frame.gene_b.isin(anchors).astype(int)
    frame['design_scope'] = np.where(frame.anchor_gene_count.gt(0), 'contains_source_screen_anchor', 'both_nonanchor_extrapolation')

    gene_table = annotation[ANNOTATION_COLUMNS].join(old_evidence, how='outer')
    gene_table = gene_table.join(primary[PRIMARY_COLUMNS].add_prefix('GSE76092_'), how='outer')
    for contrast in SECONDARY:
        part = effects.loc[effects.contrast.eq(contrast)].set_index('gene')
        gene_table = gene_table.join(part[SECONDARY_COLUMNS].add_prefix(f'GSE76092_{contrast}_'), how='outer')
    gene_table = gene_table.join(singles[['single_definition', 'perturbation_type']].add_prefix('ht29_single_'), how='outer')
    gene_table = gene_table.join(single_evidence[[
        'colo1_single_control_scaled', 'n_guide_pairs_in_single_estimate',
        'score_ht29_lfc', 'score_ht29_relative_rank']], how='outer')
    gene_table = gene_table.join(pd.read_csv(root / 'data/evidence/integrated_gene_evidence.csv').set_index('gene'), how='outer')
    frame = attach_gene_table(frame, gene_table)
    frame, gene_paths, pair_paths, shared_paths = attach_pathway_evidence(frame, disease_dir)
    derived = {}
    for arm in ['a', 'b']:
        prefix = f'gene_{arm}_GSE76092_'
        covered = frame[prefix+'log2fc'].notna()
        sig = frame[prefix+'padj_bh'].lt(0.05)
        high = sig & frame[prefix+'log2fc'].abs().ge(1)
        derived[f'gene_{arm}_gse76092_covered'] = covered
        derived[f'gene_{arm}_gse76092_bh005'] = sig
        derived[f'gene_{arm}_gse76092_high_effect_bh'] = high
        derived[f'gene_{arm}_resistance_expression_state'] = np.select([
            ~covered, sig & frame[prefix+'log2fc'].gt(0), sig & frame[prefix+'log2fc'].lt(0),
        ], ['not_measured', 'resistant_up_BH_lt_005', 'resistant_down_BH_lt_005'], default='measured_without_BH_lt_005_association')
        drug = frame[f'gene_{arm}_GSE254136_screen_log2fc']
        derived[f'gene_{arm}_drug_screen_covered'] = drug.notna()
        derived[f'gene_{arm}_drug_screen_negative'] = drug.lt(0)
    frame = pd.concat([frame, pd.DataFrame(derived, index=frame.index)], axis=1)
    for output, suffix in [
        ('gse76092_covered_gene_count', 'gse76092_covered'),
        ('gse76092_bh005_gene_count', 'gse76092_bh005'),
        ('gse76092_high_effect_bh_gene_count', 'gse76092_high_effect_bh'),
        ('drug_screen_covered_gene_count', 'drug_screen_covered'),
        ('drug_screen_negative_gene_count', 'drug_screen_negative'),
    ]:
        frame[output] = frame[[f'gene_a_{suffix}', f'gene_b_{suffix}']].sum(axis=1)
    frame['target_resistance_association_count'] = 0
    for cell in ['HCT116', 'SW480']:
        associated = np.zeros(len(frame), dtype=bool)
        for arm in ['a', 'b']:
            for dose in ['HD', 'MD']:
                prefix = f'gene_{arm}_GSE299428_{cell}_{dose}_vs_PAR_'
                associated |= frame[prefix+'padj'].lt(0.05) & frame[prefix+'log2fc'].abs().ge(1)
        frame[f'{cell}_resistance_expression_association'] = associated
        frame['target_resistance_association_count'] += associated.astype(int)
    frame['DepMap_strong_single_count'] = sum(frame[f'gene_{arm}_DepMap_HT29_gene_effect'].le(-0.5).astype(int) for arm in ['a','b'])
    frame['raw_joint_unfavorable'] = (~frame.pred_joint.lt(0)).astype(int)
    frame['observed_HT29_single_count'] = frame[['single_a','single_b']].notna().sum(axis=1)
    frame['source_relative_GI_percentile'] = 100 * np.searchsorted(np.sort(observed_ht29_proxy), frame.pred_gi, side='right') / len(observed_ht29_proxy)
    frame['single_effect_sum'] = frame.single_a + frame.single_b
    frame['minimum_single_effect'] = frame[['single_a', 'single_b']].min(axis=1)
    frame['joint_consistency_gap'] = frame.pred_joint - (frame.single_effect_sum + frame.pred_gi)
    frame['pred_joint_minus_stronger_single'] = frame.pred_joint - frame.minimum_single_effect
    frame['pi90_width'] = frame.pi90_high - frame.pi90_low
    frame['interval_sign_class'] = np.select([
        frame.pi90_high.lt(0), frame.pi90_low.gt(0),
    ], ['interval_entirely_negative', 'interval_entirely_positive'], default='interval_includes_zero')
    frame['experimental_operability_status'] = '需合作方确认'
    frame['clinical_handoff_status'] = np.where(frame.observed_HT29_single_count.eq(2), '待验证组合假设；两臂具有HT29单基因测量', '待验证组合假设；至少一臂缺少HT29单基因测量')
    frame['expression_evidence_scope'] = 'acquired-resistance association; direction is not a reversal instruction'
    frame['drug_screen_evidence_scope'] = 'single-gene depletion under oxaliplatin; one PBS control; no inferential FDR'
    frame['external_evidence_scope'] = 'independent source pair measurements; source-defined hit rules; basal growth'

    ext_match = external[external.pair_key.isin(frame.pair_key)].copy()
    ext_groups = {key: group for key, group in ext_match.groupby('pair_key')}
    ext_rows = []
    for key, group in ext_groups.items():
        hits = group.study_hit.astype(bool)
        ht29 = group.cell_line.eq('HT29')
        ext_rows.append({
            'pair_key': key, 'external_measurement_count': len(group),
            'external_study_count': group.study.nunique(), 'external_hit_count': int(hits.sum()),
            'external_measured_nonhit_count': int((~hits).sum()),
            'external_ht29_measurement_count': int(ht29.sum()),
            'external_ht29_hit_count': int((hits & ht29).sum()),
            'external_oxaliplatin_measurement_count': int(group.oxaliplatin_experiment.astype(bool).sum()),
            'external_datasets': ';'.join(sorted(group.dataset.unique())),
            'external_source_dois': ';'.join(sorted(group.source_doi.unique())),
        })
    ext_columns = ['pair_key', 'external_measurement_count', 'external_study_count',
                   'external_hit_count', 'external_measured_nonhit_count', 'external_ht29_measurement_count',
                   'external_ht29_hit_count', 'external_oxaliplatin_measurement_count',
                   'external_datasets', 'external_source_dois']
    frame = frame.merge(pd.DataFrame(ext_rows, columns=ext_columns), on='pair_key', how='left')
    for column in ext_columns[1:8]:
        frame[column] = frame[column].fillna(0).astype(int)
    frame['external_pair_state'] = np.select([
        frame.external_measurement_count.eq(0),
        frame.external_hit_count.gt(0) & frame.external_measured_nonhit_count.gt(0),
        frame.external_hit_count.gt(0),
    ], ['not_covered_unknown', 'measured_context_mixed_hit_and_nonhit', 'measured_source_hit'],
        default='measured_without_source_hit')
    # The frozen external sources contain basal-growth pair screens only and
    # establish no acquired-resistant HT29 plus oxaliplatin pair validation.
    frame['has_independent_applicable_context_pair_validation'] = False
    frame['handoff_independent_background_pair_status'] = np.select([
        frame.external_ht29_measurement_count.gt(0), frame.external_measurement_count.gt(0),
    ], ['independent_HT29_basal_pair_measurement_only', 'independent_other_cell_basal_pair_measurement_only'],
        default='no_independent_same_pair_measurement')
    nonanchor = frame.design_scope.eq('both_nonanchor_extrapolation')
    frame['handoff_requires_exploration'] = nonanchor | ~frame.has_independent_applicable_context_pair_validation
    frame['handoff_evidence_tier'] = np.where(nonanchor,
        'hypothesis_design_extrapolation_and_unvalidated_resistant_context',
        'hypothesis_unvalidated_resistant_context')
    frame['handoff_evidence_reason'] = np.where(nonanchor,
        '双非anchor组合属于源筛选设计外推，当前独立来源未建立获得性耐药HT29及奥沙利铂背景的同pair验证。',
        '含源筛选anchor，当前独立来源仍未建立获得性耐药HT29及奥沙利铂背景的同pair验证。')
    frame = frame.sort_values(['pred_gi', 'pair_key']).reset_index(drop=True)
    frame['rank_gi_all_input'] = np.arange(1, len(frame) + 1)
    frame['rank_gi_within_class'] = frame.groupby('pair_class', sort=False).cumcount() + 1
    core_mask = frame.pair_class.eq('core_metabolism_histone') & ~frame.measured_current_ht29_pair
    core_sd = frame.loc[core_mask, 'between_seed_sd']
    sd_median, sd_q75 = float(core_sd.quantile(0.5)), float(core_sd.quantile(0.75))
    frame['seed_sd_band_order'] = np.select([
        frame.between_seed_sd.le(sd_median), frame.between_seed_sd.le(sd_q75),
    ], [0, 1], default=2)
    frame['seed_sd_band'] = frame.seed_sd_band_order.map({
        0: 'lower_half_of_core_candidate_seed_spread',
        1: 'between_median_and_q75_of_core_candidate_seed_spread',
        2: 'upper_quartile_of_core_candidate_seed_spread',
    })
    frame = single_strata(frame, cutoff)
    core = frame.loc[core_mask].copy()
    top500 = core.head(policy['top500_size']).copy()
    top500['top500_rank'] = np.arange(1, len(top500) + 1)
    # First-round shortlist requires measured HT29 singles; the full ranked pool does not.
    shortlist_pool = top500[top500.observed_HT29_single_count.eq(2)].sort_values(PRIORITY_SORT, ascending=PRIORITY_ASCENDING)
    top = diversified(shortlist_pool, top_n, cap)
    top['top30_rank'] = np.arange(1, len(top) + 1)
    top['in_top30'] = True
    priority = priority_selection(top, cutoff, priority_n, cap)
    priority['in_priority8'] = True
    priority_map = priority.set_index('pair_key').priority8_rank
    top['priority8_rank'] = top.pair_key.map(priority_map)
    top['in_priority8'] = top.pair_key.isin(priority_map.index)
    frame['top30_rank'] = frame.pair_key.map(top.set_index('pair_key').top30_rank)
    frame['priority8_rank'] = frame.pair_key.map(priority_map)
    frame['in_top30'] = frame.pair_key.isin(top.pair_key)
    frame['in_priority8'] = frame.pair_key.isin(priority.pair_key)

    frame.to_csv(out / 'candidates_all.csv', index=False)
    table_names = {
        'core_metabolism_histone': 'core_candidates.csv',
        'epigenetic_extension': 'epigenetic_extension.csv',
        'signaling_extension': 'signaling_extension.csv',
    }
    for category, name in table_names.items():
        frame[frame.pair_class.eq(category) & ~frame.measured_current_ht29_pair].to_csv(out / name, index=False)
    compact = ['top500_rank', 'pair_key', 'gene_a', 'gene_b', 'metabolic_gene', 'chromatin_gene',
               'pair_class', 'pred_gi', 'pred_joint', 'single_a', 'single_b', 'source_relative_GI_percentile',
               'pi90_low', 'pi90_high', 'between_seed_sd', 'model', 'n_models', 'observed_HT29_single_count',
               'design_scope', 'handoff_evidence_tier', 'external_pair_state', 'target_resistance_association_count']
    for arm in ['a','b']:
        compact += [f'gene_{arm}_{col}' for col in ['is_direct_histone_enzyme', 'chromatin_mechanism_class', 'reviewed_source_urls',
            'DepMap_HT29_gene_effect', 'DepMap_HCT116_gene_effect', 'DepMap_SW480_gene_effect',
            'TCGA_tumor_median_log2_RSEM', 'sc_tumor_epithelial_mean_detection',
            'flux_HCT116_mapped_reaction_count', 'flux_HCT116_OxR_minus_Control_mean_log1p_abs_bound',
            'ChIP_H3K27ac_OxR_minus_parental_log1p_max_fold_enrichment',
            'GSE76092_log2fc', 'GSE76092_padj_bh', 'GSE299428_HCT116_HD_vs_PAR_log2fc', 'GSE299428_HCT116_HD_vs_PAR_padj',
            'GSE299428_SW480_HD_vs_PAR_log2fc', 'GSE299428_SW480_HD_vs_PAR_padj']]
    top500[compact].to_csv(out / 'top500.csv', index=False)
    top.to_csv(out / 'top30.csv', index=False)
    priority.to_csv(out / 'priority8.csv', index=False)
    representatives = representative_selection(priority, policy['representative_size'])
    representatives['representative_rank'] = np.arange(1, len(representatives) + 1)
    representatives['nomination_status'] = 'First-round exploratory nomination; no pair-specific resistant-cell validation'
    representatives.to_csv(out / 'representative3.csv', index=False)
    write_handoff(out, top500, representatives)
    ext_match.to_csv(out / 'candidate_external_pair_evidence.csv', index=False)
    ext_match[ext_match.pair_key.isin(top.pair_key)].to_csv(out / 'top30_external_pair_evidence.csv', index=False)
    gene_paths.to_csv(out / 'candidate_gene_pathway_evidence.csv', index=False)
    pair_paths.to_csv(out / 'candidate_pair_pathway_evidence.csv', index=False)
    shared_paths.to_csv(out / 'candidate_shared_pathway_evidence.csv', index=False)
    pair_paths[pair_paths.pair_key.isin(top.pair_key)].to_csv(out / 'top30_pathway_evidence.csv', index=False)

    sensitivity, memberships = [], []
    for threshold in policy['strong_single_sensitivity']:
        selected = priority_selection(top, threshold, priority_n, cap)
        current = single_strata(core, threshold)
        memberships.append(selected[['pair_key', 'gene_a', 'gene_b', 'priority8_rank',
                                     'strong_single_cutoff', 'strong_single_count', 'known_strong_single_count',
                                     'scaled_single_unknown_count', 'priority_band',
                                     'interval_exploration_group', 'handoff_evidence_tier',
                                     'handoff_requires_exploration']])
        sensitivity.append({
            'strong_single_cutoff': threshold, 'strict_core_candidates': len(core),
            'core_without_strong_single': int(current.strong_single_count.eq(0).sum()),
            'core_negative_GI_without_strong_single': int((current.strong_single_count.eq(0) & current.pred_gi.lt(0)).sum()),
            'core_interval_below_zero_without_strong_single': int((current.strong_single_count.eq(0) & current.pi90_high.lt(0)).sum()),
            'core_with_unknown_scaled_single': int(current.scaled_single_unknown_count.gt(0).sum()),
            'top30_with_strong_single': int(single_strata(top, threshold).strong_single_count.gt(0).sum()),
            'top30_with_unknown_scaled_single': int(single_strata(top, threshold).scaled_single_unknown_count.gt(0).sum()),
            'priority8_size': len(selected),
            'priority8_without_strong_single': int(selected.strong_single_count.eq(0).sum()),
            'priority8_with_strong_single': int(selected.strong_single_count.gt(0).sum()),
            'priority8_with_unknown_scaled_single': int(selected.scaled_single_unknown_count.gt(0).sum()),
            'priority8_interval_exploration_count': int(selected.interval_exploration_group.sum()),
            'priority8_handoff_exploration_count': int(selected.handoff_requires_exploration.sum()),
            'priority8_overlap_with_default': len(set(selected.pair_key) & set(priority.pair_key)),
            'priority8_pair_keys': ';'.join(selected.pair_key),
        })
    pd.DataFrame(sensitivity).to_csv(out / 'single_threshold_sensitivity.csv', index=False)
    pd.concat(memberships, ignore_index=True).to_csv(out / 'priority8_sensitivity_membership.csv', index=False)
    operations = []
    top_genes = sorted(set(top.gene_a) | set(top.gene_b))
    triage = single_evidence.reindex(top_genes)
    for gene in top_genes:
        scaled_value = triage.at[gene, 'colo1_single_control_scaled']
        operations.append({
            'gene': gene, 'single_effect_ht29_raw_lfc': singles.at[gene, 'single_effect'],
            'source_perturbation': singles.at[gene, 'perturbation_type'],
            'source_single_definition': singles.at[gene, 'single_definition'],
            'colo1_single_control_scaled': scaled_value,
            'n_guide_pairs_in_single_estimate': triage.at[gene, 'n_guide_pairs_in_single_estimate'],
            'is_strong_single_at_default': pd.NA if pd.isna(scaled_value) else scaled_value <= cutoff,
            'scaled_single_status': 'unknown' if pd.isna(scaled_value) else 'measured',
            'score_ht29_lfc': triage.at[gene, 'score_ht29_lfc'],
            'score_ht29_relative_rank': triage.at[gene, 'score_ht29_relative_rank'],
            'reviewed_accessions': annotation.at[gene, 'reviewed_accessions'],
            'validated_reagents_and_sequences': '需合作方确认',
            'isoform_and_target_engagement': '需合作方确认',
            'protein_expression_in_parental_and_resistant_cells': '需合作方确认',
            'available_selective_inhibitor_or_orthogonal_perturbation': '需合作方确认',
            'partial_perturbation_dynamic_range': '需合作方确认',
            'matched_drug_exposure_design': '需合作方确认',
        })
    pd.DataFrame(operations).to_csv(out / 'experimental_operability.csv', index=False)

    ranking_policy = {
        'primary_universe': 'Strict metabolic enzyme x direct histone-modifying catalytic enzyme; unmeasured in current HT29 pair screen.',
        'top500': 'First 500 core pairs ordered by predicted additive GI, then canonical pair key. Never fill using readers, cofactors or signaling extensions.',
        'top30': 'Among Top500 with both HT29 single effects observed: evidence triage using explicit sort below, with at most 3 occurrences per gene.',
        'priority8': 'First 8 diversified rows of the same Top30 evidence order.',
        'representative3': 'Lexicographically earliest feasible 3-rank tuple within Priority8 with no repeated gene; exploratory handoff only.',
        'sort_columns': PRIORITY_SORT, 'ascending': PRIORITY_ASCENDING,
        'single_cutoff': cutoff,
        'joint_use': 'Positive or missing raw joint is deprioritized in the evidence shortlist. Source raw-LFC zero is not an absolute killing threshold; this is triage, not SL validation.',
        'GI_use': 'Source-relative empirical percentile is descriptive. Neither negative GI nor an entirely negative interval is interpreted as a synthetic-lethality probability.',
        'missing_singles': 'Main GNN uses independent features, not observed HT29 single inputs. Missing raw single/joint values stay missing; candidate interval calibration may fail under gene/design shifts.',
        'resistant_context': 'HCT116/SW480 resistance expression and HT29/HTOXAR3 expression are associations. DepMap parental cells do not stand in for HCT116-R/SW480-R.',
        'extensions': 'Other epigenetic factors and metabolic signaling combinations are separate CSVs, excluded from core Top500 and core hit metrics.',
        'source_proxy_distribution': source_proxy_summary,
        'config_source': 'configs/experiment.json', 'frozen_candidate_policy': policy,
    }
    (out / 'ranking_policy.json').write_text(json.dumps(ranking_policy, ensure_ascii=False, indent=2) + '\n')
    metadata = {
        'prepared_utc': datetime.now(timezone.utc).isoformat(),
        'n_input_predictions': len(frame), 'n_current_ht29_measured_in_input': int(frame.measured_current_ht29_pair.sum()),
        'class_counts': frame.pair_class.value_counts().to_dict(),
        'n_strict_core_unmeasured': len(core), 'n_top30': len(top), 'n_priority8': len(priority), 'n_top500': len(top500), 'n_representative3': len(representatives),
        'top30_gene_max_appearances': max(Counter(top.gene_a.tolist() + top.gene_b.tolist()).values(), default=0),
        'priority8_gene_max_appearances': max(Counter(priority.gene_a.tolist() + priority.gene_b.tolist()).values(), default=0),
        'top30_strong_single_counts': {str(k): int(v) for k,v in top.strong_single_count.value_counts(dropna=False).sort_index().items()},
        'top30_unknown_scaled_single_pairs': int(top.scaled_single_unknown_count.gt(0).sum()),
        'priority8_strong_single_counts': {str(k): int(v) for k,v in priority.strong_single_count.value_counts(dropna=False).sort_index().items()},
        'priority8_unknown_scaled_single_pairs': int(priority.scaled_single_unknown_count.gt(0).sum()),
        'priority8_band_counts': priority.priority_band.value_counts().to_dict(),
        'priority8_interval_exploration_count': int(priority.interval_exploration_group.sum()),
        'priority8_handoff_exploration_count': int(priority.handoff_requires_exploration.sum()),
        'top30_handoff_exploration_count': int(top.handoff_requires_exploration.sum()),
        'priority8_handoff_evidence_tiers': priority.handoff_evidence_tier.value_counts().to_dict(),
        'priority8_positive_raw_joint_count': int(priority.pred_joint.gt(0).sum()),
        'source_proxy_distribution': source_proxy_summary,
        'top30_external_covered_pairs': int(top.external_measurement_count.gt(0).sum()),
        'top30_external_source_hit_pairs': int(top.external_hit_count.gt(0).sum()),
        'top30_gse76092_covered_genes_per_pair': top.gse76092_covered_gene_count.value_counts().sort_index().to_dict(),
        'maximum_abs_joint_consistency_gap': float(frame.joint_consistency_gap.abs().max()),
        'policy_file': 'ranking_policy.json', 'prediction_source': prediction_metadata,
        'disease_sources': ['results/disease/resistant_vs_parental_untreated_gene_effects.csv',
                            'results/disease/GSE76092_gene_effects.csv', 'results/drug_evidence/disease_evidence.csv'],
        'external_source': 'data/evidence/external/external_pairs.csv',
        'annotation_source': 'data/evidence/genes_annotation.csv',
        'pathway_contrast': PRIMARY,
        'n_candidate_gene_pathway_records': len(gene_paths),
        'n_candidate_pair_arm_pathway_records': len(pair_paths),
        'n_candidate_shared_pathway_records': len(shared_paths),
        'independent_cards': len(top), 'cards_directory': 'evidence_cards',
        'llm_api_calls': 0,
    }
    (out / 'ranking_metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n')
    card_dir = out / 'evidence_cards'
    card_dir.mkdir(parents=True, exist_ok=True)
    cards = []
    card_index = []
    gene_path_groups = {gene: block for gene, block in gene_paths.groupby('gene')}
    shared_path_groups = {key: block for key, block in shared_paths.groupby('pair_key')}
    for row in top.to_dict('records'):
        same_pair = ext_groups.get(row['pair_key'], external.iloc[:0])
        same_path = shared_path_groups.get(row['pair_key'], shared_paths.iloc[:0])
        card = evidence_card(row, same_pair, prediction_metadata, cutoff, gene_path_groups, same_path, gene_paths.iloc[:0], source_proxy_summary)
        filename = f"{int(row['top30_rank']):02d}_{row['metabolic_gene']}__{row['chromatin_gene']}.md"
        (card_dir / filename).write_text(card)
        cards.append(card)
        card_index.append({'top30_rank': row['top30_rank'], 'pair_key': row['pair_key'],
                           'priority8_rank': row['priority8_rank'], 'file': 'evidence_cards/' + filename})
    pd.DataFrame(card_index).to_csv(out / 'evidence_card_index.csv', index=False)
    (out / 'candidate_evidence_cards.md').write_text('\n\n---\n\n'.join(cards) + '\n')
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    import sys
    if '--legacy-gi' in sys.argv:
        sys.argv.remove('--legacy-gi')
        main()
    else:
        from deliver_sl import main as current_main
        current_main()
