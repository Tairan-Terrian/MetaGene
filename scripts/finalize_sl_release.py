"""Record current release status and a lightweight reproducibility inventory."""
from pathlib import Path
import json
import hashlib
import platform
import sys
from datetime import datetime, timezone
import pandas as pd
import numpy as np
import torch
from write_sl_documents import collect_fit_counts

ROOT = Path(__file__).resolve().parents[1]


def main():
    out = ROOT/'results/sl'
    review = json.loads((out/'delivery_review.json').read_text())
    cfg = json.loads((ROOT/'configs/experiment_sl.json').read_text())
    selection = json.loads((out/'model_selection.json').read_text())
    counts, count_source = collect_fit_counts(out, cfg, selection)
    total_fits = sum(counts.values())
    runtime = {'version':review['version'], 'completed_utc':datetime.now(timezone.utc).isoformat(),
               'python':sys.version, 'platform':platform.platform(), 'numpy':np.__version__,
               'pandas':pd.__version__, 'torch':torch.__version__, 'cuda':torch.version.cuda,
               'gpus':[torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())],
               'actual_fits':counts, 'fit_count_source':count_source, 'total_fits':total_fits,
               'llm_api_calls':0, 'project_root':str(ROOT)}
    for file in [out/'runtime.json', ROOT/'results/runtime.json']:
        file.write_text(json.dumps(runtime, ensure_ascii=False, indent=2)+'\n')
    annotation = pd.read_csv(ROOT/'data/evidence/genes_annotation.csv')
    n_metabolic = int(annotation.is_strict_metabolic_enzyme.sum())
    n_histone = int(annotation.is_direct_histone_enzyme.sum())
    n_candidates = int(review['candidate_universe'])
    labels = json.loads((ROOT/'data/sl/label_definition.json').read_text())
    features = json.loads((ROOT/'data/evidence/integrated_input_manifest.json').read_text())
    tiers = review['top500_single_effect_tiers']
    cal = review['calibration']
    ac = cfg['active_learning']
    slots = pd.read_csv(ROOT/'results/handoff/1000_slot_capacity_draft.csv')
    top = pd.read_csv(ROOT/'results/candidates/top500.csv')
    strong_predictions = int(top.pred_dlfc.lt(-1).sum())
    rows = [
      ('对象定义','代谢酶×直接组蛋白修饰酶','已复核',f'{n_metabolic}×{n_histone}；排除{n_metabolic*n_histone-n_candidates}种既有设计后{n_candidates}对','data/evidence/genes_annotation.csv','当前建模基因范围'),
      ('标签','重复和位置匹配的对照校正dLFC','已重建',f"{labels['model_pairs']:,}对，{labels['positive_model_pairs']}个dLFC<-1事件",'data/sl/label_definition.json',f"监督集合至少{cfg.get('minimum_replicates',2)}个有效匹配重复"),
      ('输入','独立多组学特征与图','已计算并复用',f"{features['genes']}基因×{features['features']:,}维；STRING/Reactome/BioPlex",'data/evidence/integrated_input_manifest.json','独立功能网络'),
      ('输入','DepMap/TCGA/CRC单细胞','已接入','基因依赖与表达背景','data/evidence/requested_sources/acquisition_manifest.json','基因级背景信息'),
      ('输入','HCT116/OxR通量与组蛋白峰','已接入','公开模型通量范围、三类组蛋白峰描述','data/evidence/flux_epigenome_manifest.json','作者模型范围及描述性图谱'),
      ('模型','直接强协同事件分类与排序','已评价',f"{len(cfg['models'])}种模型、{len(cfg['protocols'])}种留出、{counts['benchmark']}次拟合",'results/sl/benchmark_summary.csv',f"验证所选GNN：{selection['selected_gnn']}"),
      ('模型','概率与区间','已计算',f"独立校准集{cal['n']:,}对，{cal['positives']}阳性",'results/sl/final_calibration.json','HT29强协同事件校准'),
      ('外测','独立组合研究','已评价','包含排除HT29已测配对的迁移子集','results/sl/external_metrics.csv',f"当前Top500外部覆盖{review['top500_external_measured']}对"),
      ('耐药背景','表达与单基因筛选','已接入并标注',f"Top500有{review['top500_target_resistance_association']}对达到目标背景表达关联阈值",'results/candidates/top500.csv','候选背景证据'),
      ('主动学习',f"{ac['rounds']}轮回顾性模拟",'已评价',f"{len(ac['strategies'])}策略×{len(ac['seeds'])}种子；记录{counts['active_learning']}次拟合",'results/sl/active_learning_summary.csv','已有实测池回顾性查询'),
      ('候选','Top500排序与逐对解释','已交付',f"{tiers.get('0',0)}单基因兼容、{tiers.get('1',0)}缺测、{tiers.get('2',0)}强单基因扩展",'results/candidates/top500.csv',f'预测均值<-1的候选{strong_predictions}对'),
      ('候选','Top30、Priority8与代表组合','已交付',f"多样性约束下{review['top30']}/{review['priority8']}/{review['representatives']}对",'results/candidates/top30.csv','计算提名'),
      ('交接','基因配对与容量草案','已交付',f"{len(top)}配对；{len(slots):,}占位",'results/handoff/1000_slot_capacity_draft.csv','基因级靶点与容量安排'),
      ('交接','真实结果回传接口','已交付模板','包含单A/单B/双干预/对照及背景、重复、质控','results/handoff/measured_feedback_template.csv','合作方结果回传')]
    pd.DataFrame(rows, columns=['module','item','status','result','artifact','limitation']).to_csv(ROOT/'results/data_requirements_audit.csv', index=False)
    files = sorted([*ROOT.joinpath('metagene').glob('*.py'), *ROOT.joinpath('scripts').glob('*sl*.py'),
                    ROOT/'configs/experiment_sl.json', ROOT/'README.md', *ROOT.joinpath('docs').glob('*.md')])
    inventory = {'version':review['version'],
                 'scope':'Current release code and reading documents. Input checksums are recorded separately; model fits are counted from this run or its compact result records.',
                 'files':[{'path':str(p.relative_to(ROOT)), 'bytes':p.stat().st_size,
                           'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]}
    (out/'release_manifest.json').write_text(json.dumps(inventory, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'total_fits':total_fits, 'inventory_files':len(files)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
