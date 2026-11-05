"""Build the reviewed SL-v2 candidate handoff, with no invented interaction evidence."""
from pathlib import Path
import json
import sys
import shutil
from collections import Counter
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from metagene import sl
from metagene.reproducibility import candidate_annotations_path, ensure_run_manifest

METABOLIC={
 'ADK':'催化腺苷磷酸化，参与嘌呤补救','ADPRM':'水解ADP核糖等焦磷酸底物，参与核苷酸代谢',
 'AK6':'催化腺嘌呤核苷酸磷酸基团转移，参与核苷酸稳态','CHKA':'催化胆碱磷酸化，参与磷脂合成',
 'COQ2':'参与辅酶Q生物合成','DGUOK':'参与线粒体脱氧嘌呤核苷补救','DHODH':'连接线粒体醌代谢与嘧啶从头合成',
 'DLST':'参与α-酮戊二酸脱氢酶复合物催化','DUT':'水解dUTP，限制尿嘧啶错误掺入DNA','FDPS':'参与异戊二烯前体合成',
 'FH':'催化延胡索酸与苹果酸转化','GALK1':'催化半乳糖磷酸化','GART':'催化嘌呤从头合成反应','GNPAT':'参与醚磷脂生物合成',
 'HCCS':'催化c型细胞色素的血红素装配，支持线粒体电子传递','HK2':'催化葡萄糖磷酸化，连接糖酵解与生物合成',
 'IMPDH2':'催化IMP向XMP转化，支持鸟嘌呤核苷酸供给','MAT2A':'合成甲基供体S-腺苷甲硫氨酸（SAM）',
 'MTHFD2':'参与线粒体叶酸一碳代谢','MVK':'参与甲羟戊酸通路','NAMPT':'参与NAD补救合成','NME2':'参与核苷二磷酸与三磷酸间转化',
 'PANK1':'催化泛酸磷酸化，参与辅酶A合成','PDSS2':'参与辅酶Q侧链合成','PFAS':'参与嘌呤从头合成',
 'PGD':'参与磷酸戊糖途径并产生NADPH','PGK1':'参与糖酵解中的磷酸基团转移','PKM':'催化糖酵解末端丙酮酸生成',
 'PMVK':'参与甲羟戊酸通路','PNP':'参与嘌呤核苷分解和补救','PPAT':'参与嘌呤从头合成起始反应',
 'PPCS':'参与辅酶A合成','PSPH':'参与丝氨酸生物合成','RRM1':'参与核糖核苷酸还原和脱氧核苷酸供给',
 'TXNRD1':'维持硫氧还蛋白还原状态','TYMS':'催化dTMP合成','UCK2':'催化尿苷/胞苷磷酸化，参与嘧啶补救','UMPS':'参与UMP从头合成'}
HISTONE={
 'HDAC1':'Zn依赖组蛋白去乙酰化酶','HDAC2':'Zn依赖组蛋白去乙酰化酶','HDAC3':'Zn依赖组蛋白去乙酰化酶',
 'KAT5':'组蛋白乙酰转移酶','KAT6A':'组蛋白乙酰转移酶','KAT6B':'组蛋白乙酰转移酶','KAT7':'组蛋白乙酰转移酶',
 'TAF1':'具有组蛋白乙酰转移酶活性的多功能转录因子','KDM1A':'FAD依赖的组蛋白去甲基化酶',
 'KDM3A':'Fe(II)/α-酮戊二酸依赖的组蛋白去甲基化酶','KDM3B':'Fe(II)/α-酮戊二酸依赖的组蛋白去甲基化酶',
 'PHF8':'Fe(II)/α-酮戊二酸依赖的组蛋白去甲基化酶','PRMT5':'SAM依赖的精氨酸甲基转移酶，具有组蛋白及非组蛋白底物',
 'RNF20':'参与组蛋白H2B单泛素化的E3连接酶','BAP1':'具有组蛋白H2A去泛素化活性的酶'}


def fmt(value):return '缺测' if pd.isna(value) else f'{value:.3g}'


def mechanism(row,annotation):
    m,e=row.metabolic_gene,row.chromatin_gene
    base=f'{m}{METABOLIC[m]}；{e}是{HISTONE[e]}。'
    if m=='MAT2A' and e=='PRMT5':
        bridge='二者存在SAM供给与甲基化催化的生化联系；同一路径的双重抑制也可能呈遮蔽效应，不能由该联系直接推断协同。'
    elif m in {'PANK1','PPCS'} and e.startswith('KAT'):
        bridge='辅酶A代谢可能影响乙酰辅酶A可用性，与组蛋白乙酰化构成间接代谢联系；现有数据未测得该组合下的乙酰辅酶A变化。'
    elif m=='FH' and e in {'KDM3A','KDM3B','PHF8'}:
        bridge='延胡索酸积累可能影响α-酮戊二酸依赖双加氧酶，方向可能产生重叠或拮抗，需以负向组合残差而非通路邻近性判断协同。'
    elif m=='DLST' and e in {'KDM3A','KDM3B','PHF8'}:
        bridge='α-酮戊二酸代谢与去甲基化酶的底物可用性存在联系，但该扰动对核内辅因子浓度和组蛋白标记的净方向尚未测定。'
    elif m in {'ADK','ADPRM','AK6','DGUOK','DHODH','DUT','GART','IMPDH2','NME2','PFAS','PNP','PPAT','RRM1','TYMS','UCK2','UMPS'}:
        bridge='可检验的假说是：核苷酸供给或稳态受扰后，细胞更依赖染色质调控维持复制、转录及损伤应答；同时干预可能产生超出单基因相加预期的适应性下降。'
    elif m=='NAMPT':
        bridge='可检验的假说是NAD代谢压力与染色质调控压力产生组合脆弱性。'
        if e.startswith('HDAC'):bridge+='HDAC1/2/3为Zn依赖酶，不能将这一组合解释为NAD直接供给其催化反应。'
    else:
        bridge='可检验的假说是代谢稳态受扰后，细胞对该组蛋白修饰通路的适应性调节依赖增强；现有功能注释本身不能证明二者形成代偿。'
    s=f'当前单基因审查为{row.single_effect_status}；模型估计dLFC={row.pred_dlfc:.3f}，HT29强协同事件校准分数={row.synergy_probability:.4f}。'
    evidence=('HCT116/SW480至少一臂具有耐药表达关联，支持背景相关性。' if row.target_resistance_association_count>0 else '尚无达到预设阈值的HCT116/SW480耐药表达关联。')
    return base+bridge+s+evidence+'以上为功能事实支持的组合假说，未证明该配对在耐药细胞中合成致死或逆转耐药。'


def diversify(frame,n,cap):
    counts=Counter();chosen=[]
    for idx,row in frame.iterrows():
        if counts[row.gene_a]<cap and counts[row.gene_b]<cap:
            chosen.append(idx);counts.update([row.gene_a,row.gene_b])
        if len(chosen)==n:break
    return frame.loc[chosen].copy()


def candidate_ensemble(root,model,seeds):
    directory=root/'results/sl/runs'
    cs=[pd.read_csv(directory/f'final_{model}_{s}'/'candidates_predictions.csv') for s in seeds]
    calib=[pd.read_csv(directory/f'final_{model}_{s}'/'calibration_predictions.csv') for s in seeds]
    assert all(np.array_equal(x.pair_key,cs[0].pair_key) for x in cs)
    c=cs[0][['gene_a','gene_b','pair_key']].copy();cal=calib[0].copy()
    scores=np.stack([x.model_score for x in cs]);deltas=np.stack([x.pred_dlfc for x in cs])
    c['model_score']=scores.mean(0);c['pred_dlfc']=deltas.mean(0);c['seed_score_sd']=scores.std(0,ddof=1);c['seed_dlfc_sd']=deltas.std(0,ddof=1)
    cal['model_score']=np.mean([x.model_score for x in calib],axis=0);cal['pred_dlfc']=np.mean([x.pred_dlfc for x in calib],axis=0)
    calibration=sl.calibrate(cal.model_score,cal.synergy_label)
    c['synergy_probability']=sl.probability(calibration,c.model_score)
    errors=np.sort(np.abs(cal.dlfc-cal.pred_dlfc));q=errors[min(len(errors)-1,int(np.ceil((len(errors)+1)*.9))-1)]
    c['dlfc_pi90_lower']=c.pred_dlfc-q;c['dlfc_pi90_upper']=c.pred_dlfc+q
    c['calibration_scope']='HT29 screening dLFC < -1; only 3 calibration positives; not resistant-cell SL probability'
    for i,s in enumerate(seeds):c[f'seed{s}_score']=scores[i];c[f'seed{s}_dlfc']=deltas[i]
    calibration['conformal_abs_residual_q90']=float(q);calibration['coverage_assumption']='Exchangeability with HT29 calibration pairs; not guaranteed for two nonanchors or resistant cells'
    return c,calibration


def main():
    out=ROOT/'results/sl';cand=out/'candidates';cand.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((ROOT/'configs/experiment_sl.json').read_text());selection=json.loads((out/'model_selection.json').read_text());model=selection['selected_gnn']
    ensure_run_manifest(ROOT,cfg,out)
    old=pd.read_csv(candidate_annotations_path(ROOT))
    static=['gene_a','gene_b','pair_key','metabolic_gene','chromatin_gene','anchor_gene_count','design_scope','target_resistance_association_count','HCT116_resistance_expression_association','SW480_resistance_expression_association','external_measurement_count','external_hit_count','external_pair_state','external_source_dois']
    static += [c for c in old if c.startswith(('gene_a_','gene_b_')) and not c.endswith('strong_single')]
    frame,cal=candidate_ensemble(ROOT,model,cfg['seeds']);frame=frame.merge(old[static],on=['gene_a','gene_b','pair_key'],validate='one_to_one')
    labels=pd.read_csv(ROOT/'data/sl/ht29_labels.csv');measured=set(labels.pair_key)
    excluded=frame[frame.pair_key.isin(measured)][['gene_a','gene_b','pair_key']].copy();excluded['reason']='Already present in deposited HT29 combination designs, even if incomplete singleton/replicate support'
    excluded.to_csv(cand/'excluded_previously_assayed.csv',index=False)
    frame=frame[~frame.pair_key.isin(measured)].copy()
    if selection['best_unrestricted_model']!=model:
        base,bc=candidate_ensemble(ROOT,selection['best_unrestricted_model'],cfg['seeds'])
        base=base[['pair_key','synergy_probability','pred_dlfc','model_score']].rename(columns={c:f'baseline_{c}' for c in ['synergy_probability','pred_dlfc','model_score']})
        frame=frame.merge(base,on='pair_key');(out/'final_baseline_calibration.json').write_text(json.dumps(bc,indent=2)+'\n')
    scaled=frame[['gene_a_colo1_single_control_scaled','gene_b_colo1_single_control_scaled']]
    frame['single_effect_missing_count']=scaled.isna().sum(axis=1)
    frame['strong_single_count']=scaled.le(cfg['candidate_policy']['strong_single_cutoff']).sum(axis=1)
    frame['single_effect_tier']=np.select([frame.strong_single_count.gt(0),frame.single_effect_missing_count.gt(0)],[2,1],default=0)
    frame['single_effect_status']=frame.single_effect_tier.map({0:'两臂均未触发强单基因耗竭阈值',1:'单基因测量不完整',2:'至少一臂存在强单基因耗竭'})
    frame['model']=model
    frame['strong_synergy_point_prediction']=frame.pred_dlfc.lt(-1)
    frame['strong_synergy_interval_support']=frame.dlfc_pi90_upper.lt(-1)
    frame['high_support_sl_candidate']=frame.single_effect_tier.eq(0)&frame.strong_synergy_interval_support&frame.synergy_probability.ge(.5)
    frame['evidence_level']=np.select([frame.high_support_sl_candidate,frame.single_effect_tier.eq(2),frame.single_effect_tier.eq(1)],['A_强模型支持且单基因初筛兼容','D_强单基因背景协同探索','C_单基因证据缺失探索'],default='B_单基因初筛兼容的探索候选')
    frame['synergy_rank']=frame.model_score.rank(method='first',ascending=False).astype(int)
    frame=frame.sort_values(['single_effect_tier','synergy_probability','pred_dlfc','pair_key'],ascending=[True,False,True,True]).reset_index(drop=True)
    frame['delivery_rank']=np.arange(1,len(frame)+1)
    ann=pd.read_csv(ROOT/'data/evidence/genes_annotation.csv').set_index('gene')
    frame['biological_rationale']=[mechanism(row,ann) for row in frame.itertuples(index=False)]
    frame['metabolic_function_source']=frame.metabolic_gene.map(ann.reviewed_source_urls)
    frame['histone_function_source']=frame.chromatin_gene.map(ann.reviewed_source_urls)
    frame['hypothesis_status']='Function-supported hypothesis; pair-specific synergy and resistance reversal unverified'
    frame['prediction_context']='HT29 basal-growth KO screen; transfer hypothesis for HCT116-R/SW480-R'
    frame['pair_specific_fdr']=np.nan
    frame['pair_specific_fdr_status']='Not available for unmeasured pairs; calibration score is not FDR'
    core=frame.copy();top=frame.head(500).copy();top['top500_rank']=np.arange(1,len(top)+1)
    top30=diversify(top,30,3);top30['top30_rank']=np.arange(1,len(top30)+1)
    priority=diversify(top30,8,2);priority['priority8_rank']=np.arange(1,len(priority)+1)
    representatives=diversify(priority,3,1)
    compact=['delivery_rank','synergy_rank','gene_a','gene_b','pair_key','metabolic_gene','chromatin_gene','synergy_probability','model_score','pred_dlfc','dlfc_pi90_lower','dlfc_pi90_upper','seed_score_sd','seed_dlfc_sd','single_effect_tier','single_effect_status','strong_single_count','single_effect_missing_count','gene_a_colo1_single_control_scaled','gene_b_colo1_single_control_scaled','evidence_level','high_support_sl_candidate','strong_synergy_point_prediction','strong_synergy_interval_support','target_resistance_association_count','HCT116_resistance_expression_association','SW480_resistance_expression_association','external_pair_state','external_measurement_count','external_hit_count','anchor_gene_count','design_scope','model','calibration_scope','prediction_context','pair_specific_fdr','pair_specific_fdr_status','metabolic_function_source','histone_function_source','biological_rationale','hypothesis_status']
    compact += [c for c in frame if c.startswith('baseline_')]
    for arm in ['a','b']:
        compact += [f'gene_{arm}_{c}' for c in ['DepMap_HT29_gene_effect','DepMap_HCT116_gene_effect','DepMap_SW480_gene_effect','TCGA_tumor_median_log2_RSEM','sc_tumor_epithelial_mean_detection'] if f'gene_{arm}_{c}' in frame]
        compact += [c for c in frame if c.startswith(f'gene_{arm}_GSE299428_') and c.endswith(('_log2fc','_padj'))]
    core.to_csv(cand/'core_candidates.csv',index=False)
    for tier,name in [(0,'single_compatible_candidates'),(1,'single_evidence_missing_candidates'),(2,'strong_single_synergy_extension')]:
        core[core.single_effect_tier.eq(tier)][compact].to_csv(cand/f'{name}.csv',index=False)
    for name,f in [('top500',top),('top30',top30),('priority8',priority),('first_round_representatives',representatives)]:
        cols=[c for c in ['top500_rank','top30_rank','priority8_rank'] if c in f]+compact
        f[cols].to_csv(cand/f'{name}.csv',index=False)
    (out/'final_calibration.json').write_text(json.dumps(cal,ensure_ascii=False,indent=2)+'\n')
    fields={
      'synergy_probability':'校准后的HT29 dLFC<-1筛选事件估计；校准集仅3个阳性，绝对值不稳定，不代表耐药细胞SL概率',
      'model_score':'五种子GNN主分类头logit均值，越高越优先',
      'pred_dlfc':'辅助回归头估计的对照校正dLFC；越负代表超出加性预期的抑制越强',
      'dlfc_pi90_lower/dlfc_pi90_upper':'独立校准集绝对残差形成的名义90%预测区间；迁移背景覆盖率未保证',
      'single_effect_tier':'0=两臂有测量且均>-0.5；1=有缺测且未观测到强单效应；2=至少一臂<=-0.5。是初筛标记，不证明非必需',
      'delivery_rank':'单基因兼容层级优先，随后校准协同分数降序、预测dLFC升序、配对键；无绝对联合fitness门槛',
      'synergy_rank':'不考虑单基因层级的分类头分数排序；提供纯协同信号比较',
      'evidence_level':'A需单基因兼容、dLFC区间上界<-1且校准分数>=0.5；B/C/D为不同证据限制的探索队列',
      'pair_specific_fdr':'未测配对无实验P值/FDR，留空；不由模型分数伪造',
      'target_resistance_association_count':'HCT116/SW480的耐药表达关联细胞背景数；非双基因功能证据',
      'biological_rationale':'逐对功能事实、机制假说、预测及局限；来源列链接UniProt审查记录',
      'seed_score_sd/seed_dlfc_sd':'优化随机性，不是实验标准误、FDR或完整预测不确定性',
      'top30_rank':'按交付顺序选择，单基因最多出现3次；多样性选择不改变Top500排序'}
    pd.DataFrame([{'field':k,'definition_zh':v} for k,v in fields.items()]).to_csv(cand/'data_dictionary.csv',index=False)
    review={'version':cfg['version'],'candidate_universe':len(core),'previously_assayed_removed':len(excluded),'top500':len(top),'top30':len(top30),'priority8':len(priority),'representatives':len(representatives),
       'top500_single_effect_tiers':{str(k):int(v) for k,v in top.single_effect_tier.value_counts().items()},'top500_evidence_levels':top.evidence_level.value_counts().to_dict(),
       'top500_high_support_count':int(top.high_support_sl_candidate.sum()),'top500_predicted_dlfc_lt_minus1':int(top.strong_synergy_point_prediction.sum()),'top500_interval_upper_lt_minus1':int(top.strong_synergy_interval_support.sum()),
       'top500_target_resistance_association':int(top.target_resistance_association_count.gt(0).sum()),'top500_external_measured':int(top.external_measurement_count.gt(0).sum()),
       'model_selection':selection,'calibration':cal,'max_candidate_probability':float(top.synergy_probability.max()),'top500_is_capacity_not_hit_count':True,'absolute_joint_gate':False}
    (out/'delivery_review.json').write_text(json.dumps(review,ensure_ascii=False,indent=2)+'\n')
    from write_sl_documents import render_candidate_cards
    (cand/'candidate_evidence_cards.md').write_text(render_candidate_cards(top30.to_dict('records')))
    handoff=out/'handoff';handoff.mkdir(exist_ok=True)
    feedback=['record_id','pair_key','gene_a','gene_b','cell_model','parental_or_resistant','drug_condition_label','replicate_id','lfc_single_a','lfc_single_b','lfc_double','lfc_control','effect_scale','control_normalization','quality_pass','source_file','measured_at']
    pd.DataFrame(columns=feedback).to_csv(handoff/'measured_feedback_template.csv',index=False)
    top[['top500_rank','pair_key','metabolic_gene','chromatin_gene','synergy_probability','pred_dlfc','evidence_level']].to_csv(handoff/'500_pair_targets.csv',index=False)
    slots=[]
    for row in top.itertuples():
        slots.append({'slot_id':f'S{len(slots)+1:04d}','slot_type':'primary_candidate_capacity','pair_key':row.pair_key,'status':'target-capacity allocation only; no guide sequence'})
    for row in top.head(450).itertuples():
        slots.append({'slot_id':f'S{len(slots)+1:04d}','slot_type':'secondary_design_capacity','pair_key':row.pair_key,'status':'target-capacity allocation only; no guide sequence'})
    while len(slots)<1000:slots.append({'slot_id':f'S{len(slots)+1:04d}','slot_type':'control_or_single_capacity_reservation','pair_key':'','status':'partner to allocate'})
    pd.DataFrame(slots).to_csv(handoff/'1000_slot_capacity_draft.csv',index=False)
    # Replace reader entry points only after all new artifacts exist; retain v1 once.
    archive=ROOT/'results/legacy_v1';archive.mkdir(exist_ok=True)
    for name in ['candidates','handoff']:
        current=ROOT/'results'/name
        if current.exists() and not (archive/name).exists():shutil.copytree(current,archive/name)
        # Preserve any additional results while updating this run's reader files.
        current.mkdir(exist_ok=True)
    for name in ['core_candidates.csv','top500.csv','top30.csv','priority8.csv','candidate_evidence_cards.md','data_dictionary.csv','single_compatible_candidates.csv','single_evidence_missing_candidates.csv','strong_single_synergy_extension.csv','excluded_previously_assayed.csv']:
        shutil.copyfile(cand/name,ROOT/'results/candidates'/name)
    shutil.copyfile(cand/'first_round_representatives.csv',ROOT/'results/handoff/first_round_representatives.csv')
    for p in handoff.iterdir():shutil.copyfile(p,ROOT/'results/handoff'/p.name)
    print(json.dumps(review,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
