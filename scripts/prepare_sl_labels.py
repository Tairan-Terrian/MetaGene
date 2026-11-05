"""Reconstruct replicate- and guide-position-matched HT29 additive contrasts.

These are explicitly reconstructed screening labels, not the paper's GEMINI calls.
No double-perturbation outcome is used to estimate control baselines or singles.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    q[order] = np.minimum(1, np.minimum.accumulate((p[order]*len(order)/np.arange(1, len(order)+1))[::-1])[::-1])
    return q


def main():
    source = ROOT/'data/pairs_sources'
    out = ROOT/'data/sl'; out.mkdir(exist_ok=True)
    lfc = pd.read_csv(source/'HT29_gene_lfc.csv', index_col=0)
    reps = list(lfc)
    ann = pd.read_csv(source/'HT29_annotations.csv')
    ann['key'] = ann.Note+'_'+ann.Gene1+'_'+ann.Gene2
    frame = lfc.join(ann.drop_duplicates('key').set_index('key'))
    mapping = pd.read_csv(ROOT/'data/annotations/symbol_mapping.csv').set_index('source_symbol').gene.to_dict()
    for col in ['Gene1', 'Gene2']:
        frame[col] = frame[col].map(mapping).fillna(frame[col])
    controls = {'ADAD1', 'CYLC2', 'KLK12'}
    # Nonessential cutting controls match the two-cut structure of A+control.
    baseline = frame.loc[frame.MyNote.eq('NONESSENTIAL-NONESSENTIAL'), reps].median()
    nt = frame.loc[frame.MyNote.eq('NONTARGET-NONTARGET'), reps].median()
    s1 = frame.loc[frame.Gene2.isin(controls) & frame.MyNote.isin(['LibrarySingletons','GIControlsSingletons'])].groupby('Gene1')[reps].median()
    s2 = frame.loc[frame.Gene1.isin(controls) & frame.MyNote.isin(['AnchorSingletons','GIControlsSingletons'])].groupby('Gene2')[reps].median()
    pairs = frame.loc[frame.MyNote.isin(['AnchorCombinations','LibraryCombinations','GIControlsCombinations']) & frame.Gene1.ne(frame.Gene2)].copy()
    a = s1.reindex(pairs.Gene1).to_numpy(); b = s2.reindex(pairs.Gene2).to_numpy()
    raw = pairs[reps].to_numpy()-a-b
    delta = raw + baseline.to_numpy()
    keys = np.sort(pairs[['Gene1','Gene2']].to_numpy(), axis=1)
    detail = pd.DataFrame({'gene_a':keys[:,0], 'gene_b':keys[:,1], 'source_key':pairs.index})
    for i,rep in enumerate(reps):
        detail[f'dlfc_r{i+1}'] = delta[:,i]
        detail[f'raw_dlfc_r{i+1}'] = raw[:,i]
        detail[f'single_pos1_r{i+1}'] = a[:,i]
        detail[f'single_pos2_r{i+1}'] = b[:,i]
        detail[f'joint_r{i+1}'] = pairs[rep].to_numpy()
    labels = detail.groupby(['gene_a','gene_b'], as_index=False)[[f'{p}_r{i}' for p in ['dlfc','raw_dlfc'] for i in [1,2,3]]].mean()
    ds = labels[[f'dlfc_r{i}' for i in [1,2,3]]]
    labels['pair_key'] = labels.gene_a+'|'+labels.gene_b
    labels['n_replicates'] = ds.notna().sum(axis=1)
    labels['dlfc'] = ds.mean(axis=1)
    labels['raw_dlfc'] = labels[[f'raw_dlfc_r{i}' for i in [1,2,3]]].mean(axis=1)
    labels['replicate_sd'] = ds.std(axis=1, ddof=1)
    labels['replicate_se'] = labels.replicate_sd/np.sqrt(labels.n_replicates)
    labels['negative_replicate_fraction'] = ds.lt(0).sum(axis=1)/labels.n_replicates.replace(0,np.nan)
    # One-sided paired replicate t test; orientations are averaged, not pseudoreplicated.
    stat = labels.dlfc/labels.replicate_se
    labels['replicate_p_less_zero'] = student_t.cdf(stat, labels.n_replicates-1)
    labels.loc[labels.n_replicates.lt(2),'replicate_p_less_zero'] = np.nan
    labels['replicate_q_bh'] = bh(labels.replicate_p_less_zero)
    labels['synergy_label'] = np.where(labels.n_replicates.ge(2), labels.dlfc.lt(-1).astype(float), np.nan)
    labels['significant_synergy_support'] = labels.dlfc.lt(-1) & labels.replicate_q_bh.lt(.1)
    for threshold in [.5, 1, 1.5]:
        labels[f'hit_dlfc_lt_minus_{str(threshold).replace(".","p")}'] = labels.dlfc.lt(-threshold)
    genes = set(pd.read_csv(ROOT/'data/processed/model_genes.csv').gene)
    labels['in_feature_domain'] = labels.gene_a.isin(genes) & labels.gene_b.isin(genes)
    labels['used_for_model'] = labels.synergy_label.notna() & labels.in_feature_domain
    labels.to_csv(out/'ht29_labels.csv',index=False)
    detail.to_csv(out/'orientation_replicate_contrasts.csv',index=False)
    s1.assign(position=1).rename_axis('gene').reset_index().to_csv(out/'single_position1.csv',index=False)
    s2.assign(position=2).rename_axis('gene').reset_index().to_csv(out/'single_position2.csv',index=False)
    meta = {'version':'sl_v2', 'threshold':-1, 'event':'control-centered dLFC < -1 with at least 2 measured biological replicates',
            'formula':'dLFC_r = LFC_AB,r - median(LFC_A+control,r at position1) - median(LFC_control+B,r at position2) + median(LFC_nonessential+nonessential,r)',
            'baseline_nonessential_by_replicate':baseline.to_dict(), 'baseline_nontarget_sensitivity':nt.to_dict(),
            'baseline_pair_count':int(frame.MyNote.eq('NONESSENTIAL-NONESSENTIAL').sum()),
            'position1_single_genes':len(s1),'position2_single_genes':len(s2),
            'distinct_pairs':len(labels),'model_pairs':int(labels.used_for_model.sum()),
            'positive_model_pairs':int(labels.loc[labels.used_for_model,'synergy_label'].sum()),
            'replicate_significant_model_pairs':int(labels.loc[labels.used_for_model,'significant_synergy_support'].sum()),
            'excluded_insufficient_replicates':int(labels.n_replicates.lt(2).sum()),
            'missing_label_policy':'Never zero-fill missing singles or missing pair outcomes; missing or one-replicate labels excluded from supervised fit.',
            'statistical_scope':'Replicate t-test and BH are exploratory supporting annotations on aggregated gene-level LFC. They are not author GEMINI FDR, and do not propagate guide/control estimation error.',
            'biological_scope':'Negative synergy is the primary learning event. Single-gene depletion is audited separately to distinguish SL-compatible candidates from synthetic-sick/strong-single combinations.',
            'sources':['https://www.nature.com/articles/s41467-025-67256-9','https://www.nature.com/articles/s41467-024-47795-3','https://pmc.ncbi.nlm.nih.gov/articles/PMC11370406/']}
    (out/'label_definition.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(meta,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
