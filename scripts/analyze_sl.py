"""Post-fit sensitivity, subgroup and independent-study analyses; no model choice."""
from pathlib import Path
import json
import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score,roc_auc_score
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from metagene import sl
from metagene.reproducibility import ensure_run_manifest, resolve_devices


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--devices',nargs='+',help='Execution devices; the first device is used for inference.')
    parser.add_argument('--device',help='Single inference device.')
    args=parser.parse_args()
    if args.device and args.devices:parser.error('Choose --device or --devices, not both.')
    device=resolve_devices(args.devices or ([args.device] if args.device else None))[0]
    out=ROOT/'results/sl';cfg=json.loads((ROOT/'configs/experiment_sl.json').read_text());selection=json.loads((out/'model_selection.json').read_text());selected=selection['selected_gnn']
    ensure_run_manifest(ROOT,cfg,out)
    table,genes,pairs,features,graphs,anchors=sl.load_inputs(ROOT)
    meta=json.loads((ROOT/'data/sl/label_definition.json').read_text())
    control=np.array(list(meta['baseline_nonessential_by_replicate'].values()));nt=np.array(list(meta['baseline_nontarget_sensitivity'].values()))
    labels=pd.read_csv(ROOT/'data/sl/ht29_labels.csv');labels=labels[labels.used_for_model].copy();values=labels[[f'dlfc_r{i}' for i in [1,2,3]]].to_numpy()
    records=[]
    for name,shift in [('nonessential_cutting_control',np.zeros(3)),('nontarget_control_sensitivity',nt-control),('uncentered_sensitivity',-control)]:
        delta=np.nanmean(values+shift,axis=1)
        for threshold in [-1.5,-1.25,-1,-.75,-.5,0]:
            records.append({'baseline':name,'threshold':threshold,'n':len(delta),'positives':int((delta<threshold).sum()),'fraction':float((delta<threshold).mean())})
    pd.DataFrame(records).to_csv(out/'label_threshold_sensitivity.csv',index=False)
    ann=pd.read_csv(ROOT/'data/evidence/genes_annotation.csv').set_index('gene')
    all_labels=pd.read_csv(ROOT/'data/sl/ht29_labels.csv')
    m=set(ann.index[ann.is_strict_metabolic_enzyme.eq(1)]);h=set(ann.index[ann.is_direct_histone_enzyme.eq(1)])
    core=all_labels[(all_labels.gene_a.isin(m)&all_labels.gene_b.isin(h))|(all_labels.gene_b.isin(m)&all_labels.gene_a.isin(h))]
    core.to_csv(out/'core_measured_label_coverage.csv',index=False)
    subset=[]
    for protocol in cfg['protocols']:
        t=pd.read_csv(out/f'ensemble_{protocol}_{selected}_test.csv')
        ma=t.gene_a.map(ann.is_strict_metabolic_enzyme).eq(1);mb=t.gene_b.map(ann.is_strict_metabolic_enzyme).eq(1)
        ha=t.gene_a.map(ann.is_direct_histone_enzyme).eq(1);hb=t.gene_b.map(ann.is_direct_histone_enzyme).eq(1)
        counts=t.gene_a.isin(anchors).astype(int)+t.gene_b.isin(anchors).astype(int)
        for name,mask in [('core_metabolic_histone',(ma&hb)|(mb&ha)),('two_nonanchors',counts.eq(0)),('at_least_one_anchor',counts.gt(0)),('three_replicates',t.n_replicates.eq(3))]:
            frame=t[mask]
            subset.append({'protocol':protocol,'subgroup':name,**(sl.metrics(frame,frame.model_score,frame.pred_dlfc) if len(frame) else {'n':0,'positives':0,'auprc':np.nan})})
    pd.DataFrame(subset).to_csv(out/'subgroup_metrics.csv',index=False)
    external=pd.read_csv(ROOT/'data/evidence/external/external_pairs_model_domain.csv')
    external=external[external.gene_a.isin(genes)&external.gene_b.isin(genes)].copy()
    external['pair_seen_in_ht29_deposit']=external.pair_key.isin(set(pd.read_csv(ROOT/'data/sl/ht29_labels.csv').pair_key))
    gi={g:i for i,g in enumerate(genes)};ep=external[['gene_a','gene_b']].apply(lambda s:s.map(gi)).to_numpy(np.int64)
    outputs=[]
    for seed in cfg['seeds']:
        fit=sl.load(out/'runs'/f'final_{selected}_{seed}',device);outputs.append(sl.predict(fit,features,ep,graphs)['score'])
    external['ht29_transfer_score']=np.mean(outputs,axis=0)
    external['prediction_scope']='Uncalibrated transfer ranking from HT29; source-defined external labels retained'
    external.to_csv(out/'external_predictions.csv',index=False)
    rows=[]
    for dataset,group in external.groupby('dataset'):
        for scope,frame in [('all_available',group),('pair_not_in_ht29',group[~group.pair_seen_in_ht29_deposit])]:
            if not len(frame):continue
            y=frame.study_hit.eq(True).to_numpy(int);score=frame.ht29_transfer_score
            rows.append({'dataset':dataset,'scope':scope,'n':len(frame),'positives':int(y.sum()),'prevalence':float(y.mean()),'auprc':float(average_precision_score(y,score)) if y.sum() else np.nan,'auroc':float(roc_auc_score(y,score)) if 0<y.sum()<len(y) else np.nan,'source_doi':frame.source_doi.iloc[0],'hit_definition':frame.hit_definition.iloc[0]})
    pd.DataFrame(rows).to_csv(out/'external_metrics.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
