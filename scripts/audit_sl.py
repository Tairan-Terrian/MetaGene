"""Fail-fast invariants for actual labels, partition isolation and delivery files."""
from pathlib import Path
import json
import hashlib
import sys
import argparse
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from metagene import sl
from metagene.reproducibility import ensure_run_manifest, resolve_devices


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--devices',nargs='+',help='Execution devices; the first device is used for the model audit.')
    parser.add_argument('--device',help='Single audit device.')
    args=parser.parse_args()
    if args.device and args.devices:parser.error('Choose --device or --devices, not both.')
    device=resolve_devices(args.devices or ([args.device] if args.device else None))[0]
    out=ROOT/'results/sl';checks=[]
    def check(name,condition):
        checks.append({'check':name,'passed':bool(condition)})
        if not condition:raise AssertionError(name)
    table,genes,pairs,features,graphs,anchors=sl.load_inputs(ROOT)
    cfg=json.loads((ROOT/'configs/experiment_sl.json').read_text());meta=json.loads((ROOT/'data/sl/label_definition.json').read_text())
    ensure_run_manifest(ROOT,cfg,out)
    detail=pd.read_csv(ROOT/'data/sl/orientation_replicate_contrasts.csv');baseline=list(meta['baseline_nonessential_by_replicate'].values())
    for i,c in enumerate(baseline,1):
        check(f'replicate_{i}_control_centered_contrast',np.allclose(detail[f'dlfc_r{i}'],detail[f'raw_dlfc_r{i}']+c,equal_nan=True))
    # A common additive offset must cancel in AB-A-B+control.
    check('contrast_invariant_to_lfc_origin',np.isclose((-2+7)-(-.2+7)-(-.3+7)+(.1+7),-2-(-.2)-(-.3)+.1))
    check('no_one_replicate_supervised_labels',table.n_replicates.ge(2).all())
    check('event_exactly_dlfc_below_minus1',np.array_equal(table.synergy_label,table.dlfc.lt(-1).astype(float)))
    check('finite_supervised_labels',np.isfinite(table.dlfc).all())
    check('unique_unordered_pairs',not table.pair_key.duplicated().any() and table.gene_a.lt(table.gene_b).all())
    for protocol in cfg['protocols']:
        part=sl.partition(table,protocol,cfg,anchors)
        used=np.concatenate([part[k] for k in ['train','val','calibration','test']])
        check(f'{protocol}_disjoint_rows',len(used)==len(np.unique(used)))
        if protocol!='pair_holdout':
            train=set(table.iloc[part['train']][['gene_a','gene_b']].values.ravel())
            test=table.iloc[part['test']]
            check(f'{protocol}_every_test_pair_contains_held_gene',((~test.gene_a.isin(train))|(~test.gene_b.isin(train))).all())
        for model in cfg['models']:
            seeds=cfg['seeds'] if protocol=='pair_holdout' else cfg['secondary_seeds']
            for seed in seeds:
                saved=np.load(out/'runs'/f'{protocol}_{model}_{seed}'/'partitions.npz')
                check(f'{protocol}_{model}_{seed}_fixed_partitions',all(np.array_equal(saved[k],part[k]) for k in part))
    selection=json.loads((out/'model_selection.json').read_text());selected=selection['selected_gnn']
    fit=sl.load(out/'runs'/f'final_{selected}_{cfg["seeds"][0]}',device)
    pred=sl.predict(fit,features,pairs[:32],graphs)['score'];swapped=sl.predict(fit,features,pairs[:32,::-1].copy(),graphs)['score']
    check('prediction_gene_order_symmetry',np.allclose(pred,swapped,atol=1e-6))
    top=pd.read_csv(out/'candidates/top500.csv');core=pd.read_csv(out/'candidates/core_candidates.csv');ann=pd.read_csv(ROOT/'data/evidence/genes_annotation.csv').set_index('gene')
    check('exactly500_unique_predictions',len(top)==500 and not top.pair_key.duplicated().any())
    check('only_strict_metabolic_and_direct_histone_enzymes',top.metabolic_gene.map(ann.is_strict_metabolic_enzyme).eq(1).all() and top.chromatin_gene.map(ann.is_direct_histone_enzyme).eq(1).all())
    measured=set(pd.read_csv(ROOT/'data/sl/ht29_labels.csv').pair_key)
    check('no_previously_assayed_pairs_in_candidates',not set(core.pair_key)&measured)
    check('finite_predictions_probabilities',np.isfinite(top[['pred_dlfc','model_score','synergy_probability']]).all().all() and top.synergy_probability.between(0,1).all())
    expected=core.sort_values(['single_effect_tier','synergy_probability','pred_dlfc','pair_key'],ascending=[True,False,True,True]).head(500).pair_key
    check('top500_follows_declared_order',np.array_equal(expected,top.pair_key))
    check('no_fabricated_pair_fdr',top.pair_specific_fdr.isna().all())
    check('every_pair_has_function_sources_and_rationale',top[['metabolic_function_source','histone_function_source','biological_rationale']].notna().all().all())
    shortlist=pd.read_csv(out/'candidates/top30.csv')
    check('top30_diversity_cap',pd.concat([shortlist.gene_a,shortlist.gene_b]).value_counts().max()<=3)
    slots=pd.read_csv(out/'handoff/1000_slot_capacity_draft.csv')
    check('capacity_plan_covers_all500',len(slots)==1000 and set(top.pair_key)<=set(slots.pair_key.dropna()))
    check('empty_real_measurement_template',len(pd.read_csv(out/'handoff/measured_feedback_template.csv'))==0)
    part=sl.partition(table,'pair_holdout',cfg,anchors)
    for seed in cfg['active_learning']['seeds']:
        for strategy in cfg['active_learning']['strategies']:
            q=pd.read_csv(out/'runs'/f'active_{strategy}_{seed}'/'queries.csv')
            check(f'active_{strategy}_{seed}_queries_only_training_pool',not q.row_id.duplicated().any() and set(q.row_id)<=set(part['train']))
    manifest=json.loads((out/'input_manifest.json').read_text())
    check('frozen_inputs_unchanged',all(hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'] for r in manifest['files']))
    result={'version':cfg['version'],'passed':all(x['passed'] for x in checks),'checks':checks,'count':len(checks)}
    for path in [out/'integrity_audit.json',ROOT/'results/integrity_audit.json']:path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'passed':result['passed'],'checks':len(checks)}))

if __name__=='__main__':main()
