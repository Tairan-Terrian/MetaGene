"""Run registered SL-v2 benchmarks, final ensembles and retrospective acquisition.

Independent GPU workers are processes, with fixed partition manifests and seed lists.
All heavyweight outputs stay under results/sl/runs on the compute server.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import subprocess
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from metagene import sl
from metagene.reproducibility import (candidate_annotations_path, device_name,
                                      ensure_run_manifest, resolve_devices)


def write_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')


def jobs(config):
    return [(p,m,s) for p in config['protocols'] for m in config['models']
            for s in (config['seeds'] if p=='pair_holdout' else config['secondary_seeds'])]


def benchmark(args,config,inputs,out):
    table,genes,pairs,features,graphs,anchors=inputs
    for i,(protocol,model,seed) in enumerate(jobs(config)):
        if i%args.workers!=args.worker:continue
        directory=out/'runs'/f'{protocol}_{model}_{seed}'
        if (directory/'run.json').exists():continue
        part=sl.partition(table,protocol,config,anchors)
        directory.mkdir(parents=True,exist_ok=True);np.savez_compressed(directory/'partitions.npz',**part)
        result=sl.fit(features,pairs,table,graphs,part['train'],part['val'],model,seed,config,args.device)
        sl.save(result,directory)
        # Persist predictions before scoring; architecture selection later reads validation only.
        metrics={}
        for name in ['val','calibration','test']:
            pred=sl.predict(result,features,pairs[part[name]],graphs)
            frame=table.iloc[part[name]].copy();frame['model_score']=pred['score'];frame['pred_dlfc']=pred['dlfc']
            frame['row_id']=part[name];frame.to_csv(directory/f'{name}_predictions.csv',index=False)
            if name=='val':metrics=sl.metrics(frame,pred['score'],pred['dlfc'])
        write_json(directory/'run.json',{'protocol':protocol,'model':model,'seed':seed,'completed':True,'validation':metrics,
                   'best_epoch':result.best_epoch,'elapsed_seconds':result.elapsed,'device':args.device,
                   'gpu':device_name(args.device),
                   'partition_counts':{k:len(v) for k,v in part.items()},'test_used_for_training_or_selection':False})
        print(json.dumps({'completed':directory.name,'val_ap':metrics['auprc'],'seconds':result.elapsed}),flush=True)
        del result;torch.cuda.empty_cache()


def bootstrap(table,score,replicates=500):
    rng=np.random.default_rng(2203);genes=np.unique(table[['gene_a','gene_b']].values);gindex={g:i for i,g in enumerate(genes)}
    a=table.gene_a.map(gindex).to_numpy();b=table.gene_b.map(gindex).to_numpy();y=table.synergy_label.to_numpy()
    values=[]
    for _ in range(replicates):
        w=rng.multinomial(len(genes),np.full(len(genes),1/len(genes)));weights=w[a]*w[b]
        if np.sum(weights*y)>0 and np.sum(weights*(1-y))>0:
            values.append(average_precision_score(y,score,sample_weight=weights))
    return {'metric':'AUPRC','unit':'dyadic gene cluster bootstrap, product of endpoint multiplicities',
            'requested':replicates,'valid':len(values),'lower_95':float(np.quantile(values,.025)) if values else None,
            'upper_95':float(np.quantile(values,.975)) if values else None}


def evaluate(config,inputs,out):
    records=[json.loads((out/'runs'/f'{p}_{m}_{s}'/'run.json').read_text()) for p,m,s in jobs(config)]
    vals=pd.DataFrame([{'protocol':r['protocol'],'model':r['model'],'seed':r['seed'],**r['validation']} for r in records])
    vals.to_csv(out/'validation_metrics.csv',index=False)
    ranks=vals[vals.protocol.eq('pair_holdout')].groupby('model').auprc.mean().sort_values(ascending=False)
    selected=ranks.loc[['string','multigraph']].idxmax();best=ranks.index[0]
    selection={'selected_gnn':selected,'best_unrestricted_model':best,'criterion':'mean validation AUPRC, primary dLFC < -1',
               'validation_scores':ranks.to_dict(),'test_or_calibration_used':False}
    write_json(out/'model_selection.json',selection)
    evaluated=[];ensembles=[];cis=[]
    for protocol in config['protocols']:
        seeds=config['seeds'] if protocol=='pair_holdout' else config['secondary_seeds']
        for model in config['models']:
            cals=[];tests=[]
            for seed in seeds:
                directory=out/'runs'/f'{protocol}_{model}_{seed}'
                c=pd.read_csv(directory/'calibration_predictions.csv');t=pd.read_csv(directory/'test_predictions.csv')
                cal=sl.calibrate(c.model_score,c.synergy_label);prob=sl.probability(cal,t.model_score)
                evaluated.append({'protocol':protocol,'model':model,'seed':seed,**sl.metrics(t,t.model_score,t.pred_dlfc,prob)})
                cals.append(c);tests.append(t)
            c=cals[0].copy();t=tests[0].copy()
            assert all(np.array_equal(x.pair_key,t.pair_key) for x in tests)
            c['model_score']=np.mean([x.model_score for x in cals],axis=0)
            t['model_score']=np.mean([x.model_score for x in tests],axis=0);t['pred_dlfc']=np.mean([x.pred_dlfc for x in tests],axis=0)
            cal=sl.calibrate(c.model_score,c.synergy_label);t['synergy_probability']=sl.probability(cal,t.model_score)
            t.to_csv(out/f'ensemble_{protocol}_{model}_test.csv',index=False)
            write_json(out/f'calibration_{protocol}_{model}.json',cal)
            ensembles.append({'protocol':protocol,'model':model,**sl.metrics(t,t.model_score,t.pred_dlfc,t.synergy_probability)})
            if model==selected:cis.append({'protocol':protocol,**bootstrap(t,t.model_score)})
    pd.DataFrame(evaluated).to_csv(out/'test_metrics.csv',index=False)
    pd.DataFrame(ensembles).to_csv(out/'ensemble_metrics.csv',index=False)
    write_json(out/'bootstrap_intervals.json',cis)
    summary=pd.DataFrame(evaluated).groupby(['protocol','model']).agg({c:['mean','std'] for c in ['auprc','auroc','precision_at_50','enrichment_at_50','recall_at_500','dlfc_rmse','brier']})
    summary.columns=['_'.join(x) for x in summary.columns];summary.reset_index().to_csv(out/'benchmark_summary.csv',index=False)
    print(json.dumps(selection),flush=True)


def final(args,config,inputs,out):
    table,genes,pairs,features,graphs,anchors=inputs
    selection=json.loads((out/'model_selection.json').read_text());part=sl.partition(table,'pair_holdout',config,anchors)
    models=list(dict.fromkeys([selection['selected_gnn'],selection['best_unrestricted_model']]))
    cj=[(m,s) for m in models for s in config['seeds']]
    core=pd.read_csv(candidate_annotations_path(ROOT))
    mapping={g:i for i,g in enumerate(genes)};cpairs=core[['gene_a','gene_b']].apply(lambda c:c.map(mapping)).to_numpy(np.int64)
    for i,(model,seed) in enumerate(cj):
        if i%args.workers!=args.worker:continue
        directory=out/'runs'/f'final_{model}_{seed}'
        if (directory/'run.json').exists():continue
        train=np.concatenate([part['train'],part['test']])
        result=sl.fit(features,pairs,table,graphs,train,part['val'],model,seed,config,args.device);sl.save(result,directory)
        np.savez_compressed(directory/'partitions.npz',train=train,val=part['val'],calibration=part['calibration'])
        for name,ps,frame in [('candidates',cpairs,core[['gene_a','gene_b','pair_key']]),('calibration',pairs[part['calibration']],table.iloc[part['calibration']])]:
            pred=sl.predict(result,features,ps,graphs)
            data=frame.copy();data['model_score']=pred['score'];data['pred_dlfc']=pred['dlfc'];data.to_csv(directory/f'{name}_predictions.csv',index=False)
        write_json(directory/'run.json',{'model':model,'seed':seed,'best_epoch':result.best_epoch,'elapsed_seconds':result.elapsed,'n_train':len(train),'n_calibration':len(part['calibration'])})
        print(json.dumps({'completed':directory.name}),flush=True)
        del result;torch.cuda.empty_cache()


def active(args,config,inputs,out):
    table,genes,pairs,features,graphs,anchors=inputs
    selected=json.loads((out/'model_selection.json').read_text())['selected_gnn']
    part=sl.partition(table,'pair_holdout',config,anchors);pool=part['train'];ac=config['active_learning']
    initial=max(20,round(len(pool)*ac['initial_fraction']));batch=max(1,round(len(pool)*ac['query_fraction']))
    aljobs=[(s,strategy) for s in ac['seeds'] for strategy in ac['strategies']]
    for i,(seed,strategy) in enumerate(aljobs):
        if i%args.workers!=args.worker:continue
        directory=out/'runs'/f'active_{strategy}_{seed}';directory.mkdir(parents=True,exist_ok=True)
        if (directory/'run.json').exists():continue
        rng=np.random.default_rng(seed);seen=list(rng.permutation(pool)[:initial]);curves=[];queries=[]
        cfg={**config,'max_epochs':min(config['max_epochs'],60),'patience':10}
        for rnd in range(ac['rounds']+1):
            result=sl.fit(features,pairs,table,graphs,np.array(seen),part['val'],selected,seed+1009*rnd,cfg,args.device)
            pred=sl.predict(result,features,pairs[part['test']],graphs)
            truth=table.iloc[seen]
            curves.append({'strategy':strategy,'seed':seed,'round':rnd,'n_measured':len(seen),'strong_hits':int(truth.synergy_label.sum()),'tail_hits_dlfc_lt_minus_075':int(truth.dlfc.lt(-.75).sum()),'best_observed_dlfc':float(truth.dlfc.min()),'test_auprc':sl.metrics(table.iloc[part['test']],pred['score'])['auprc']})
            pd.DataFrame(curves).to_csv(directory/'curves.csv',index=False)
            if rnd==ac['rounds']:break
            remaining=pool[~np.isin(pool,seen)]
            if strategy=='random':chosen=rng.choice(remaining,min(batch,len(remaining)),replace=False)
            else:
                pred=sl.predict(result,features,pairs[remaining],graphs,passes=ac['mc_passes'],seed=seed+31*rnd)
                p=1/(1+np.exp(-pred['score']))
                value={'probability':p,'uncertainty':pred['uncertainty'],'upper_confidence':p+pred['uncertainty']}[strategy]
                chosen=remaining[np.lexsort((remaining,-value))[:batch]]
            for row in chosen:queries.append({'strategy':strategy,'seed':seed,'round':rnd+1,'row_id':int(row),'pair_key':table.pair_key.iloc[row],'revealed_dlfc':float(table.dlfc.iloc[row]),'revealed_strong_hit':int(table.synergy_label.iloc[row])})
            seen.extend(chosen.tolist());pd.DataFrame(queries).to_csv(directory/'queries.csv',index=False)
            del result;torch.cuda.empty_cache()
            print(json.dumps({'active':strategy,'seed':seed,'round':rnd,'hits':curves[-1]['strong_hits']}),flush=True)
        write_json(directory/'run.json',{'completed':True,'strategy':strategy,'seed':seed,'rounds':ac['rounds'],'model':selected,'pool_size':len(pool),'pool_hits':int(table.synergy_label.iloc[pool].sum()),'calibration_and_test_labels_queried':False})


def summarize_active(config,out):
    ac=config['active_learning'];frames=[pd.read_csv(out/'runs'/f'active_{strategy}_{seed}'/'curves.csv') for seed in ac['seeds'] for strategy in ac['strategies']]
    curves=pd.concat(frames,ignore_index=True);curves.to_csv(out/'active_learning_curves.csv',index=False)
    agg=curves.groupby(['strategy','round']).agg(n_measured=('n_measured','mean'),strong_hits_mean=('strong_hits','mean'),strong_hits_sd=('strong_hits','std'),tail_hits_mean=('tail_hits_dlfc_lt_minus_075','mean'),test_auprc_mean=('test_auprc','mean')).reset_index()
    agg.to_csv(out/'active_learning_summary.csv',index=False)


def orchestrate(config,out,devices):
    out.mkdir(parents=True,exist_ok=True)
    for phase in ['benchmark','evaluate','final','active','summarize']:
        processes=[];handles=[]
        for worker in range(1 if phase in ['evaluate','summarize'] else len(devices)):
            log=(out/f'{phase}_worker{worker}.log').open('a');handles.append(log)
            command=[sys.executable,__file__,'--phase',phase,'--worker',str(worker),'--workers',str(len(devices)),'--device',devices[worker]]
            processes.append(subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT))
        codes=[p.wait() for p in processes]
        for h in handles:h.close()
        if any(codes):raise RuntimeError(f'{phase} failed: {codes}')
        print(json.dumps({'phase_completed':phase}),flush=True)
    write_json(out/'completed.json',{'completed_unix':time.time(),'benchmark_fits':len(jobs(config)),'active_learning_fits':len(config['active_learning']['strategies'])*len(config['active_learning']['seeds'])*(config['active_learning']['rounds']+1)})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=['all','benchmark','evaluate','final','active','summarize'],default='all')
    parser.add_argument('--worker',type=int,default=0);parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--device',help='Single device for a worker or serial execution.')
    parser.add_argument('--devices',nargs='+',help='Execution devices; defaults to METAGENE_DEVICES or automatic CUDA/CPU selection.')
    args=parser.parse_args();config=json.loads((ROOT/'configs/experiment_sl.json').read_text());out=ROOT/'results/sl';out.mkdir(parents=True,exist_ok=True)
    if args.device and args.devices:parser.error('Choose --device or --devices, not both.')
    if args.workers<1 or not 0<=args.worker<args.workers:parser.error('Require workers >= 1 and 0 <= worker < workers.')
    devices=resolve_devices(args.devices or ([args.device] if args.device else None))
    args.device=devices[args.worker%len(devices)]
    ensure_run_manifest(ROOT,config,out)
    if args.phase=='all':orchestrate(config,out,devices);return
    inputs=sl.load_inputs(ROOT)
    if args.phase=='benchmark':benchmark(args,config,inputs,out)
    elif args.phase=='evaluate':evaluate(config,inputs,out)
    elif args.phase=='final':final(args,config,inputs,out)
    elif args.phase=='active':active(args,config,inputs,out)
    else:summarize_active(config,out)


if __name__=='__main__':main()
