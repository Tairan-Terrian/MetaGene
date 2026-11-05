"""Frozen-model external evaluation, final fits and calibrated candidate predictions."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from run_benchmarks import (load_inputs, make_partitions, effective_split, train_one,
                             load_fitted, model_predict, registered_jobs)
from metagene.evaluation import evaluate_predictions, compute_calibration_radius


def diagnostic_audit(root, data, table, config):
    output = root / "results"
    keys = data.pairs[:, 0] * len(data.gene_features) + data.pairs[:, 1]
    records = []
    for job in registered_jobs(config):
        directory = output / "benchmarks/runs" / job["name"]
        parts = dict(np.load(directory / "partition_indices.npz"))
        train = parts["effective_train"]
        for role in ["val", "calibration", "test"]:
            held = parts[role]
            records.append({**job, "held_role": role, "n_train": len(train), "n_held": len(held),
                "shared_record_ids": len(np.intersect1d(train, held)),
                "shared_pair_identities": len(np.intersect1d(keys[train], keys[held])),
                "shared_gene_identities": len(np.intersect1d(data.pairs[train], data.pairs[held])),
                "scope": "Cross-cell transfer permits paired identities across different source contexts"
                    if job["protocol"] == "transfer_ht29" else "Pair-disjoint double-perturbation labels"})
    pd.DataFrame(records).to_csv(output / "benchmarks/leakage_audit.csv", index=False)
    twohead = pd.concat([pd.read_csv(output / f"benchmarks/runs/pair_holdout_twohead_gnn_{seed}/test_predictions.csv")
                    for seed in config["primary_seeds"]], ignore_index=True)
    ht = table[table.dataset.eq("ht29_colo1")]
    # The reconstructed column equals predicted GI plus measured singles.
    diagnostic = {"twohead_inconsistency_rms": float(np.sqrt(np.mean((twohead.pred_joint-twohead.pred_joint_from_gi)**2)))}
    diagnostic.update(twohead_gi_rmse=float(np.sqrt(np.mean((twohead.pred_gi-twohead.additive_gi_proxy)**2))),
        ht29_n_pairs=len(ht), ht29_n_anchor_pairs=int((ht.anchor_count>0).sum()),
        ht29_anchor_pair_fraction=float((ht.anchor_count>0).mean()),
        record_overlap_max=max(r["shared_record_ids"] for r in records),
        nontransfer_pair_overlap_max=max(r["shared_pair_identities"] for r in records if r["protocol"]!="transfer_ht29"),
        source="Five primary-seed test predictions from the two-head GNN comparison")
    (output / "diagnostics.json").write_text(json.dumps(diagnostic, indent=2)+"\n")
    print(json.dumps(diagnostic), flush=True)


def external_data(root, data):
    table = pd.read_csv(root / "data/evidence/external/external_pairs_model_domain.csv")
    genes = pd.read_csv(root / "data/processed/model_genes.csv").gene
    index = dict(zip(genes, range(len(genes))))
    datasets = sorted(table.dataset.unique())
    features = []
    for name in datasets:
        cell = table.loc[table.dataset.eq(name), "cell_line"].iloc[0]
        features.append([int(cell=="K562"), 0, int(cell=="HT29"), 0, 1])
    contexts = np.concatenate([data.context_features, np.asarray(features, dtype=np.float32)])
    dataset_index = {name:i for i,name in enumerate(datasets)}
    n = len(table)
    prepared = replace(data,
        pairs=table[["gene_a","gene_b"]].apply(lambda c:c.map(index)).to_numpy(np.int64),
        context_ids=table.dataset.map(dataset_index).to_numpy(np.int64)+len(data.context_features),
        source_ids=table.dataset.map(dataset_index).to_numpy(np.int64)+100,
        joint=np.full(n,np.nan,np.float32), gi_proxy=np.full(n,np.nan,np.float32),
        single_effects=table[["single_a","single_b"]].to_numpy(np.float32), context_features=contexts)
    return prepared, table


def evaluate_external(root, data, graphs, config, worker, workers, device):
    output = root / "results/external"
    output.mkdir(parents=True, exist_ok=True)
    prepared, table = external_data(root, data)
    jobs = [j for j in registered_jobs(config) if j["protocol"]=="pair_holdout"]
    for job in jobs:
        if job["job_index"] % workers != worker:
            continue
        engine, fitted = load_fitted(root / "results/benchmarks/runs" / job["name"], device)
        prediction = model_predict(engine, fitted, prepared, np.arange(len(table)), graphs)
        frame = table.copy()
        frame["pred_gi_score"] = prediction["gi_score"]
        frame["prediction_units"] = "Frozen training-standardized GI; external source scale uncalibrated"
        frame.to_csv(output / f"predictions_{job['model']}_{job['seed']}.csv",index=False)
        rows = []
        for dataset, group in frame.groupby("dataset"):
            for subset, part in [("all_covered",group),
                                  ("pair_unseen",group[~group.pair_seen_in_any_training_study])]:
                ix = part.index.to_numpy()
                metrics = evaluate_predictions(part.author_gi,part.pred_gi_score,prepared.pairs[ix],device=device,
                    seed=job["seed"],bootstrap_replicates=config["uncertainty"]["bootstrap_replicates"])
                for name in ["rmse","mae","rmse_ci_low","rmse_ci_high","mae_ci_low","mae_ci_high"]:
                    metrics[name]=None
                hits = part.study_hit.to_numpy(bool)
                metrics.update(dataset=dataset,subset=subset,model=job["model"],seed=job["seed"],
                    n_hits=int(hits.sum()), baseline_prevalence=float(hits.mean()),
                    auprc=float(average_precision_score(hits,-part.pred_gi_score)) if hits.any() else None,
                    singleton_policy="Observed external singletons normalized by frozen training-global moments",
                    context_policy="Actual myeloid/colorectal indicator when applicable; other tissue indicators zero; knockout=1",
                    pair_labels_used_in_fitting=False)
                rows.append(metrics)
        pd.DataFrame(rows).to_csv(output / f"metrics_{job['model']}_{job['seed']}.csv",index=False)
        print(json.dumps({"external_completed":job["name"],"n_pairs":len(table)}),flush=True)
        del fitted


def fit_final(root,data,table,graphs,config,worker,workers,device):
    output=root/"results"
    selection=json.loads((output/"benchmarks/selection.json").read_text())
    model=selection["selected_model"]
    for i,seed in enumerate(config["primary_seeds"]):
        if i%workers!=worker:
            continue
        parts=make_partitions(data,table,"final_fit",seed,config)
        split=effective_split(parts,model,"final_fit",seed)
        directory=output/f"final_models/{model}_{seed}"
        directory.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(directory/"partition_indices.npz",**parts,effective_train=split.train)
        engine,fitted=train_one(data,split,model,config,device,graphs,directory,
            progress=lambda row:print(json.dumps({"final_model":model,"seed":seed,**row}),flush=True))
        details=fitted[1] if engine=="xgb" else fitted
        metadata={"model":model,"seed":seed,"engine":engine,"device":device,"completed":True,
                  "best_epoch":details["best_epoch"] if engine=="xgb" else details.best_epoch,
                  "elapsed_seconds":details["elapsed_seconds"] if engine=="xgb" else details.elapsed_seconds,
                  "n_train":len(split.train),"n_val":len(split.val),"n_calibration":len(parts["calibration"]),
                  "partition_seed":config["partition_seed"],"selection_uses_primary_validation":True}
        (directory/"run.json").write_text(json.dumps(metadata,indent=2)+"\n")
        cal=model_predict(engine,fitted,data,parts["calibration"],graphs)
        pd.DataFrame({"row_id":parts["calibration"],"observed_gi":data.gi_proxy[parts["calibration"]],
                      "pred_gi":cal["gi"]}).to_csv(directory/"calibration_predictions.csv",index=False)
        print(json.dumps({"final_fit_complete":metadata}),flush=True)


def predict_candidates(root,data,observed,graphs,config,device):
    output=root/"results"
    model=json.loads((output/"benchmarks/selection.json").read_text())["selected_model"]
    annotation=pd.read_csv(root/"data/evidence/genes_annotation.csv").set_index("gene")
    model_genes=pd.read_csv(root/"data/processed/model_genes.csv").gene
    index=dict(zip(model_genes,range(len(model_genes))))
    single=pd.read_csv(root/"data/processed/singles.csv")
    single=single[single.dataset.eq("ht29_colo1")].set_index("gene").single_effect
    supported=annotation.index
    active=annotation.loc[supported]
    metabolic=active.index[(active.is_strict_metabolic_enzyme==1)|(active.is_extended_metabolic_signaling==1)]
    chromatin=active.index[active.is_chromatin_regulator==1]
    pairs=sorted({tuple(sorted((a,b))) for a in metabolic for b in chromatin if a!=b})
    measured={tuple(row) for row in observed.loc[observed.dataset.eq("ht29_colo1"),["gene_a","gene_b"]].to_numpy()}
    frame=pd.DataFrame([pair for pair in pairs if pair not in measured],columns=["gene_a","gene_b"])
    frame["pair_key"]=frame.gene_a+"|"+frame.gene_b
    frame["single_a"]=frame.gene_a.map(single)
    frame["single_b"]=frame.gene_b.map(single)
    n=len(frame)
    candidate=replace(data,pairs=frame[["gene_a","gene_b"]].apply(lambda c:c.map(index)).to_numpy(np.int64),
        context_ids=np.full(n,2,np.int64),source_ids=np.full(n,2,np.int64),
        joint=np.full(n,np.nan,np.float32),gi_proxy=np.full(n,np.nan,np.float32),
        single_effects=frame[["single_a","single_b"]].to_numpy(np.float32))
    gi,calibration=[],[]
    for seed in config["primary_seeds"]:
        directory=output/f"final_models/{model}_{seed}"
        engine,fitted=load_fitted(directory,device)
        prediction=model_predict(engine,fitted,candidate,np.arange(n),graphs)
        gi.append(prediction["gi"])
        frame[f"pred_gi_seed{seed}"]=prediction["gi"]
        calibration.append(pd.read_csv(directory/"calibration_predictions.csv"))
    cal=calibration[0][["row_id","observed_gi"]].copy()
    cal["ensemble_pred_gi"]=np.mean([f.pred_gi.to_numpy() for f in calibration],axis=0)
    radius=compute_calibration_radius(abs(cal.observed_gi-cal.ensemble_pred_gi),config["uncertainty"]["prediction_interval_level"])
    frame["pred_gi"]=np.mean(gi,axis=0)
    frame["pred_joint"]=frame.pred_gi+frame.single_a+frame.single_b
    frame["between_seed_sd"]=np.std(gi,axis=0,ddof=1)
    frame["pi90_low"]=frame.pred_gi-radius
    frame["pi90_high"]=frame.pred_gi+radius
    frame["ht29_single_coverage"] = frame[["single_a", "single_b"]].notna().sum(axis=1)
    frame["n_models"]=len(gi)
    frame["model"]=model
    directory=output/"candidates"
    directory.mkdir(parents=True,exist_ok=True)
    frame.to_csv(directory/"predicted_pairs.csv",index=False)
    cal.to_csv(directory/"ensemble_calibration.csv",index=False)
    (directory/"prediction_metadata.json").write_text(json.dumps({
        "model":model,"seeds":config["primary_seeds"],"n_candidates":n,"n_calibration":len(cal),
        "nominal_coverage":config["uncertainty"]["prediction_interval_level"],"calibration_radius":radius,
        "single_partition_seed":config["partition_seed"],"candidate_is_ht29_unmeasured":True,
        "joint_definition":"measured singleton sum + ensemble predicted GI; missing when either HT29 singleton is unavailable",
        "primary_scope":"strict metabolic enzymes x direct histone-modifying enzymes; broader pairs exported only as extensions",
        "candidate_features":"Independent covariates for all model genes; no missing HT29 single-effect imputation in the GI predictor",
        "interval_definition":"Absolute-residual split calibration of the five-fit ensemble; the same calibration labels excluded from primary architecture selection, all five final fits and early stopping",
        "interval_scope":"Marginal calibration on observed screen design; candidate nonanchor and disease shifts can invalidate coverage",
        "seed_sd_scope":"Training randomness conditional on the same final split; not a calibrated interval",
        "llm_api_calls":0},indent=2)+"\n")
    print(json.dumps({"candidate_predictions":n,"ensemble_radius":radius}),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--phase",choices=["audit","external","external-summary","fit-final","predict"],required=True)
    parser.add_argument("--worker",type=int,default=0)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--device",default="cuda:0")
    args=parser.parse_args()
    torch.set_num_threads(4)
    config=json.loads((args.root/"configs/experiment.json").read_text())
    if args.phase=="external-summary":
        out=args.root/"results/external"
        pd.concat([pd.read_csv(p) for p in sorted(out.glob("metrics_*.csv"))],ignore_index=True).to_csv(out/"metrics.csv",index=False)
        return
    data,table,graphs=load_inputs(args.root)
    if args.phase=="audit":diagnostic_audit(args.root,data,table,config)
    elif args.phase=="external":evaluate_external(args.root,data,graphs,config,args.worker,args.workers,args.device)
    elif args.phase=="fit-final":fit_final(args.root,data,table,graphs,config,args.worker,args.workers,args.device)
    else:predict_candidates(args.root,data,table,graphs,config,args.device)


if __name__=="__main__":
    main()
