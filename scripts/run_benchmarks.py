"""Registered current fits, validation-only model choice and held-out evaluation."""
import argparse
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from base_data_io import load_data
from metagene.models import ModelConfig
from metagene.training import Preprocessing, TrainConfig, fit_model, load_fit, predict, save_fit
from metagene.splits import Split
from metagene.residual_models import ResidualModelConfig
from metagene.residual_training import ResidualTrainConfig, fit_residual, load_residual, predict_residual, save_residual


def load_inputs(root):
    data, pairs = load_data(root)
    genes = pd.read_csv(root / "data/evidence/genes_annotation.csv", keep_default_na=False)
    lookup = genes.set_index("gene")
    ma = pairs.gene_a.map(lookup.is_strict_metabolic_enzyme).eq(1)
    mb = pairs.gene_b.map(lookup.is_strict_metabolic_enzyme).eq(1)
    ea = pairs.gene_a.map(lookup.is_chromatin_regulator).eq(1)
    eb = pairs.gene_b.map(lookup.is_chromatin_regulator).eq(1)
    pairs["is_strict_enzyme_chromatin"] = (ma & eb) | (mb & ea)
    ca = pairs.gene_a.map(lookup.is_direct_histone_enzyme).eq(1)
    cb = pairs.gene_b.map(lookup.is_direct_histone_enzyme).eq(1)
    pairs["is_core_metabolism_histone"] = (ma & cb) | (mb & ca)
    mapping = pd.read_csv(root / "data/annotations/symbol_mapping.csv", keep_default_na=False).set_index("source_symbol").gene.to_dict()
    annotation = pd.read_csv(root / "data/pairs_sources/HT29_annotations.csv", keep_default_na=False)
    raw_anchors = annotation.loc[annotation.MyNote.eq("AnchorSingletons"), "Gene2"]
    anchors = set(raw_anchors.map(mapping).fillna(raw_anchors))
    pairs["anchor_a"] = pairs.gene_a.isin(anchors)
    pairs["anchor_b"] = pairs.gene_b.isin(anchors)
    pairs["anchor_count"] = pairs.anchor_a.astype(int) + pairs.anchor_b.astype(int)
    graph_edges = dict(np.load(root / "data/evidence/graphs.npz"))
    return data, pairs, graph_edges


def allocations(values, seed, fractions):
    values = np.random.default_rng(seed).permutation(values)
    nt = round(len(values) * fractions["test"])
    nc = round(len(values) * fractions["calibration"])
    nv = round(len(values) * fractions["validation"])
    return {"test": values[:nt], "calibration": values[nt:nt+nc],
            "val": values[nt+nc:nt+nc+nv], "train": values[nt+nc+nv:]}


def make_partitions(data, table, protocol, seed, config):
    ht = np.flatnonzero(data.context_ids == 2)
    blood = np.flatnonzero(data.context_ids != 2)
    pair_keys = data.pairs[:, 0] * len(data.gene_features) + data.pairs[:, 1]
    partition_seed = config["partition_seed"]
    if protocol == "pair_holdout":
        parts = allocations(ht, partition_seed, config["fractions"])
        held = np.concatenate([parts[k] for k in ["val", "calibration", "test"]])
        blood_train = blood[~np.isin(pair_keys[blood], pair_keys[held])]
    elif protocol == "gene_holdout":
        groups = allocations(np.unique(data.pairs[ht]), partition_seed, config["fractions"])
        parts = {name: ht[np.isin(data.pairs[ht], group).all(axis=1)] for name, group in groups.items()}
        held_genes = np.concatenate([groups[k] for k in ["val", "calibration", "test"]])
        blood_train = blood[~np.isin(data.pairs[blood], held_genes).any(axis=1)]
    elif protocol == "anchor_holdout":
        anchor_genes = np.unique(data.pairs[ht][table.iloc[ht][["anchor_a", "anchor_b"]].to_numpy()])
        groups = allocations(anchor_genes, partition_seed, config["fractions"])
        assignment = np.full(len(data.gene_features), -1, dtype=int)
        for number, name in enumerate(["train", "val", "calibration", "test"]):
            assignment[groups[name]] = number
        labels = assignment[data.pairs[ht]]
        parts = {}
        for number, name in enumerate(["train", "val", "calibration", "test"]):
            mask = ((labels == number) | (labels == -1)).all(axis=1) & (labels == number).any(axis=1)
            if name == "train":
                mask |= (labels == -1).all(axis=1)
            parts[name] = ht[mask]
        held_genes = np.concatenate([groups[k] for k in ["val", "calibration", "test"]])
        blood_train = blood[~np.isin(data.pairs[blood], held_genes).any(axis=1)]
    elif protocol == "transfer_ht29":
        unique = np.unique(pair_keys[blood])
        fractions = {"test": 0, "calibration": 0.10, "validation": 0.15}
        groups = allocations(unique, partition_seed, fractions)
        parts = {name: blood[np.isin(pair_keys[blood], keys)] for name, keys in groups.items()}
        parts["test"] = ht
        blood_train = parts["train"]
    elif protocol == "final_fit":
        parts = allocations(ht, partition_seed, config["fractions"])
        parts["train"] = np.concatenate([parts["train"], parts["test"]])
        parts["test"] = np.empty(0, dtype=np.int64)
        held = np.concatenate([parts["val"], parts["calibration"]])
        blood_train = blood[~np.isin(pair_keys[blood], pair_keys[held])]
    parts["pooled_train"] = (blood_train if protocol == "transfer_ht29"
                              else np.concatenate([parts["train"], blood_train]))
    used = np.concatenate([parts[k] for k in ["train", "val", "calibration", "test"]])
    parts["discarded_ht29"] = ht[~np.isin(ht, used)]
    return parts


def effective_split(parts, model, protocol, seed):
    train = parts["pooled_train"] if model.startswith("pooled_") or protocol == "transfer_ht29" else parts["train"]
    return Split(train, parts["val"], parts["test"], protocol, seed, parts["discarded_ht29"])


def registered_jobs(config):
    jobs = []
    for protocol, models in config["protocols"].items():
        seeds = config["primary_seeds"] if protocol in {"pair_holdout", "gene_holdout"} else config["secondary_seeds"]
        for seed in seeds:
            for model in models:
                jobs.append({"protocol": protocol, "seed": seed, "model": model,
                             "name": f"{protocol}_{model}_{seed}", "job_index": len(jobs)})
    return jobs


def sparse_pair_features(data, preprocessing):
    feature = sparse.csr_matrix(data.gene_features)
    a, b = feature[data.pairs[:, 0]], feature[data.pairs[:, 1]]
    difference = a - b
    difference.data = np.abs(difference.data)
    singles = preprocessing.singles.transform(data.singles().reshape(-1, 1), np.repeat(data.source_ids, 2)).reshape(-1, 2)
    single_features = np.column_stack([singles.sum(axis=1), np.abs(singles[:, 0]-singles[:, 1]),
                                       singles[:, 0]*singles[:, 1]])
    context = data.context_features[data.context_ids]
    return sparse.hstack([a+b, difference, a.multiply(b), sparse.csr_matrix(single_features),
                          sparse.csr_matrix(context)], format="csr", dtype=np.float32)


def fit_xgb(data, split, config, device, output):
    import xgboost as xgb
    preprocessing = Preprocessing.fit(data, split.train)
    features = sparse_pair_features(data, preprocessing)
    targets = preprocessing.labels.transform(data.targets(), data.source_ids)[:, 1]
    params = dict(config["xgboost"])
    rounds = params.pop("num_boost_round")
    stopping = params.pop("early_stopping_rounds")
    params.update(device=device, objective="reg:squarederror", eval_metric="rmse", seed=split.seed, nthread=4)
    train = xgb.QuantileDMatrix(features[split.train], label=targets[split.train], max_bin=params["max_bin"])
    val = xgb.QuantileDMatrix(features[split.val], label=targets[split.val], ref=train, max_bin=params["max_bin"])
    history = {}
    start = time.monotonic()
    booster = xgb.train(params, train, rounds, evals=[(val, "validation")],
                        early_stopping_rounds=stopping, evals_result=history, verbose_eval=50)
    booster.save_model(output / "model.ubj")
    result = {"preprocessing": preprocessing, "split": split, "best_epoch": booster.best_iteration+1,
              "elapsed_seconds": time.monotonic()-start, "history": history, "device": device}
    torch.save(result, output / "preprocessing.pt")
    pd.DataFrame({"epoch": np.arange(1, len(history["validation"]["rmse"])+1),
                  "val_gi_loss": history["validation"]["rmse"]}).to_csv(output / "history.csv", index=False)
    return (booster, result)


def train_one(data, split, model, config, device, graph_edges, output, progress=None):
    common = {key: config[key] for key in ["max_epochs", "patience", "batch_size", "learning_rate", "weight_decay"]}
    if model == "xgb":
        return "xgb", fit_xgb(data, split, config, device, output)
    if model in {"linear", "naiad", "twohead_gnn", "twohead_graphless"}:
        variant = {"linear": "additive", "naiad": "naiad_reimplementation", "twohead_gnn": "gnn", "twohead_graphless": "graphless"}[model]
        result = fit_model(data, split, ModelConfig(variant=variant, hidden_dim=config["hidden_dim"], dropout=config["dropout"]),
                           TrainConfig(seed=split.seed, device=device, **common), progress=progress)
        save_fit(result, output / "model.pt")
        pd.DataFrame(result.history).to_csv(output / "history.csv", index=False)
        return "twohead", result
    graph_names = {"target_graphless": (), "target_string": ("string",), "pooled_graphless": (),
                   "pooled_reactome": ("reactome",), "pooled_string": ("string",),
                   "pooled_multigraph": ("reactome", "string", "bioplex")}[model]
    result = fit_residual(data, split, ResidualModelConfig(hidden_dim=config["hidden_dim"], dropout=config["dropout"], graph_names=graph_names, use_single_effects=False),
                    ResidualTrainConfig(seed=split.seed, device=device, validation_source=None if split.name=="transfer_ht29" else 2,
                        source_weights={int(k): v for k,v in config["source_weights"].items()}, **common),
                    graph_edges=graph_edges, progress=progress)
    save_residual(result, output / "model.pt")
    pd.DataFrame(result.history).to_csv(output / "history.csv", index=False)
    return "residual", result


def load_fitted(directory, device):
    meta = json.loads((directory / "run.json").read_text())
    engine = meta["engine"]
    if engine == "xgb":
        import xgboost as xgb
        booster = xgb.Booster()
        booster.load_model(directory / "model.ubj")
        booster.set_param({"device": device})
        result = (booster, torch.load(directory / "preprocessing.pt", weights_only=False, map_location="cpu"))
        result[1]["device"] = device
    elif engine == "twohead":
        result = load_fit(directory / "model.pt", device)
    else:
        result = load_residual(directory / "model.pt", device)
    return engine, result


def model_predict(engine, result, data, indices, graph_edges, mc_passes=1, seed=None, return_samples=False):
    if engine == "xgb":
        import xgboost as xgb
        booster, metadata = result
        pre = metadata["preprocessing"]
        features = sparse_pair_features(data, pre)
        dm = xgb.DMatrix(features[indices])
        score = booster.predict(dm, iteration_range=(0, metadata["best_epoch"]))
        rows = pre.labels.rows(data.source_ids[indices])
        seen = rows > 0
        gi = np.where(seen, score*pre.labels.std[rows,1]+pre.labels.mean[rows,1], np.nan)
        return {"gi": gi, "gi_score": score, "joint": gi+data.singles()[indices].sum(axis=1),
                "source_scale_seen": seen, "gi_std": np.zeros(len(indices)), "gi_score_std": np.zeros(len(indices))}
    if engine == "twohead":
        pred = predict(result, data, indices, mc_passes=mc_passes, seed=seed)
        pred["gi_score"] = pred["gi_standardized"]
        pred["gi_score_std"] = pred["gi_std_standardized"]
        pred["joint_from_gi"] = pred["gi"] + data.singles()[indices].sum(axis=1)
        return pred
    return predict_residual(result, data, indices, graph_edges=graph_edges, mc_passes=mc_passes,
                      seed=seed, return_samples=return_samples)


def prediction_table(data, table, indices, predicted):
    result = table.iloc[indices].copy()
    result["pred_gi"] = predicted["gi"]
    result["gi_score"] = predicted["gi_score"]
    result["pred_joint"] = np.where(predicted["source_scale_seen"], predicted["joint"], np.nan)
    result["pred_joint_from_gi"] = result.pred_gi + result.single_a + result.single_b
    result["source_scale_seen"] = predicted["source_scale_seen"]
    return result


def fit_job(data, table, graphs, job, config, device, output_root):
    output = output_root / "benchmarks/runs" / job["name"]
    output.mkdir(parents=True, exist_ok=True)
    parts = make_partitions(data, table, job["protocol"], job["seed"], config)
    split = effective_split(parts, job["model"], job["protocol"], job["seed"])
    np.savez_compressed(output / "partition_indices.npz", **parts, effective_train=split.train)
    engine, result = train_one(data, split, job["model"], config, device, graphs, output,
        progress=lambda row: print(json.dumps({"job": job["name"], **row}), flush=True))
    val_pred = model_predict(engine, result, data, parts["val"], graphs)
    cal_pred = model_predict(engine, result, data, parts["calibration"], graphs)
    val_table = prediction_table(data, table, parts["val"], val_pred)
    prediction_table(data, table, parts["calibration"], cal_pred).to_csv(output / "calibration_predictions.csv", index=False)
    val_table.to_csv(output / "validation_predictions.csv", index=False)
    if job["protocol"] == "transfer_ht29":
        pre = result[1]["preprocessing"] if engine == "xgb" else result.preprocessing
        target = pre.labels.transform(data.targets()[parts["val"]], data.source_ids[parts["val"]])[:,1]
        val_rmse = float(np.sqrt(np.mean((target-val_pred["gi_score"])**2)))
    else:
        val_rmse = float(np.sqrt(np.mean((val_table.additive_gi_proxy-val_table.pred_gi)**2)))
    training = result[1] if engine == "xgb" else result
    epoch = training["best_epoch"] if engine == "xgb" else training.best_epoch
    seconds = training["elapsed_seconds"] if engine == "xgb" else training.elapsed_seconds
    metadata = {**job, "engine": engine, "completed": True, "device": device,
                "gpu": torch.cuda.get_device_name(torch.device(device)), "torch": torch.__version__,
                "n_train": len(split.train), "n_val": len(split.val), "n_calibration": len(parts["calibration"]),
                "n_test": len(split.test), "n_discarded_ht29": len(parts["discarded_ht29"]),
                "validation_rmse": val_rmse, "best_epoch": epoch, "elapsed_seconds": seconds,
                "test_evaluated_during_fit": False}
    if engine == "residual":
        metadata.update(learned_graph_gates=result.learned_graph_gates, source_statistics=result.source_statistics)
    (output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"completed_job": job["name"], "validation_rmse": val_rmse}), flush=True)


def select_models(output):
    records = [json.loads(path.read_text()) for path in sorted((output / "benchmarks/runs").glob("*/run.json"))]
    fits = pd.DataFrame(records)
    primary = fits[fits.protocol.eq("pair_holdout")]
    summary = primary.groupby("model").validation_rmse.agg(["mean", "std", "size"]).sort_values("mean")
    candidates = [{"model": name, "mean_validation_rmse": float(row["mean"]),
                   "sd_validation_rmse": float(row["std"]), "n_seeds": int(row["size"])} for name,row in summary.iterrows()]
    al = summary.loc[["target_graphless", "target_string"]].sort_values("mean").index[0]
    gnn_models = ['target_string', 'pooled_string', 'pooled_reactome', 'pooled_multigraph']
    task_model = summary.loc[gnn_models].sort_values('mean').index[0]
    selection = {"selected_model": task_model, "best_overall_model": summary.index[0], "selected_al_model": al,
                 "deployment_scope": "Validation-best GNN without source-specific observed single effects; matches the requested graph model and predicts all core gene pairs. Best unrestricted baseline reported separately.",
                 "validation_metric": "Mean HT29 validation GI RMSE; smaller is better",
                 "candidate_scores": candidates, "selection_time_utc": datetime.now(timezone.utc).isoformat(),
                 "test_metrics_used": False}
    (output / "benchmarks/selection.json").write_text(json.dumps(selection, indent=2)+"\n")
    fits.to_csv(output / "benchmarks/fit_registry.csv", index=False)
    print(json.dumps(selection, indent=2))


def evaluate_job(data, table, graphs, job, config, device, output):
    from metagene.evaluation import evaluate_predictions, compute_calibration_radius, interval_metrics
    directory = output / "benchmarks/runs" / job["name"]
    engine, result = load_fitted(directory, device)
    parts = dict(np.load(directory / "partition_indices.npz"))
    indices = parts["test"]
    pred = model_predict(engine, result, data, indices, graphs)
    test = prediction_table(data, table, indices, pred)
    test.to_csv(directory / "test_predictions.csv", index=False)
    rows = []
    key = data.pairs[:, 0] * len(data.gene_features) + data.pairs[:, 1]
    pair_unseen = ~np.isin(key[indices], key[parts["effective_train"]])
    subsets = [("all", np.ones(len(test), dtype=bool)),
               ("strict_enzyme_chromatin", test.is_strict_enzyme_chromatin.to_numpy()),
               ("core_metabolism_histone", test.is_core_metabolism_histone.to_numpy())]
    if job["protocol"] == "transfer_ht29":
        subsets.append(("pair_unseen", pair_unseen))
    pd.DataFrame([{"protocol":job["protocol"],"model":job["model"],"seed":job["seed"],
                   "subset":name,"n_test":int(chosen.sum()),"metrics_evaluated":int(chosen.sum())>=5}
                  for name,chosen in subsets]).to_csv(directory / "subset_counts.csv", index=False)
    for subset, chosen in subsets:
        ix = np.flatnonzero(chosen)
        if len(ix) < 5:
            continue
        targets = [("additive_gi_proxy", test.additive_gi_proxy.to_numpy(),
                    test.pred_gi.to_numpy() if test.source_scale_seen.all() else test.gi_score.to_numpy())]
        if test.source_scale_seen.all():
            targets.append(("joint", test.joint_effect.to_numpy(), test.pred_joint.to_numpy()))
            if engine == "twohead":
                targets.append(("joint_from_gi", test.joint_effect.to_numpy(), test.pred_joint_from_gi.to_numpy()))
        for target, truth, prediction in targets:
            metrics = evaluate_predictions(truth[ix], prediction[ix], data.pairs[indices[ix]], device=device,
                                 seed=job["seed"], bootstrap_replicates=config["uncertainty"]["bootstrap_replicates"])
            if not test.source_scale_seen.all():
                for name in ["rmse", "mae", "rmse_ci_low", "rmse_ci_high", "mae_ci_low", "mae_ci_high"]:
                    metrics[name] = None
            metrics.update(protocol=job["protocol"], model=job["model"], seed=job["seed"],
                           dataset="ht29_colo1", target=target, subset=subset,
                           source_scale_seen=bool(test.source_scale_seen.all()))
            rows.append(metrics)
    pd.DataFrame(rows).to_csv(directory / "metrics.csv", index=False)
    cal = pd.read_csv(directory / "calibration_predictions.csv")
    if test.source_scale_seen.all():
        radius = compute_calibration_radius(np.abs(cal.pred_gi-cal.additive_gi_proxy),
            coverage=config["uncertainty"]["prediction_interval_level"])
        coverage = interval_metrics(test.additive_gi_proxy, test.pred_gi, radius)
        coverage.update(protocol=job["protocol"], model=job["model"], seed=job["seed"],
                        nominal_coverage=config["uncertainty"]["prediction_interval_level"],
                        radius=radius, n_calibration=len(cal))
        pd.DataFrame([coverage]).to_csv(directory / "interval_metrics.csv", index=False)
    print(json.dumps({"evaluated_job": job["name"]}), flush=True)


def summarize(output):
    runs = pd.concat([pd.read_csv(path) for path in sorted((output / "benchmarks/runs").glob("*/metrics.csv"))], ignore_index=True)
    runs.to_csv(output / "benchmarks/metrics.csv", index=False)
    keys = ["protocol", "dataset", "subset", "target", "model"]
    columns = ["pearson", "spearman", "rmse", "mae", "top_k_overlap", "top_k_enrichment", "n"]
    summary = runs.groupby(keys)[columns].agg(["mean", "std", "count"])
    summary.columns = [f"{metric}_{stat}" for metric,stat in summary.columns]
    summary.reset_index().to_csv(output / "benchmarks/summary.csv", index=False)
    pd.concat([pd.read_csv(path) for path in sorted((output / "benchmarks/runs").glob("*/subset_counts.csv"))], ignore_index=True).to_csv(output / "benchmarks/subset_counts.csv", index=False)
    pd.concat([pd.read_csv(path) for path in sorted((output / "benchmarks/runs").glob("*/interval_metrics.csv"))], ignore_index=True).to_csv(output / "benchmarks/calibration.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--phase", choices=["fit", "select", "evaluate", "summarize"], required=True)
    parser.add_argument("--worker", type=int, default=0)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--start-job", type=int, default=0)
    parser.add_argument("--stop-job", type=int, default=100000)
    parser.add_argument("--only-model")
    parser.add_argument("--only-protocol")
    parser.add_argument("--skip-model")
    args = parser.parse_args()
    config = json.loads((args.root / "configs/experiment.json").read_text())
    output = args.root / "results"
    output.mkdir(parents=True, exist_ok=True)
    if args.phase == "select":
        select_models(output)
        return
    if args.phase == "summarize":
        summarize(output)
        return
    torch.set_num_threads(4)
    data, table, graphs = load_inputs(args.root)
    jobs = registered_jobs(config)
    (output / f"registered_jobs_worker{args.worker}.json").write_text(json.dumps(jobs, indent=2)+"\n")
    if args.phase == "fit":
        (output / f"frozen_config_worker{args.worker}.json").write_text(json.dumps(config, indent=2)+"\n")
    else:
        selection = json.loads((output / "benchmarks/selection.json").read_text())
        print(json.dumps({"selection_before_test": selection["selected_model"]}), flush=True)
    for job in jobs:
        if args.only_protocol is not None and job["protocol"] != args.only_protocol:
            continue
        if (args.only_model is not None and job["model"] != args.only_model) or job["model"] == args.skip_model:
            continue
        if job["job_index"] % args.workers == args.worker and args.start_job <= job["job_index"] < args.stop_job:
            print(json.dumps({"starting_job": job, "phase": args.phase}), flush=True)
            if args.phase == "fit":
                fit_job(data, table, graphs, job, config, args.device, output)
            else:
                evaluate_job(data, table, graphs, job, config, args.device, output)


if __name__ == "__main__":
    main()
