"""Registered 12-round retrospective GI discovery with a validation-chosen model.

Run two workers after primary model selection:
  python scripts/run_active.py --worker 0 --workers 2 --device cuda:0
  python scripts/run_active.py --worker 1 --workers 2 --device cuda:1
Then collect completed sequences:
  python scripts/run_active.py --phase summarize

The acquisition code receives predictions and already revealed labels only.
Unqueried pool/calibration/test double-perturbation targets are explicitly NaN
in the PairData passed to fitting and prediction. An evaluation-only oracle
tracks measured-pool hits and fixed-test metrics; neither selects a strategy,
an epoch or the next query. Singleton covariates remain available throughout.
"""

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from run_benchmarks import load_inputs, make_partitions, effective_split, train_one, model_predict
from metagene.evaluation import evaluate_predictions


def active_jobs(config, model):
    jobs = []
    for seed in config["primary_seeds"]:
        for strategy in config["active_learning"]["strategies"]:
            jobs.append({"job_index": len(jobs), "seed": seed, "strategy": strategy,
                         "model": model, "name": f"{strategy}_{seed}"})
    return jobs


def add_core_flag(table, root):
    genes = pd.read_csv(root / "data/evidence/genes_annotation.csv", keep_default_na=False).set_index("gene")
    ma = table.gene_a.map(genes.is_strict_metabolic_enzyme).eq(1)
    mb = table.gene_b.map(genes.is_strict_metabolic_enzyme).eq(1)
    ca = table.gene_a.map(genes.is_direct_histone_enzyme).eq(1)
    cb = table.gene_b.map(genes.is_direct_histone_enzyme).eq(1)
    table["is_core_metabolism_histone"] = (ma & cb) | (mb & ca)
    return table


def oracle_targets(data, table, pool):
    """Evaluation-only indices; these targets never enter an acquisition score."""
    targets = {"all": pool[np.argsort(data.gi_proxy[pool], kind="stable")[:100]]}
    metadata = {"all": {"pool_n": len(pool), "hit_k": len(targets["all"])}}
    for name, column in [("strict", "is_strict_enzyme_chromatin"),
                         ("core", "is_core_metabolism_histone")]:
        subset = pool[table.iloc[pool][column].to_numpy()]
        k = min(100, max(1, round(0.1 * len(subset)))) if len(subset) >= 20 else 0
        targets[name] = subset[np.argsort(data.gi_proxy[subset], kind="stable")[:k]]
        metadata[name] = {"pool_n": len(subset), "hit_k": k,
                          "minimum_pool_n": 20, "role": "secondary; not used for acquisition"}
    return targets, metadata


def revealed_table(data, table, indices, reveal_round):
    identities = table.iloc[indices]
    return pd.DataFrame({
        "rowid": indices, "gene_a": identities.gene_a.to_numpy(), "gene_b": identities.gene_b.to_numpy(),
        "gi_truth": data.gi_proxy[indices], "joint_truth": data.joint[indices],
        "label_source": identities.dataset.to_numpy(), "label_type": identities.label_type.to_numpy(),
        "label_definition": "measured joint minus measured single_a minus measured single_b",
        "reveal_round": reveal_round[indices],
    })


def visible_data(data, observed, validation):
    visible = np.concatenate((observed, validation))
    joint = np.full(len(data.pairs), np.nan, dtype=np.float32)
    gi = np.full(len(data.pairs), np.nan, dtype=np.float32)
    joint[visible] = np.asarray(data.joint)[visible]
    gi[visible] = np.asarray(data.gi_proxy)[visible]
    return replace(data, joint=joint, gi_proxy=gi)


@torch.no_grad()
def query_decision(engine, fitted, data, table, available, observed, graphs,
                   strategy, query_size, settings, device, seed, random_generator):
    """Score on GPU without reading any unqueried double-perturbation label."""
    stochastic = strategy in {"uncertainty", "expected_improvement"}
    prediction = model_predict(
        engine, fitted, data, available, graphs,
        mc_passes=settings["mc_passes"] if stochastic else 1,
        seed=seed, return_samples=strategy == "expected_improvement",
    )
    mean = torch.as_tensor(prediction["gi"], dtype=torch.float32, device=device)
    spread = (torch.as_tensor(prediction["gi_std"], dtype=torch.float32, device=device)
              if stochastic else torch.full_like(mean, torch.nan))
    threshold = None
    if strategy == "random":
        score = torch.rand(len(available), device=device, generator=random_generator)
    elif strategy == "effect":
        score = -mean
    elif strategy == "uncertainty":
        score = spread
    else:
        observed_gi = torch.as_tensor(data.gi_proxy[observed], dtype=torch.float32, device=device)
        tau = torch.quantile(observed_gi, settings["ei_observed_quantile"])
        # Gaussian EI approximation; log space avoids finite-MC zero ties.
        sigma = spread.double().clamp_min(1e-12)
        z = (tau.double() - mean.double()) / sigma
        log_phi = -0.5 * z.square() - 0.5 * math.log(2 * math.pi)
        direct = torch.exp(log_phi) + z * torch.special.ndtr(z)
        log_h = torch.log(direct.clamp_min(1e-300))
        negative = z < -10
        u = (-z[negative]).clamp_min(10)
        inv2 = u.reciprocal().square()
        log_h[negative] = log_phi[negative] - 2 * torch.log(u) + torch.log1p(-3 * inv2 + 15 * inv2.square() - 105 * inv2.pow(3))
        score = torch.log(sigma) + log_h
        deterministic = spread.le(1e-12)
        score[deterministic] = torch.log((tau.double() - mean.double())[deterministic].clamp_min(1e-300))
        threshold = float(tau)
    tie_order = torch.argsort(mean, stable=True) if strategy == 'expected_improvement' else torch.arange(len(score), device=device)
    chosen_positions = tie_order[torch.argsort(score[tie_order], descending=True, stable=True)[:query_size]].cpu().numpy()
    selected = available[chosen_positions]
    ranks = np.zeros(len(available), dtype=np.int64)
    ranks[chosen_positions] = np.arange(1, len(selected) + 1)
    identities = table.iloc[available]
    decisions = pd.DataFrame({
        "rowid": available, "gene_a": identities.gene_a.to_numpy(), "gene_b": identities.gene_b.to_numpy(),
        "acquisition_score": score.cpu().numpy(), "pred_gi": mean.cpu().numpy(),
        "spread": spread.cpu().numpy(), "selected_rank": ranks,
        "acquisition_units": 'log_Gaussian_EI' if strategy == 'expected_improvement' else strategy,
    })
    return selected, decisions, threshold


def hit_counts(observed, targets, metadata):
    hits = len(np.intersect1d(observed, targets["all"]))
    result = {"discovered_top100": hits, "recall_top100": hits / len(targets["all"])}
    for name in ("strict", "core"):
        k = metadata[name]["hit_k"]
        n_hits = len(np.intersect1d(observed, targets[name]))
        result.update({f"{name}_pool_n": metadata[name]["pool_n"], f"{name}_hit_k": k,
                       f"{name}_discovered_topk": n_hits if k else None,
                       f"{name}_recall_topk": n_hits / k if k else None})
    return result


def run_sequence(data, table, graphs, job, config, selection, device, output_root):
    seed, strategy, model = job["seed"], job["strategy"], job["model"]
    settings = config["active_learning"]
    output = output_root / "runs" / job["name"]
    output.mkdir(parents=True, exist_ok=True)
    parts = make_partitions(data, table, "pair_holdout", seed, config)
    query_identities = table[["gene_a", "gene_b"]]
    pool = np.asarray(parts["train"])
    initial_size = round(len(pool) * settings["initial_fraction"])
    requested_query_size = round(len(pool) * settings["query_fraction"])
    observed = np.random.default_rng(seed).permutation(pool)[:initial_size]
    initial = observed.copy()
    reveal_round = np.full(len(table), -1, dtype=np.int64)
    reveal_round[initial] = 0
    random_generator = torch.Generator(device=device).manual_seed(seed)
    targets, target_metadata = oracle_targets(data, table, pool)
    np.savez_compressed(output / "partition_indices.npz", train=pool, val=parts["val"],
                        calibration=parts["calibration"], test=parts["test"], initial=initial)
    np.savez_compressed(output / "monitor_oracle_indices.npz", **targets)
    (output / "monitor_oracle_definition.json").write_text(json.dumps({
        "targets": target_metadata, "access": "evaluation only; not passed to acquisition",
        "secondary_hit_rule": "Lowest GI k=min(100,max(1,round(0.1*n_subset))); require n_subset>=20",
    }, indent=2) + "\n")
    revealed_table(data, table, initial, reveal_round).to_csv(output / "initial_revealed_labels.csv", index=False)
    metadata = {
        **job, "completed": False, "device": device,
        "gpu": torch.cuda.get_device_name(torch.device(device)),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "selection": selection, "pool_size": len(pool), "initial_size": initial_size,
        "query_size": requested_query_size, "query_rounds": settings["query_rounds"],
        "validation_size": len(parts["val"]), "calibration_size": len(parts["calibration"]),
        "test_size": len(parts["test"]), "gi_target": "measured additive_gi_proxy",
        "label_source": "ht29_colo1 public measured double-perturbation table",
        "unknown_double_targets_masked": True, "singleton_covariates_available": True,
        "test_metrics_used_for_training_or_acquisition": False,
        "calibration_labels_used": False, "gene_frequency_penalty": None,
        "query_decision_schema": {
            "selected_rank": "1..query_size for selected pairs; zero means unselected",
            "spread": "MC dropout standard deviation; NaN when not computed for random/effect",
            "truth": "No truth in query_decisions; selected truth appears only in revealed_labels",
        },
    }
    (output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    curve, round_records, incoming_query_size = [], [], 0
    for round_number in range(settings["query_rounds"] + 1):
        round_output = output / f"round_{round_number:02d}"
        round_output.mkdir(parents=True, exist_ok=True)
        current_data = visible_data(data, observed, parts["val"])
        current_parts = {**parts, "train": observed.copy(), "pooled_train": observed.copy()}
        split = effective_split(current_parts, model, "pair_holdout", seed + round_number)
        engine, fitted = train_one(
            current_data, split, model, config, device, graphs, round_output,
            progress=lambda row: print(json.dumps({"job": job["name"], "round": round_number, **row}), flush=True),
        )
        test_prediction = model_predict(engine, fitted, current_data, parts["test"], graphs)
        test_metrics = evaluate_predictions(data.gi_proxy[parts["test"]], test_prediction["gi"],
                                  data.pairs[parts["test"]], device, seed + round_number, point_only=True)
        (round_output / "holdout_metrics.json").write_text(json.dumps(test_metrics, indent=2) + "\n")
        revealed_table(data, table, observed, reveal_round).to_csv(round_output / "observed_labels.csv", index=False)
        row = {
            "strategy": strategy, "seed": seed, "round": round_number, "model": model,
            "n_observed": len(observed), **hit_counts(observed, targets, target_metadata),
            "holdout_spearman": test_metrics["spearman"], "holdout_mae": test_metrics["mae"],
            "holdout_pearson": test_metrics["pearson"], "holdout_rmse": test_metrics["rmse"],
            "val_size": len(parts["val"]), "cal_size": len(parts["calibration"]),
            "test_size": len(parts["test"]), "query_size": incoming_query_size,
            "pool_size": len(pool), "initial_size": initial_size,
            "labels_observed_plus_validation": len(observed) + len(parts["val"]),
            "best_epoch": fitted.best_epoch, "elapsed_seconds": fitted.elapsed_seconds,
        }
        available = pool[~np.isin(pool, observed)]
        selected = np.empty(0, dtype=pool.dtype)
        threshold = None
        decisions = pd.DataFrame(columns=["rowid", "gene_a", "gene_b", "acquisition_score",
                                          "pred_gi", "spread", "selected_rank"])
        if round_number < settings["query_rounds"]:
            selected, decisions, threshold = query_decision(
                engine, fitted, current_data, query_identities, available, observed, graphs, strategy,
                min(requested_query_size, len(available)), settings, device, seed + round_number,
                random_generator,
            )
        # Persist the complete decision before accessing selected oracle labels.
        decisions.to_csv(round_output / "query_decisions.csv", index=False)
        reveal_round[selected] = round_number + 1
        revealed_table(data, table, selected, reveal_round).to_csv(round_output / "revealed_labels.csv", index=False)
        np.savez_compressed(round_output / "round_indices.npz", observed=observed, query=selected)
        round_records.append({
            "round": round_number, "initialization_seed": seed + round_number,
            "observed_indices": observed.tolist(), "query_indices": selected.tolist(),
            "query_reveal_round": round_number + 1 if len(selected) else None,
            "label_source": "ht29_colo1 measured additive_gi_proxy; revealed after query ranking",
            "ei_threshold": threshold, "holdout_metrics_role": "monitor only",
        })
        (output / "round_records.json").write_text(json.dumps(round_records, indent=2) + "\n")
        row["next_query_size"] = len(selected)
        row["ei_threshold"] = threshold
        curve.append(row)
        pd.DataFrame(curve).to_csv(output / "learning_curve.csv", index=False)
        print(json.dumps({"job": job["name"], **row}), flush=True)
        observed = np.concatenate((observed, selected))
        incoming_query_size = len(selected)
        del fitted
    metadata.update(completed=True, completed_utc=datetime.now(timezone.utc).isoformat(),
                    final_observed=len(observed), rounds_saved=len(curve))
    (output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"completed_job": job["name"], "final_observed": len(observed)}), flush=True)


def summarize(output):
    completed = [path.parent for path in sorted((output / "runs").glob("*/run.json"))
                 if json.loads(path.read_text())["completed"]]
    curves = pd.concat([pd.read_csv(path / "learning_curve.csv") for path in completed], ignore_index=True)
    curves = curves.sort_values(["strategy", "seed", "round"])
    curves.to_csv(output / "curves.csv", index=False)
    metrics = ["n_observed", "discovered_top100", "recall_top100", "holdout_spearman", "holdout_mae",
               "strict_discovered_topk", "strict_recall_topk", "core_discovered_topk", "core_recall_topk"]
    summary = curves.groupby(["strategy", "round"])[metrics].agg(["mean", "std", "count"])
    summary.columns = [f"{metric}_{statistic}" for metric, statistic in summary.columns]
    summary.reset_index().to_csv(output / "summary.csv", index=False)
    print(json.dumps({"completed_sequences": len(completed), "curve_rows": len(curves),
                      "curves": str(output / "curves.csv")}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--phase", choices=["run", "summarize"], default="run")
    parser.add_argument("--worker", type=int, default=0)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--start-job", type=int, default=0)
    parser.add_argument("--stop-job", type=int, default=100000)
    args = parser.parse_args()
    output = args.root / "results/active_learning"
    output.mkdir(parents=True, exist_ok=True)
    if args.phase == "summarize":
        summarize(output)
        return
    config = json.loads((args.root / "configs/experiment.json").read_text())
    selection = json.loads((args.root / "results/benchmarks/selection.json").read_text())
    model = selection["selected_al_model"]
    jobs = active_jobs(config, model)
    (output / f"registered_jobs_worker{args.worker}.json").write_text(json.dumps(jobs, indent=2) + "\n")
    (output / f"frozen_config_worker{args.worker}.json").write_text(json.dumps(config, indent=2) + "\n")
    torch.set_num_threads(4)
    data, table, graphs = load_inputs(args.root)
    table = add_core_flag(table, args.root)
    for job in jobs:
        if job["job_index"] % args.workers == args.worker and args.start_job <= job["job_index"] < args.stop_job:
            print(json.dumps({"starting_job": job, "selected_al_model": model}), flush=True)
            run_sequence(data, table, graphs, job, config, selection, args.device, output)


if __name__ == "__main__":
    main()
