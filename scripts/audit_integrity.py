"""Audit completed scientific artifacts without fitting or executing models.

Read the registered benchmark, final-fit, active-learning and candidate
artifacts; write only results/integrity_audit.json. Missing
artifacts fail at their ordinary read operation. This is an audit of actual
delivery records, not a unit/smoke-test suite or a claim of biological validity.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import itertools
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RTOL = 1e-6
ATOL = 1e-7


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def add(checks, section, name, passed, **details):
    checks.append({"section": section, "check": name, "passed": bool(passed), **details})


def same_rows(checks, section, name, left, right):
    add(checks, section, name, np.array_equal(left, right),
        left_n=len(left), right_n=len(right), order_sensitive=True)


def disjoint(checks, section, name, left, right):
    overlap = np.intersect1d(left, right)
    add(checks, section, name, len(overlap) == 0,
        overlap_n=len(overlap), example_ids=overlap[:10].tolist())


def labels_match(checks, section, name, frame, indices, raw, expected_rounds):
    same_rows(checks, section, f"{name}/rowids", frame.rowid.to_numpy(), indices)
    source = raw.loc[frame.rowid.to_numpy()]
    add(checks, section, f"{name}/measured_gi",
        np.allclose(frame.gi_truth.to_numpy(dtype=float), source.additive_gi_proxy.to_numpy(dtype=float),
                    rtol=RTOL, atol=ATOL))
    add(checks, section, f"{name}/measured_joint",
        np.allclose(frame.joint_truth.to_numpy(dtype=float), source.joint_effect.to_numpy(dtype=float),
                    rtol=RTOL, atol=ATOL))
    add(checks, section, f"{name}/identity_and_source",
        np.array_equal(frame.gene_a.to_numpy(), source.gene_a.to_numpy())
        and np.array_equal(frame.gene_b.to_numpy(), source.gene_b.to_numpy())
        and np.array_equal(frame.label_source.to_numpy(), source.dataset.to_numpy())
        and np.array_equal(frame.label_type.to_numpy(), source.label_type.to_numpy()))
    same_rows(checks, section, f"{name}/reveal_round",
              frame.reveal_round.to_numpy(), expected_rounds[indices])


def registration(config):
    jobs = []
    for protocol, models in config["protocols"].items():
        seeds = config["primary_seeds"] if protocol in {"pair_holdout", "gene_holdout"} else config["secondary_seeds"]
        for seed in seeds:
            for model in models:
                jobs.append({"protocol": protocol, "seed": seed, "model": model,
                             "name": f"{protocol}_{model}_{seed}", "job_index": len(jobs)})
    return jobs


def audit_benchmarks(root, config, pair_codes, checks):
    output = root / "results"
    directory = output / "benchmarks/runs"
    jobs = registration(config)
    add(checks, "benchmarks", "registered_run_count", len(jobs) == sum(len(m) * (len(config["primary_seeds"]) if p in {"pair_holdout","gene_holdout"} else len(config["secondary_seeds"])) for p,m in config["protocols"].items()), expected=len(jobs), actual=len(jobs))
    for worker in (0, 1):
        saved = read_json(output / f"registered_jobs_worker{worker}.json")
        add(checks, "benchmarks", f"worker_{worker}_registration", saved == jobs, actual=len(saved))
    saved_names = sorted(path.parent.name for path in directory.glob("*/run.json"))
    expected_names = sorted(job["name"] for job in jobs)
    add(checks, "benchmarks", "artifact_run_names", saved_names == expected_names,
        actual=len(saved_names), expected=len(expected_names),
        unexpected=sorted(set(saved_names) - set(expected_names)))
    partitions, metadata, protocol_reference, model_reference = {}, [], {}, {}
    for job in jobs:
        name = job["name"]
        run = directory / name
        meta = read_json(run / "run.json")
        parts = dict(np.load(run / "partition_indices.npz"))
        partitions[name] = parts
        metadata.append(meta)
        add(checks, "benchmarks", f"{name}/completed_cuda",
            meta["completed"] and meta["device"].startswith("cuda") and bool(meta["gpu"]),
            device=meta["device"], gpu=meta["gpu"])
        add(checks, "benchmarks", f"{name}/registered_identity",
            all(meta[key] == job[key] for key in ("name", "model", "seed", "protocol", "job_index")))
        add(checks, "benchmarks", f"{name}/test_not_used_during_fit", meta["test_evaluated_during_fit"] is False)
        reference = protocol_reference.setdefault(job["protocol"], parts)
        for role in ("train", "val", "calibration", "test", "pooled_train"):
            same_rows(checks, "fixed_partitions", f"{name}/{role}", parts[role], reference[role])
            add(checks, "fixed_partitions", f"{name}/{role}_unique", len(np.unique(parts[role])) == len(parts[role]))
        key = (job["protocol"], job["model"])
        reference_train = model_reference.setdefault(key, parts["effective_train"])
        same_rows(checks, "fixed_partitions", f"{name}/effective_train_across_seeds",
                  parts["effective_train"], reference_train)
        train_role = "pooled_train" if job["model"].startswith("pooled_") or job["protocol"] == "transfer_ht29" else "train"
        same_rows(checks, "benchmarks", f"{name}/effective_training_role", parts["effective_train"], parts[train_role])
        if job["protocol"] == "pair_holdout":
            for role in ("val", "calibration", "test"):
                disjoint(checks, "primary_separation", f"{name}/train_vs_{role}_records",
                         parts["effective_train"], parts[role])
                disjoint(checks, "primary_separation", f"{name}/train_vs_{role}_pairs",
                         pair_codes[parts["effective_train"]], pair_codes[parts[role]])
    primary_jobs = [job for job in jobs if job["protocol"] == "pair_holdout"]
    primary = partitions[primary_jobs[0]["name"]]
    for a, b in itertools.combinations(("train", "val", "calibration", "test"), 2):
        disjoint(checks, "primary_separation", f"base_{a}_vs_{b}_records", primary[a], primary[b])
        disjoint(checks, "primary_separation", f"base_{a}_vs_{b}_pairs", pair_codes[primary[a]], pair_codes[primary[b]])
    development = np.unique(np.concatenate([
        partitions[job["name"]][role] for job in primary_jobs for role in ("effective_train", "val")
    ]))
    reserved = np.unique(np.concatenate([
        partitions[job["name"]][role] for job in primary_jobs for role in ("test", "calibration")
    ]))
    disjoint(checks, "primary_separation", "union_development_vs_reserved_records", development, reserved)
    disjoint(checks, "primary_separation", "union_development_vs_reserved_pairs", pair_codes[development], pair_codes[reserved])
    selection = read_json(output / "benchmarks/selection.json")
    fit_table = pd.DataFrame(metadata)
    scores = fit_table[fit_table.protocol.eq("pair_holdout")].groupby("model").validation_rmse.mean()
    add(checks, "selection", "primary_validation_choice", selection["selected_model"] == scores.loc[["target_string","pooled_string","pooled_reactome","pooled_multigraph"]].idxmin(),
        selected=selection["selected_model"], recomputed=scores.loc[["target_string","pooled_string","pooled_reactome","pooled_multigraph"]].idxmin())
    al_choice = scores.loc[["target_graphless", "target_string"]].idxmin()
    add(checks, "selection", "target_only_al_choice", selection["selected_al_model"] == al_choice,
        selected=selection["selected_al_model"], recomputed=al_choice)
    add(checks, "selection", "test_metrics_not_used", selection["test_metrics_used"] is False)
    return primary, selection


def audit_final(root, config, primary, selection, raw, pair_codes, checks):
    output = root / "results"
    model = selection["selected_model"]
    expected_target_train = np.concatenate((primary["train"], primary["test"]))
    reference, calibration_tables = None, []
    for seed in config["primary_seeds"]:
        name = f"{model}_{seed}"
        directory = output / "final_models" / name
        meta = read_json(directory / "run.json")
        parts = dict(np.load(directory / "partition_indices.npz"))
        cal = pd.read_csv(directory / "calibration_predictions.csv")
        calibration_tables.append(cal)
        add(checks, "final_models", f"{name}/completed_cuda",
            meta["completed"] and meta["device"].startswith("cuda"), device=meta["device"])
        add(checks, "final_models", f"{name}/selected_model_and_partition",
            meta["model"] == model and meta["seed"] == seed and meta["partition_seed"] == config["partition_seed"])
        same_rows(checks, "final_models", f"{name}/calibration_is_primary", parts["calibration"], primary["calibration"])
        same_rows(checks, "final_models", f"{name}/validation_is_primary", parts["val"], primary["val"])
        same_rows(checks, "final_models", f"{name}/target_train_is_primary_train_plus_test", parts["train"], expected_target_train)
        add(checks, "final_models", f"{name}/no_reused_test_partition", len(parts["test"]) == 0)
        train = parts["effective_train"]
        for role in ("val", "calibration"):
            disjoint(checks, "final_models", f"{name}/train_vs_{role}_records", train, parts[role])
            disjoint(checks, "final_models", f"{name}/train_vs_{role}_pairs", pair_codes[train], pair_codes[parts[role]])
        disjoint(checks, "final_models", f"{name}/val_vs_cal_records", parts["val"], parts["calibration"])
        disjoint(checks, "final_models", f"{name}/val_vs_cal_pairs", pair_codes[parts["val"]], pair_codes[parts["calibration"]])
        same_rows(checks, "final_models", f"{name}/calibration_prediction_rows", cal.row_id.to_numpy(), primary["calibration"])
        add(checks, "final_models", f"{name}/calibration_measured_labels",
            np.allclose(cal.observed_gi, raw.loc[cal.row_id, "additive_gi_proxy"], rtol=RTOL, atol=ATOL))
        if reference is None:
            reference = parts
        for role in ("train", "effective_train", "val", "calibration", "test"):
            same_rows(checks, "final_models", f"{name}/{role}_across_seeds", parts[role], reference[role])
    add(checks, "final_models", "five_final_fits", len(calibration_tables) == 5, actual=len(calibration_tables))
    ensemble = pd.read_csv(output / "candidates/ensemble_calibration.csv")
    same_rows(checks, "ensemble_calibration", "primary_calibration_rows", ensemble.row_id.to_numpy(), primary["calibration"])
    add(checks, "ensemble_calibration", "measured_labels",
        np.allclose(ensemble.observed_gi, raw.loc[ensemble.row_id, "additive_gi_proxy"], rtol=RTOL, atol=ATOL))
    expected = np.mean([frame.pred_gi.to_numpy() for frame in calibration_tables], axis=0)
    add(checks, "ensemble_calibration", "mean_of_five_frozen_predictions",
        np.allclose(ensemble.ensemble_pred_gi, expected, rtol=RTOL, atol=ATOL))
    residual = np.sort(np.abs(ensemble.observed_gi.to_numpy() - ensemble.ensemble_pred_gi.to_numpy()))
    rank = math.ceil((len(residual) + 1) * config["uncertainty"]["prediction_interval_level"])
    radius = float(np.append(residual, np.inf)[rank - 1])
    prediction_meta = read_json(output / "candidates/prediction_metadata.json")
    add(checks, "ensemble_calibration", "finite_sample_radius",
        np.isclose(prediction_meta["calibration_radius"], radius, rtol=RTOL, atol=ATOL),
        reported=prediction_meta["calibration_radius"], recomputed=radius, order_statistic=rank)
    return radius


def audit_active(root, config, primary, selection, raw, checks):
    directory = root / "results/active_learning"
    seeds, strategies = config["primary_seeds"], config["active_learning"]["strategies"]
    expected_names = sorted(f"{strategy}_{seed}" for seed in seeds for strategy in strategies)
    saved_names = sorted(path.parent.name for path in (directory / "runs").glob("*/run.json"))
    add(checks, "active_learning", "twenty_sequences", len(expected_names) == 20 and saved_names == expected_names,
        expected=len(expected_names), actual=len(saved_names), unexpected=sorted(set(saved_names) - set(expected_names)))
    initial_by_seed, all_curves = {}, []
    decision_columns = {"rowid", "gene_a", "gene_b", "acquisition_score", "pred_gi", "spread", "selected_rank"}
    n_rounds = config["active_learning"]["query_rounds"] + 1
    for seed in seeds:
        for strategy in strategies:
            name = f"{strategy}_{seed}"
            run = directory / "runs" / name
            meta = read_json(run / "run.json")
            parts = dict(np.load(run / "partition_indices.npz"))
            curve = pd.read_csv(run / "learning_curve.csv")
            records = read_json(run / "round_records.json")
            all_curves.append(curve)
            add(checks, "active_learning", f"{name}/completed_cuda_model",
                meta["completed"] and meta["device"].startswith("cuda") and bool(meta["gpu"])
                and meta["model"] == selection["selected_al_model"])
            add(checks, "active_learning", f"{name}/thirteen_rounds",
                len(curve) == len(records) == n_rounds == 13 and curve["round"].tolist() == list(range(n_rounds)),
                curve_rows=len(curve), record_rows=len(records))
            add(checks, "active_learning", f"{name}/label_access_policy",
                meta["unknown_double_targets_masked"] is True
                and meta["test_metrics_used_for_training_or_acquisition"] is False
                and meta["calibration_labels_used"] is False)
            for role in ("train", "val", "calibration", "test"):
                same_rows(checks, "active_learning", f"{name}/{role}_is_primary", parts[role], primary[role])
            pool, initial = parts["train"], parts["initial"]
            reference_initial = initial_by_seed.setdefault(seed, initial)
            same_rows(checks, "active_learning", f"{name}/same_initial_across_strategies", initial, reference_initial)
            add(checks, "active_learning", f"{name}/initial_budget",
                len(initial) == round(len(pool) * config["active_learning"]["initial_fraction"]))
            for role in ("val", "calibration", "test"):
                disjoint(checks, "active_learning", f"{name}/pool_vs_{role}", pool, parts[role])
            reveal_round = np.full(len(raw), -1, dtype=np.int64)
            reveal_round[initial] = 0
            labels_match(checks, "active_labels", f"{name}/initial",
                         pd.read_csv(run / "initial_revealed_labels.csv"), initial, raw, reveal_round)
            expected_observed = initial.copy()
            previous_query = np.empty(0, dtype=np.int64)
            for round_number in range(n_rounds):
                prefix = f"{name}/round_{round_number:02d}"
                current = run / f"round_{round_number:02d}"
                indices = dict(np.load(current / "round_indices.npz"))
                observed, query = indices["observed"], indices["query"]
                entry = records[round_number]
                same_rows(checks, "active_indices", f"{prefix}/observed_growth", observed, expected_observed)
                same_rows(checks, "active_indices", f"{prefix}/observed_json", observed, np.asarray(entry["observed_indices"]))
                same_rows(checks, "active_indices", f"{prefix}/query_json", query, np.asarray(entry["query_indices"]))
                add(checks, "active_indices", f"{prefix}/seed_and_round",
                    entry["round"] == round_number and entry["initialization_seed"] == seed + round_number)
                add(checks, "active_indices", f"{prefix}/unique_observations_and_queries",
                    len(np.unique(observed)) == len(observed) and len(np.unique(query)) == len(query))
                disjoint(checks, "active_indices", f"{prefix}/query_not_observed", query, observed)
                add(checks, "active_indices", f"{prefix}/pool_membership",
                    np.isin(observed, pool).all() and np.isin(query, pool).all())
                for role in ("val", "calibration", "test"):
                    disjoint(checks, "active_indices", f"{prefix}/observed_vs_{role}", observed, parts[role])
                    disjoint(checks, "active_indices", f"{prefix}/query_vs_{role}", query, parts[role])
                row = curve.iloc[round_number]
                expected_query_n = (min(round(len(pool) * config["active_learning"]["query_fraction"]), len(pool) - len(observed))
                                    if round_number < n_rounds - 1 else 0)
                add(checks, "active_indices", f"{prefix}/budget_and_curve",
                    len(query) == expected_query_n and row.n_observed == len(observed)
                    and row.query_size == len(previous_query) and row.next_query_size == len(query)
                    and row.val_size == len(parts["val"]) and row.cal_size == len(parts["calibration"])
                    and row.test_size == len(parts["test"]))
                decisions = pd.read_csv(current / "query_decisions.csv")
                expected_columns = decision_columns | {'acquisition_units'} if round_number < n_rounds - 1 else decision_columns
                add(checks, "active_decisions", f"{prefix}/prediction_only_columns",
                    set(decisions.columns) == expected_columns, columns=decisions.columns.tolist())
                if round_number < n_rounds - 1:
                    expected_units = 'log_Gaussian_EI' if strategy == 'expected_improvement' else strategy
                    add(checks, 'active_decisions', f'{prefix}/declared_score_units',
                        decisions.acquisition_units.eq(expected_units).all(), expected=expected_units)
                available = pool[~np.isin(pool, observed)] if round_number < n_rounds - 1 else np.empty(0, dtype=np.int64)
                same_rows(checks, "active_decisions", f"{prefix}/available_pool", decisions.rowid.to_numpy(), available)
                selected = decisions[decisions.selected_rank > 0].sort_values("selected_rank")
                same_rows(checks, "active_decisions", f"{prefix}/selected_query", selected.rowid.to_numpy(), query)
                same_rows(checks, "active_decisions", f"{prefix}/selected_rank", selected.selected_rank.to_numpy(), np.arange(1, len(query) + 1))
                labels_match(checks, "active_labels", f"{prefix}/observed",
                             pd.read_csv(current / "observed_labels.csv"), observed, raw, reveal_round)
                reveal_round[query] = round_number + 1
                labels_match(checks, "active_labels", f"{prefix}/revealed",
                             pd.read_csv(current / "revealed_labels.csv"), query, raw, reveal_round)
                expected_observed = np.concatenate((observed, query))
                previous_query = query
            add(checks, "active_learning", f"{name}/final_budget",
                meta["final_observed"] == len(expected_observed) == int(curve.iloc[-1].n_observed))
    merged = pd.read_csv(directory / "curves.csv").sort_values(["strategy", "seed", "round"]).reset_index(drop=True)
    per_run = pd.concat(all_curves, ignore_index=True).sort_values(["strategy", "seed", "round"]).reset_index(drop=True)
    add(checks, "active_learning", "merged_curve_keys",
        merged[["strategy", "seed", "round"]].equals(per_run[["strategy", "seed", "round"]]),
        expected_rows=20 * 13, actual_rows=len(merged))
    for column in ("n_observed", "discovered_top100", "recall_top100", "holdout_spearman", "holdout_mae"):
        add(checks, "active_learning", f"merged_curve_{column}",
            np.allclose(merged[column], per_run[column], rtol=RTOL, atol=ATOL, equal_nan=True))


def canonical_keys(frame):
    return pd.Series(["|".join(sorted(pair)) for pair in frame[["gene_a", "gene_b"]].to_numpy()], index=frame.index)


def audit_candidates(root, config, selection, raw, radius, checks):
    directory = root / "results/candidates"
    predicted = pd.read_csv(directory / "predicted_pairs.csv")
    core = pd.read_csv(directory / "core_candidates.csv")
    top = pd.read_csv(directory / "top30.csv")
    priority = pd.read_csv(directory / "priority8.csv")
    ranking = read_json(directory / "ranking_metadata.json")
    annotation = pd.read_csv(root / "data/evidence/genes_annotation.csv").set_index("gene")
    measured = set(canonical_keys(raw[raw.dataset.eq("ht29_colo1")]))
    metabolism = annotation.index[annotation.is_strict_metabolic_enzyme.eq(1)]
    histone = annotation.index[annotation.is_direct_histone_enzyme.eq(1)]
    expected_core = {'|'.join(sorted((a,b))) for a in metabolism for b in histone if a != b} - measured
    top500 = pd.read_csv(directory / 'top500.csv')
    representative = pd.read_csv(directory / 'representative3.csv')
    add(checks, 'candidates', 'full_core_space_reconstructed', set(core.pair_key) == expected_core)
    for name, frame, expected_n in (("strict_core", core, len(expected_core)), ('top500', top500, 500),
                                   ("top30", top, 30), ("priority8", priority, 8), ('representative3', representative, 3)):
        keys = canonical_keys(frame)
        add(checks, "candidates", f"{name}/count", len(frame) == expected_n, expected=expected_n, actual=len(frame))
        add(checks, "candidates", f"{name}/unique_unmeasured_distinct_pairs",
            keys.nunique() == len(frame) and not keys.isin(measured).any() and frame.gene_a.ne(frame.gene_b).all())
        ma = frame.gene_a.map(annotation.is_strict_metabolic_enzyme).eq(1)
        mb = frame.gene_b.map(annotation.is_strict_metabolic_enzyme).eq(1)
        ca = frame.gene_a.map(annotation.is_direct_histone_enzyme).eq(1)
        cb = frame.gene_b.map(annotation.is_direct_histone_enzyme).eq(1)
        add(checks, "candidates", f"{name}/strict_metabolic_core_chromatin", ((ma & cb) | (mb & ca)).all())
        add(checks, "candidates", f"{name}/canonical_pair_key", np.array_equal(keys.to_numpy(), frame.pair_key.to_numpy()))
        if name in {"top30", "priority8"}:
            counts = Counter(frame.gene_a.tolist() + frame.gene_b.tolist())
            maximum = max(counts.values(), default=0)
            add(checks, "candidates", f"{name}/gene_frequency_cap", maximum <= config["candidate_policy"]["maximum_pairs_per_gene"],
                maximum=maximum, cap=config["candidate_policy"]["maximum_pairs_per_gene"])
    add(checks, "candidates", "priority_subset_top30_subset_core",
        set(priority.pair_key) <= set(top.pair_key) <= set(core.pair_key) <= set(predicted.pair_key))
    add(checks, "candidates", "ranking_metadata_counts",
        ranking["n_strict_core_unmeasured"] == len(core) and ranking["n_top30"] == len(top)
        and ranking["n_priority8"] == len(priority))
    for name, frame in (("predicted", predicted), ("strict_core", core), ("top30", top), ("priority8", priority)):
        expected_joint = frame.single_a.to_numpy() + frame.single_b.to_numpy() + frame.pred_gi.to_numpy()
        gap = np.abs(frame.pred_joint.to_numpy() - expected_joint)
        add(checks, "candidate_predictions", f"{name}/joint_reconstruction",
            np.allclose(frame.pred_joint, expected_joint, rtol=RTOL, atol=ATOL, equal_nan=True), maximum_absolute_gap=float(np.nanmax(gap)))
        add(checks, "candidate_predictions", f"{name}/selected_ensemble",
            frame.model.eq(selection["selected_model"]).all() and frame.n_models.eq(len(config["primary_seeds"])).all())
        add(checks, "candidate_predictions", f"{name}/nominal_calibration_band",
            np.allclose(frame.pi90_low, frame.pred_gi - radius, rtol=RTOL, atol=ATOL)
            and np.allclose(frame.pi90_high, frame.pred_gi + radius, rtol=RTOL, atol=ATOL))
    add(checks, 'candidates', 'top500_exact_score_order', top500.pair_key.tolist() == core.sort_values(['pred_gi','pair_key']).head(500).pair_key.tolist())
    add(checks, 'candidates', 'shortlist_known_single_effects', top[['single_a','single_b']].notna().all().all())
    add(checks, 'candidates', 'representatives_diverse', pd.concat([representative.gene_a,representative.gene_b]).nunique() == 6)
    add(checks, 'handoff', 'gene_pair_manifest_count', len(pd.read_csv(root/'results/handoff/mini_pool_gene_pairs.csv')) == 500)
    add(checks, 'handoff', 'capacity_slots_count', len(pd.read_csv(root/'results/handoff/library_capacity_plan.csv')) == 1000)
    index = pd.read_csv(directory / "evidence_card_index.csv")
    card_paths = sorted((directory / "evidence_cards").glob("*.md"))
    add(checks, "evidence_cards", "thirty_indexed_and_physical_cards", len(index) == len(card_paths) == 30,
        index_n=len(index), physical_n=len(card_paths))
    add(checks, "evidence_cards", "one_card_per_top30_pair",
        index.pair_key.nunique() == 30 and set(index.pair_key) == set(top.pair_key))
    add(checks, "evidence_cards", "no_stale_extra_files",
        {str(path.relative_to(directory)) for path in card_paths} == set(index.file))
    cards = []
    top_by_key = top.set_index("pair_key")
    for row in index.itertuples(index=False):
        text = (directory / row.file).read_text(encoding="utf-8")
        cards.append(text)
        pair = top_by_key.loc[row.pair_key]
        header = f"# {int(pair.top30_rank):02d}. {pair.metabolic_gene} + {pair.chromatin_gene}"
        add(checks, "evidence_cards", f"{row.pair_key}/identity", text.splitlines()[0] == header)
    combined = (directory / "candidate_evidence_cards.md").read_text(encoding="utf-8")
    add(checks, "evidence_cards", "combined_cards_match_individual_cards", combined == "\n\n---\n\n".join(cards) + "\n")


def audit_independent_inputs(root, checks):
    evidence = root / 'data/evidence'
    manifest = read_json(evidence / 'integrated_input_manifest.json')
    features = np.load(evidence / 'integrated_features.npz')['features']
    names = pd.read_csv(evidence / 'integrated_feature_names.csv')
    genes = pd.read_csv(root / 'data/processed/model_genes.csv').gene.tolist()
    values = pd.read_csv(evidence / 'integrated_gene_evidence.csv', index_col='gene')
    add(checks, 'independent_inputs', 'gene_order', values.index.tolist() == genes)
    add(checks, 'independent_inputs', 'feature_dimensions', features.shape == (len(genes), len(names))
        and manifest['features'] == features.shape[1])
    add(checks, 'independent_inputs', 'no_pair_labels', manifest['pair_label_use'] is False)
    n = len(values.columns)
    add(checks, 'independent_inputs', 'raw_covariates_preserved',
        np.allclose(features[:, -2*n:-n], values.to_numpy(float), equal_nan=True, rtol=RTOL, atol=ATOL))
    add(checks, 'independent_inputs', 'explicit_missing_indicators',
        np.array_equal(features[:, -n:], (~np.isfinite(values.to_numpy(float))).astype(float)))
    for name, edges in dict(np.load(evidence / 'graphs.npz')).items():
        add(checks, 'independent_inputs', f'{name}/valid_gene_indices',
            edges.shape[0] == 2 and np.isfinite(edges).all() and edges.min() >= 0 and edges.max() < len(genes))
    mapped = pd.read_csv(evidence / 'hct116_flux_gene_reaction_map.csv')
    add(checks, 'independent_inputs', 'flux_mapping_count', mapped.gene.nunique() == manifest['metabolic_flux']['covered_genes'])
    peaks = pd.read_csv(evidence / 'histone_gene_peak_evidence.csv')
    add(checks, 'independent_inputs', 'six_histone_tracks', peaks['sample'].nunique() == 6)
    for sidecar in (evidence / 'requested_sources').glob('*.provenance.json'):
        record = read_json(sidecar)
        path = sidecar.parent / record['file']
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(2**20), b''):
                digest.update(block)
        add(checks, 'source_provenance', record['file'],
            path.stat().st_size == record['bytes'] and digest.hexdigest() == record['sha256'])
    expected = read_json(root / 'results/input_fingerprints.json')
    for file, digest in expected.items():
        add(checks, 'frozen_inputs', file, hashlib.sha256((root / file).read_bytes()).hexdigest() == digest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    config = read_json(args.root / "configs/experiment.json")
    raw = pd.read_csv(args.root / "data/processed/model_pairs.csv").set_index("row_id", drop=False)
    arrays = np.load(args.root / "data/processed/model_arrays.npz")
    pairs = np.sort(arrays["pairs"], axis=1)
    pair_codes = pairs[:, 0] * len(arrays["gene_features"]) + pairs[:, 1]
    checks = []
    audit_independent_inputs(args.root, checks)
    same_rows(checks, "source", "measurement_row_order", raw.row_id.to_numpy(), np.arange(len(pairs)))
    primary, selection = audit_benchmarks(args.root, config, pair_codes, checks)
    radius = audit_final(args.root, config, primary, selection, raw, pair_codes, checks)
    audit_active(args.root, config, primary, selection, raw, checks)
    audit_candidates(args.root, config, selection, raw, radius, checks)
    failed = [check for check in checks if not check["passed"]]
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(), "passed": len(failed) == 0,
        "check_count": len(checks), "passed_count": len(checks) - len(failed), "failed_count": len(failed),
        "partition_seed": config["partition_seed"], "primary_seeds": config["primary_seeds"],
        "selected_model": selection["selected_model"], "selected_al_model": selection["selected_al_model"],
        "numeric_tolerance": {"rtol": RTOL, "atol": ATOL, "reason": "float32 model labels serialized through CSV"},
        "scope": "Completed artifact integrity; no model execution, new fitting or biological validation",
        "primary_union_scope": "Only primary runs influence architecture selection; secondary protocols have their own fixed partitions",
        "limitations": "CUDA is verified from runtime metadata. Index separation and label files do not alone prove every implementation access. Shared-gene dependence and candidate distribution shift still limit interval interpretation.",
        "failed_checks": failed, "checks": checks,
    }
    output = args.root / "results/integrity_audit.json"
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": report["passed"], "checks": len(checks), "failed": len(failed),
                      "output": str(output)}, ensure_ascii=False), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    import sys
    if '--legacy-gi' in sys.argv:
        sys.argv.remove('--legacy-gi')
        main()
    else:
        from audit_sl import main as current_main
        current_main()
