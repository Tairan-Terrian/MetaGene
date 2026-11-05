"""Regression/ranking metrics, dyadic gene bootstrap and split conformal bands.

Bootstrap replicates draw an independent Poisson(1) count for each distinct
gene. An observed pair (a,b) receives count w_a*w_b, preserving shared-gene
dependence. Weighted Spearman recomputes exact weighted midranks in every
replicate, including ties; it does not correlate fixed original ranks.

Percentile confidence intervals are conditional on this fitted model and
test collection. They do not include model-training variation, which is
reported separately across seeds. A small number of genes can yield unstable
or degenerate intervals. Zero-total-weight replicates and mathematically
undefined correlations are omitted, with their counts retained in the output.

TopK is evaluated once on the original finite test ranking, selecting lower
values and breaking ties by original row order. It is not bootstrapped by
duplicating pairs. Conformal bands are marginal split-conformal intervals;
their usual guarantee additionally requires exchangeable calibration/test
observations, which is not established merely by using dyadic bootstrap.
"""

import math
import numpy as np
from scipy.stats import rankdata
import torch


def _pearson(x, y):
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    centered_x, centered_y = x - x.mean(), y - y.mean()
    value = np.dot(centered_x, centered_y) / (np.linalg.norm(centered_x) * np.linalg.norm(centered_y))
    return float(np.clip(value, -1.0, 1.0))


def _rank_plan(values):
    order = torch.argsort(values, stable=True)
    sorted_values = values[order]
    starts = torch.cat((torch.ones(1, dtype=torch.bool, device=values.device),
                        sorted_values[1:] != sorted_values[:-1]))
    group_ids = starts.long().cumsum(dim=0) - 1
    ends = torch.cat((torch.nonzero(starts[1:], as_tuple=True)[0],
                      torch.tensor([len(values) - 1], device=values.device)))
    return order, torch.argsort(order), group_ids, ends


def _weighted_midranks(weights, plan):
    order, inverse, group_ids, ends = plan
    cumulative = weights[:, order].cumsum(dim=1)
    group_end_mass = cumulative[:, ends]
    group_start_mass = torch.cat((torch.zeros_like(group_end_mass[:, :1]),
                                  group_end_mass[:, :-1]), dim=1)
    group_mass = group_end_mass - group_start_mass
    # A common +0.5 rank offset is unnecessary for a correlation.
    group_midrank = (group_start_mass + group_end_mass) / 2
    ranks = group_midrank[:, group_ids][:, inverse]
    distinct_supported_values = (group_mass > 0).sum(dim=1)
    return ranks, distinct_supported_values


def _weighted_correlation(x, y, weights):
    total = weights.sum(dim=1)
    centered_x = x - (weights * x).sum(dim=1, keepdim=True) / total[:, None]
    centered_y = y - (weights * y).sum(dim=1, keepdim=True) / total[:, None]
    covariance = (weights * centered_x * centered_y).sum(dim=1)
    variance_x = (weights * centered_x.square()).sum(dim=1)
    variance_y = (weights * centered_y.square()).sum(dim=1)
    return (covariance / (variance_x * variance_y).sqrt()).clamp(-1.0, 1.0)


@torch.no_grad()
def evaluate_predictions(y_true, y_pred, pairs, device, seed, bootstrap_replicates=500,
               point_only=False):
    """Return flat point estimates and 95% dyadic-bootstrap percentile CIs.

    Bootstrap computation uses the requested torch device, without fallback.
    Point estimates use the finite prediction/label pairs. Per-metric valid
    replicate counts distinguish zero-weight draws from constant-support
    correlation draws. point_only=True omits all bootstrap computation.
    """
    truth, prediction = np.asarray(y_true, dtype=np.float64), np.asarray(y_pred, dtype=np.float64)
    measured = np.isfinite(truth) & np.isfinite(prediction)
    truth, prediction = truth[measured], prediction[measured]
    pairs = np.asarray(pairs)[measured]
    genes, gene_rows = np.unique(pairs, return_inverse=True)
    gene_rows = gene_rows.reshape(-1, 2)
    n = len(truth)
    output = {
        "n": n, "n_genes": len(genes), "pearson": None, "spearman": None,
        "rmse": None, "mae": None, "top_k": 0, "top_k_hits": 0,
        "top_k_overlap": None, "top_k_random_expected": None, "top_k_enrichment": None,
        "bootstrap_requested": 0 if point_only else bootstrap_replicates,
        "bootstrap_completed": 0, "bootstrap_nonzero": 0, "bootstrap_zero_weight": 0,
        "ci_level": 0.95, "ci_method": "dyadic_gene_poisson_percentile",
        "bootstrap_spearman": "exact_weighted_midranks",
    }
    for metric in ("pearson", "spearman", "rmse", "mae"):
        output[f"{metric}_ci_low"] = None
        output[f"{metric}_ci_high"] = None
        output[f"{metric}_bootstrap_valid"] = 0
    if n == 0:
        return output
    output["pearson"] = _pearson(truth, prediction)
    output["spearman"] = _pearson(rankdata(truth), rankdata(prediction))
    error = prediction - truth
    output["rmse"] = float(np.sqrt(np.mean(error ** 2)))
    output["mae"] = float(np.mean(np.abs(error)))
    k = min(100, max(1, round(0.1 * n)))
    true_top = np.argsort(truth, kind="stable")[:k]
    predicted_top = np.argsort(prediction, kind="stable")[:k]
    hits = len(np.intersect1d(true_top, predicted_top))
    output.update(top_k=k, top_k_hits=hits, top_k_overlap=hits / k,
                  top_k_random_expected=k / n, top_k_enrichment=(hits / k) / (k / n))
    if point_only or bootstrap_replicates == 0:
        return output

    truth_gpu = torch.as_tensor(truth, dtype=torch.float64, device=device)
    prediction_gpu = torch.as_tensor(prediction, dtype=torch.float64, device=device)
    gene_rows_gpu = torch.as_tensor(gene_rows, dtype=torch.long, device=device)
    error_gpu = prediction_gpu - truth_gpu
    truth_plan, prediction_plan = _rank_plan(truth_gpu), _rank_plan(prediction_gpu)
    generator = torch.Generator(device=device).manual_seed(seed)
    replicates = {metric: [] for metric in ("pearson", "spearman", "rmse", "mae")}
    nonzero_count = torch.zeros((), dtype=torch.long, device=device)
    for start in range(0, bootstrap_replicates, 64):
        batch_size = min(64, bootstrap_replicates - start)
        gene_counts = torch.poisson(
            torch.ones((batch_size, len(genes)), dtype=torch.float64, device=device),
            generator=generator,
        )
        weights = gene_counts[:, gene_rows_gpu[:, 0]] * gene_counts[:, gene_rows_gpu[:, 1]]
        nonzero = weights.sum(dim=1) > 0
        nonzero_count += nonzero.sum()
        weights = weights[nonzero]
        total = weights.sum(dim=1)
        truth_ranks, truth_support = _weighted_midranks(weights, truth_plan)
        prediction_ranks, prediction_support = _weighted_midranks(weights, prediction_plan)
        correlation_defined = (truth_support > 1) & (prediction_support > 1)
        pearson = _weighted_correlation(truth_gpu, prediction_gpu, weights)
        spearman = _weighted_correlation(truth_ranks, prediction_ranks, weights)
        replicates["pearson"].append(pearson.masked_fill(~correlation_defined, torch.nan))
        replicates["spearman"].append(spearman.masked_fill(~correlation_defined, torch.nan))
        replicates["rmse"].append(((weights * error_gpu.square()).sum(dim=1) / total).sqrt())
        replicates["mae"].append((weights * error_gpu.abs()).sum(dim=1) / total)
    output["bootstrap_completed"] = bootstrap_replicates
    output["bootstrap_nonzero"] = int(nonzero_count)
    output["bootstrap_zero_weight"] = bootstrap_replicates - output["bootstrap_nonzero"]
    probabilities = torch.tensor([0.025, 0.975], dtype=torch.float64, device=device)
    for metric, batches in replicates.items():
        values = torch.cat(batches)
        values = values[torch.isfinite(values)]
        output[f"{metric}_bootstrap_valid"] = len(values)
        if len(values):
            bounds = torch.quantile(values, probabilities).cpu().tolist()
            output[f"{metric}_ci_low"], output[f"{metric}_ci_high"] = bounds
    return output


def compute_calibration_radius(abs_residual, coverage=0.9):
    """Finite-sample split-conformal radius at ceil((n+1)*coverage).

    Append +infinity as the (n+1)th order statistic. When the calibration
    sample is too small for a finite radius, return infinity instead of
    clipping the rank and silently claiming unsupported coverage.
    """
    residual = np.asarray(abs_residual, dtype=np.float64)
    residual = np.sort(residual[np.isfinite(residual)])
    rank = math.ceil((len(residual) + 1) * coverage)
    return float(np.append(residual, np.inf)[rank - 1])


def interval_metrics(y_true, y_pred, radius):
    """Empirical held-out coverage and mean width for symmetric intervals."""
    truth, prediction = np.asarray(y_true, dtype=np.float64), np.asarray(y_pred, dtype=np.float64)
    observed = np.isfinite(truth) & np.isfinite(prediction)
    radii = np.broadcast_to(np.asarray(radius, dtype=np.float64), prediction.shape)[observed]
    truth, prediction = truth[observed], prediction[observed]
    n = len(truth)
    return {
        "interval_n": n,
        "interval_coverage": float(np.mean(np.abs(truth - prediction) <= radii)) if n else None,
        "interval_mean_width": float(np.mean(2 * radii)) if n else None,
    }
