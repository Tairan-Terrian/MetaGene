"""Fitness and GI are evaluated as separate regression/ranking tasks."""

import numpy as np
from scipy.stats import rankdata


def regression_metrics(y_true, y_pred, top_k=100, direction=-1):
    """Top-k overlap uses the requested effect direction; -1 selects lower values."""
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    valid = np.isfinite(y_true) & np.isfinite(y_pred)
    truth, prediction = y_true[valid], y_pred[valid]
    n = len(truth)
    if n == 0:
        return {"n": 0, "rmse": None, "mae": None, "spearman": None,
                "top_k": 0, "top_k_overlap": None}
    k = min(top_k, max(1, round(0.10 * n)))
    rho = None
    if n > 1 and np.ptp(truth) > 0 and np.ptp(prediction) > 0:
        rho = float(np.corrcoef(rankdata(truth), rankdata(prediction))[0, 1])
    true_top = np.argsort(-direction * truth, kind="stable")[:k]
    pred_top = np.argsort(-direction * prediction, kind="stable")[:k]
    return {
        "n": n, "rmse": float(np.sqrt(np.mean((prediction - truth) ** 2))),
        "mae": float(np.mean(np.abs(prediction - truth))), "spearman": rho,
        "top_k": k, "top_k_overlap": (float(len(np.intersect1d(true_top, pred_top)) / k)
                                        if np.ptp(prediction) > 0 and np.ptp(truth) > 0 else None),
    }


def evaluate_predictions(joint_true, joint_pred, gi_true=None, gi_pred=None,
                         source_ids=None, top_k=100, joint_direction=-1,
                         gi_direction=-1, gi_target_name="gi_proxy", source_scale_seen=None):
    """Return task-specific metrics, with source-specific metrics for mixed assays.

    Raw RMSE pooled across assays with different units is intentionally omitted.
    Call this function on training-standardized values for a pooled comparison.
    """
    joint_true, joint_pred = np.asarray(joint_true), np.asarray(joint_pred)
    if gi_true is None:
        gi_true = np.full(len(joint_true), np.nan)
        gi_pred = np.full(len(joint_true), np.nan)
    gi_true, gi_pred = np.asarray(gi_true), np.asarray(gi_pred)

    def score(mask):
        result = {
            "joint": regression_metrics(joint_true[mask], joint_pred[mask], top_k, joint_direction),
            gi_target_name: regression_metrics(gi_true[mask], gi_pred[mask], top_k, gi_direction),
        }
        if source_scale_seen is not None and not np.asarray(source_scale_seen)[mask].all():
            for metrics in result.values():
                metrics["rmse"] = None
                metrics["mae"] = None
        return result

    if source_ids is None or len(np.unique(source_ids)) == 1:
        return score(np.ones(len(joint_true), dtype=bool))
    source_ids = np.asarray(source_ids)
    return {"by_source": {str(source): score(source_ids == source)
                          for source in np.unique(source_ids)}}
