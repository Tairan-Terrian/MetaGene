"""GI-only training with source-balanced loss and structurally consistent joint.

Reuse the two-head PairData and Split containers. graph_edges maps a graph name to
an external edge_index[2,E] using the gene_features row index. Only configured
graphs are used. All feature and label transforms are fitted on train rows;
all graph topologies remain fixed between training, validation and prediction.

Known source: predicted_joint = measured_single_a + measured_single_b + GI.
Unseen source: gi_score is a ranking score; raw GI/joint are not calibrated
and are returned as NaN. Missing singleton values likewise leave joint NaN.
"""

from dataclasses import dataclass, field, replace
import time
import numpy as np
import torch
from torch import nn

from .models import normalized_graph
from .metrics import regression_metrics
from .training import PairData, Preprocessing
from .splits import Split
from .residual_models import ResidualModelConfig, ResidualGIModel


@dataclass
class ResidualTrainConfig:
    seed: int = 17
    device: str = "cuda"
    max_epochs: int = 100
    patience: int = 15
    batch_size: int = 4096
    learning_rate: float = 0.001
    weight_decay: float = 0.0001
    min_delta: float = 0.0001
    amp: bool = True
    gradient_clip: float = 5.0
    log_every: int = 10
    source_weights: dict[int, float] = field(default_factory=lambda: {0: 0.25, 1: 0.25, 2: 0.5})
    validation_source: int | None = 2


@dataclass
class ResidualPreparedData:
    gene_features: torch.Tensor
    graphs: dict[str, torch.Tensor]
    pairs: torch.Tensor
    context_ids: torch.Tensor
    context_features: torch.Tensor | None
    single_effects: torch.Tensor
    single_missing: torch.Tensor
    raw_single_sum: torch.Tensor
    gi_target: torch.Tensor
    source_rows: torch.Tensor


@dataclass
class ResidualFitResult:
    model: ResidualGIModel
    preprocessing: Preprocessing
    model_config: ResidualModelConfig
    train_config: ResidualTrainConfig
    split: Split
    history: list
    best_epoch: int
    best_validation_rmse: float
    validation_scale: str
    elapsed_seconds: float
    source_statistics: dict
    learned_graph_gates: dict
    gi_target_name: str = "additive_gi_proxy"
    graph_scope: str = "fixed_external_annotation_transductive"


def _prepare(data, preprocessing, device, graph_names, graph_edges):
    def floating(values):
        return torch.as_tensor(values, dtype=torch.float32, device=device)

    singles = data.singles()
    missing = ~np.isfinite(singles)
    normalized_singles = preprocessing.singles.transform(
        singles.reshape(-1, 1), np.repeat(data.source_ids, 2)
    ).reshape(-1, 2)
    normalized_singles[missing] = 0
    contexts = (None if preprocessing.contexts is None
                else floating(preprocessing.contexts.transform(data.context_features)))
    graphs = {name: normalized_graph(graph_edges[name], len(data.gene_features), device)
              for name in graph_names}
    return ResidualPreparedData(
        floating(preprocessing.genes.transform(data.gene_features)), graphs,
        torch.as_tensor(data.pairs, dtype=torch.long, device=device),
        torch.as_tensor(preprocessing.context_rows(data), dtype=torch.long, device=device),
        contexts, floating(normalized_singles), floating(missing), floating(singles.sum(axis=1)),
        floating(preprocessing.labels.transform(data.targets(), data.source_ids)[:, 1]),
        torch.as_tensor(preprocessing.labels.rows(data.source_ids), dtype=torch.long, device=device),
    )


def _source_weights(preprocessing, config):
    weights = [0.0] + [config.source_weights[int(source)] for source in preprocessing.labels.sources]
    return torch.tensor(weights, dtype=torch.float32, device=config.device)


def _source_sums(prediction, target, source_rows, n_sources):
    measured = torch.isfinite(target)
    squared_error = torch.where(measured, prediction.float() - target, 0).square()
    errors = torch.zeros(n_sources, device=prediction.device).scatter_add_(0, source_rows, squared_error)
    counts = torch.zeros_like(errors).scatter_add_(0, source_rows, measured.float())
    return errors, counts


def _weighted_mean(errors, counts, source_weights):
    weights = source_weights * (counts > 0)
    return ((errors / counts.clamp_min(1)) * weights).sum() / weights.sum()


def _forward(model, prepared, indices, embeddings=None):
    if embeddings is None:
        embeddings = model.encode_genes(prepared.gene_features, prepared.graphs)
    return model(embeddings, prepared.pairs[indices], prepared.context_ids[indices],
                 prepared.context_features, prepared.single_effects[indices], prepared.single_missing[indices])


@torch.no_grad()
def _validation(model, prepared, indices, preprocessing, source_weights, config):
    model.eval()
    errors = torch.zeros_like(source_weights)
    counts = torch.zeros_like(source_weights)
    predictions = []
    amp = config.amp and torch.device(config.device).type == "cuda"
    with torch.autocast(device_type=torch.device(config.device).type, enabled=amp):
        embeddings = model.encode_genes(prepared.gene_features, prepared.graphs)
        for batch in indices.split(config.batch_size):
            prediction = _forward(model, prepared, batch, embeddings).float()
            batch_errors, batch_counts = _source_sums(
                prediction, prepared.gi_target[batch], prepared.source_rows[batch], len(source_weights)
            )
            errors += batch_errors
            counts += batch_counts
            predictions.append(prediction)
    loss = _weighted_mean(errors, counts, source_weights)
    prediction = torch.cat(predictions).cpu().numpy()
    truth = prepared.gi_target[indices].cpu().numpy()
    rho = regression_metrics(truth, prediction)["spearman"]
    if config.validation_source is None:
        rmse = float(loss.sqrt())
    else:
        row = preprocessing.labels.rows(np.array([config.validation_source]))[0]
        rmse = float(loss.sqrt()) * float(preprocessing.labels.std[row, 1])
    return float(loss), rmse, rho


def _source_statistics(data, split, preprocessing):
    labels = preprocessing.labels
    statistics = {
        "global_training_reference": {
            "joint_mean": float(labels.mean[0, 0]), "joint_std": float(labels.std[0, 0]),
            "gi_mean": float(labels.mean[0, 1]), "gi_std": float(labels.std[0, 1]),
        },
        "by_source": {},
        "active_gene_features": int(preprocessing.genes.active.sum()),
        "total_gene_features": len(preprocessing.genes.active),
    }
    for row, source in enumerate(labels.sources, start=1):
        selected = split.train[np.asarray(data.source_ids)[split.train] == source]
        single_row = preprocessing.singles.rows(np.array([source]))[0]
        statistics["by_source"][str(int(source))] = {
            "n_train": len(selected), "n_gi": int(np.isfinite(data.targets()[selected, 1]).sum()),
            "joint_mean": float(labels.mean[row, 0]), "joint_std": float(labels.std[row, 0]),
            "gi_mean": float(labels.mean[row, 1]), "gi_std": float(labels.std[row, 1]),
            "single_mean": float(preprocessing.singles.mean[single_row, 0]),
            "single_std": float(preprocessing.singles.std[single_row, 0]),
        }
    return statistics


def fit_residual(data: PairData, split: Split, model_config=None, train_config=None,
           graph_edges=None, progress=None):
    """Train standardized GI; select on fixed-source validation GI RMSE.

    Source means are averaged with source_weights, renormalized over sources
    present in each batch. validation_source=2 selects only HT29 validation
    rows. For a protocol without HT29 validation, set validation_source=None;
    selection then uses the weighted, source-standardized validation RMSE.
    The supplied GI labels are observed-minus-additive residuals, not author
    interaction scores on a different scale.
    """
    model_config = ResidualModelConfig() if model_config is None else model_config
    config = ResidualTrainConfig() if train_config is None else train_config
    start = time.monotonic()
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    if torch.device(config.device).type == "cuda":
        torch.cuda.manual_seed_all(config.seed)
        torch.backends.cuda.matmul.allow_tf32 = True
    preprocessing = Preprocessing.fit(data, split.train)
    context_dim = 0 if data.context_features is None else data.context_features.shape[1]
    model = ResidualGIModel(
        data.gene_features.shape[1], len(preprocessing.known_contexts), context_dim, model_config
    ).to(config.device)
    prepared = _prepare(data, preprocessing, config.device, model_config.graph_names, graph_edges)
    source_weights = _source_weights(preprocessing, config)
    train_indices = torch.as_tensor(split.train, dtype=torch.long, device=config.device)
    validation_rows = (split.val if config.validation_source is None else
                       split.val[np.asarray(data.source_ids)[split.val] == config.validation_source])
    val_indices = torch.as_tensor(validation_rows, dtype=torch.long, device=config.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
    amp = config.amp and torch.device(config.device).type == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    history, best_loss, best_epoch, stale = [], float("inf"), 0, 0
    best_state, best_validation_rmse = None, float("inf")
    for epoch in range(1, config.max_epochs + 1):
        model.train()
        shuffled = train_indices[torch.randperm(len(train_indices), device=config.device)]
        train_loss = torch.zeros((), device=config.device)
        n_batches = 0
        for batch in shuffled.split(config.batch_size):
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type=torch.device(config.device).type, enabled=amp):
                prediction = _forward(model, prepared, batch)
                errors, counts = _source_sums(
                    prediction, prepared.gi_target[batch], prepared.source_rows[batch], len(source_weights)
                )
                loss = _weighted_mean(errors, counts, source_weights)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            nn.utils.clip_grad_norm_(model.parameters(), config.gradient_clip)
            scaler.step(optimizer)
            scaler.update()
            train_loss += loss.detach()
            n_batches += 1
        val_loss, val_rmse, val_rho = _validation(
            model, prepared, val_indices, preprocessing, source_weights, config
        )
        entry = {"epoch": epoch, "train_gi_loss": float(train_loss / n_batches),
                 "val_gi_loss": val_loss, "val_gi_rmse": val_rmse, "val_gi_spearman": val_rho}
        history.append(entry)
        stale = 0 if val_loss < best_loss - config.min_delta else stale + 1
        if val_loss < best_loss:
            best_loss, best_epoch, best_validation_rmse = val_loss, epoch, val_rmse
            best_state = {name: value.detach().clone() for name, value in model.state_dict().items()}
        if progress is not None and (epoch == 1 or epoch % config.log_every == 0):
            progress(entry)
        if stale >= config.patience:
            break
    model.load_state_dict(best_state)
    model.eval()
    return ResidualFitResult(
        model, preprocessing, model_config, config, split, history, best_epoch,
        best_validation_rmse, "source_standardized" if config.validation_source is None else "source_raw",
        time.monotonic() - start, _source_statistics(data, split, preprocessing),
        model.gene_encoder.gate_values(), data.gi_target_name,
    )


@torch.no_grad()
def predict_residual(result, data, indices=None, graph_edges=None, mc_passes=1,
               seed=None, return_samples=False):
    """GPU prediction and MC aggregation; exported arrays are ready for ranking.

    gi_score is in the train-fitted standardized residual space. For seen
    sources, gi and joint are in training-fitted source units and share exactly
    the same residual error. MC spread is a heuristic, not a confidence interval.
    """
    model, config = result.model, result.train_config
    indices = np.arange(len(data.pairs)) if indices is None else np.asarray(indices)
    if seed is not None:
        torch.manual_seed(seed)
    prepared = _prepare(data, result.preprocessing, config.device,
                           result.model_config.graph_names, graph_edges)
    rows = torch.as_tensor(indices, dtype=torch.long, device=config.device)
    model.train(mc_passes > 1)
    model.context_encoder.eval()
    amp = config.amp and torch.device(config.device).type == "cuda"
    samples = []
    for _ in range(mc_passes):
        batches = []
        with torch.autocast(device_type=torch.device(config.device).type, enabled=amp):
            embeddings = model.encode_genes(prepared.gene_features, prepared.graphs)
            for batch in rows.split(config.batch_size):
                batches.append(_forward(model, prepared, batch, embeddings).float())
        samples.append(torch.cat(batches))
    model.eval()
    samples = torch.stack(samples)
    score, score_std = samples.mean(dim=0), samples.std(dim=0, correction=0)
    label_rows = prepared.source_rows[rows]
    seen = label_rows > 0
    means = torch.as_tensor(result.preprocessing.labels.mean, device=config.device)[label_rows, 1]
    scales = torch.as_tensor(result.preprocessing.labels.std, device=config.device)[label_rows, 1]
    gi = (score * scales + means).masked_fill(~seen, torch.nan)
    gi_std = (score_std * scales).masked_fill(~seen, torch.nan)
    joint = gi + prepared.raw_single_sum[rows]
    singles_measured = torch.isfinite(prepared.raw_single_sum[rows])
    joint_std = gi_std.masked_fill(~singles_measured, torch.nan)
    output = {
        "joint": joint.cpu().numpy(), "gi": gi.cpu().numpy(),
        "gi_from_joint": gi.masked_fill(~singles_measured, torch.nan).cpu().numpy(),
        "joint_std": joint_std.cpu().numpy(),
        "gi_std": gi_std.cpu().numpy(), "gi_score": score.cpu().numpy(),
        "gi_score_std": score_std.cpu().numpy(), "gi_standardized": score.cpu().numpy(),
        "gi_std_standardized": score_std.cpu().numpy(),
        "source_scale_seen": seen.cpu().numpy(), "gi_target_name": result.gi_target_name,
    }
    output[result.gi_target_name] = output["gi"]
    if return_samples:
        output["gi_score_samples"] = samples.cpu().numpy()
        output["gi_samples"] = (samples * scales + means).masked_fill(~seen, torch.nan).cpu().numpy()
    return output


def save_residual(result, path):
    """Persist the fitted model and its training-only source statistics."""
    torch.save(result, path)


def load_residual(path, device="cuda"):
    """Read a checkpoint produced by save_residual."""
    result = torch.load(path, map_location=device, weights_only=False)
    result.model.to(device).eval()
    result.train_config = replace(result.train_config, device=device)
    return result
