"""Training, prediction, comparisons and measured-pool active learning.

Public interface
----------------
PairData uses a fixed gene index shared by gene_features[G,F], edge_index[2,E]
and pairs[N,2]. Context IDs index context_features[C,D] when supplied. Every
observation has a source ID denoting its assay scale. Missing targets/single
effects are NaN. Joint and GI labels are distinct columns. Supplying a GI
residual does not make it an author-reported interaction score.

fit_model(data, split, model_config, train_config) returns fitted preprocessing
and a trained model. predict(fit, data, indices) returns raw-unit predictions
and training-standardized predictions. No preprocessing is fit on validation,
test or candidate data. CUDA is the default; there is no automatic CPU fallback.
"""

from dataclasses import dataclass, field, replace
import time
import numpy as np
from scipy.stats import rankdata
import torch
from torch import nn

from .models import ModelConfig, build_model, normalized_graph
from .splits import Split
from .metrics import evaluate_predictions, regression_metrics


@dataclass
class PairData:
    gene_features: np.ndarray
    edge_index: np.ndarray
    pairs: np.ndarray
    context_ids: np.ndarray
    joint: np.ndarray
    source_ids: np.ndarray
    single_effects: np.ndarray | None = None
    gi_proxy: np.ndarray | None = None
    context_features: np.ndarray | None = None
    gi_target_name: str = "additive_gi_proxy"
    metadata: dict = field(default_factory=dict)

    def targets(self):
        gi = np.full(len(self.pairs), np.nan) if self.gi_proxy is None else self.gi_proxy
        return np.column_stack((self.joint, gi)).astype(np.float32)

    def singles(self):
        if self.single_effects is None:
            return np.full((len(self.pairs), 2), np.nan, dtype=np.float32)
        return np.asarray(self.single_effects, dtype=np.float32)


@dataclass
class TrainConfig:
    seed: int = 0
    device: str = "cuda"
    max_epochs: int = 80
    patience: int = 10
    batch_size: int = 8192
    learning_rate: float = 0.001
    weight_decay: float = 0.0001
    gi_weight: float = 1.0
    min_delta: float = 0.0001
    amp: bool = True
    gradient_clip: float = 5.0
    log_every: int = 10


def _moments(values):
    valid = np.isfinite(values)
    count = valid.sum(axis=0)
    clean = np.where(valid, values, 0)
    mean = clean.sum(axis=0) / np.maximum(count, 1)
    variance = np.where(valid, (values - mean) ** 2, 0).sum(axis=0) / np.maximum(count, 1)
    std = np.sqrt(variance)
    std = np.where(std > 1e-8, std, 1.0)
    return mean.astype(np.float32), std.astype(np.float32), count


@dataclass
class FeatureNormalizer:
    mean: np.ndarray
    std: np.ndarray
    active: np.ndarray

    @classmethod
    def fit(cls, values, indices):
        selected = np.asarray(values)[indices]
        mean, std, _ = _moments(selected)
        active = np.nanstd(selected, axis=0) > 1e-8
        return cls(mean, std, active)

    def transform(self, values):
        standardized = (np.asarray(values) - self.mean) / self.std
        # Only coordinates supported by training variation have learned weights.
        standardized[:, ~self.active] = 0
        return np.where(np.isfinite(standardized), standardized, 0).astype(np.float32)


@dataclass
class SourceNormalizer:
    sources: np.ndarray
    mean: np.ndarray
    std: np.ndarray
    active: np.ndarray

    @classmethod
    def fit(cls, values, source_ids, indices):
        values, source_ids = np.asarray(values), np.asarray(source_ids)
        selected = values[indices]
        global_mean, global_std, count = _moments(selected)
        sources = np.unique(source_ids[indices])
        means, stds = [global_mean], [global_std]
        for source in sources:
            mean, std, source_count = _moments(selected[source_ids[indices] == source])
            means.append(np.where(source_count > 0, mean, global_mean))
            stds.append(np.where(source_count > 0, std, global_std))
        return cls(sources, np.stack(means), np.stack(stds), count > 0)

    def rows(self, source_ids):
        rows = np.zeros(len(source_ids), dtype=np.int64)
        for i, source in enumerate(self.sources, start=1):
            rows[np.asarray(source_ids) == source] = i
        return rows

    def transform(self, values, source_ids):
        rows = self.rows(source_ids)
        return ((np.asarray(values) - self.mean[rows]) / self.std[rows]).astype(np.float32)

    def inverse(self, values, source_ids):
        rows = self.rows(source_ids)
        return np.asarray(values) * self.std[rows] + self.mean[rows]


@dataclass
class Preprocessing:
    genes: FeatureNormalizer
    contexts: FeatureNormalizer | None
    labels: SourceNormalizer
    singles: SourceNormalizer
    known_genes: np.ndarray
    known_contexts: np.ndarray

    @classmethod
    def fit(cls, data, indices):
        pairs = np.asarray(data.pairs)
        genes = np.unique(pairs[indices])
        contexts = np.unique(np.asarray(data.context_ids)[indices])
        single_sources = np.repeat(data.source_ids, 2)
        single_indices = (np.asarray(indices)[:, None] * 2 + np.arange(2)).reshape(-1)
        return cls(
            FeatureNormalizer.fit(data.gene_features, genes),
            None if data.context_features is None else FeatureNormalizer.fit(data.context_features, contexts),
            SourceNormalizer.fit(data.targets(), data.source_ids, indices),
            SourceNormalizer.fit(data.singles().reshape(-1, 1), single_sources, single_indices),
            genes, contexts,
        )

    def context_rows(self, data):
        if self.contexts is not None:
            return np.asarray(data.context_ids)
        rows = np.zeros(len(data.pairs), dtype=np.int64)
        for i, context in enumerate(self.known_contexts, start=1):
            rows[np.asarray(data.context_ids) == context] = i
        return rows


@dataclass
class PreparedData:
    gene_features: torch.Tensor
    pairs: torch.Tensor
    context_ids: torch.Tensor
    context_features: torch.Tensor | None
    single_effects: torch.Tensor
    single_missing: torch.Tensor
    targets: torch.Tensor
    graph: torch.Tensor | None


def _prepare(data, preprocessing, device, use_graph, included_genes=None):
    def floating(values):
        return torch.as_tensor(values, dtype=torch.float32, device=device)

    singles = data.singles()
    standardized_singles = preprocessing.singles.transform(
        singles.reshape(-1, 1), np.repeat(data.source_ids, 2)
    ).reshape(-1, 2)
    missing = ~np.isfinite(singles)
    standardized_singles[missing] = 0
    context_features = (None if preprocessing.contexts is None
                        else floating(preprocessing.contexts.transform(data.context_features)))
    graph = (normalized_graph(data.edge_index, len(data.gene_features), device, included_genes)
             if use_graph else None)
    targets = preprocessing.labels.transform(data.targets(), data.source_ids)
    targets[:, ~preprocessing.labels.active] = np.nan
    return PreparedData(
        floating(preprocessing.genes.transform(data.gene_features)),
        torch.as_tensor(data.pairs, dtype=torch.long, device=device),
        torch.as_tensor(preprocessing.context_rows(data), dtype=torch.long, device=device),
        context_features, floating(standardized_singles), floating(missing),
        floating(targets), graph,
    )


@dataclass
class FitResult:
    model: nn.Module
    preprocessing: Preprocessing
    model_config: ModelConfig
    train_config: TrainConfig
    split: Split
    history: list
    best_epoch: int
    elapsed_seconds: float
    gi_target_name: str


def _forward(model, prepared, indices, gene_embeddings=None):
    if gene_embeddings is None:
        gene_embeddings = model.encode_genes(prepared.gene_features, prepared.graph)
    return model(gene_embeddings, prepared.pairs[indices], prepared.context_ids[indices],
                 prepared.context_features, prepared.single_effects[indices],
                 prepared.single_missing[indices])


def _task_loss(prediction, targets, gi_weight):
    mask = torch.isfinite(targets)
    error = torch.where(mask, prediction.float() - targets, 0).square()
    counts = mask.sum(dim=0)
    means = error.sum(dim=0) / counts.clamp_min(1)
    weights = prediction.new_tensor([1.0, gi_weight]) * (counts > 0)
    return (means * weights).sum() / weights.sum().clamp_min(1)


@torch.no_grad()
def _validation_loss(model, prepared, indices, config):
    model.eval()
    errors = torch.zeros(2, device=config.device)
    counts = torch.zeros(2, device=config.device)
    with torch.autocast(device_type=torch.device(config.device).type,
                        enabled=config.amp and torch.device(config.device).type == "cuda"):
        embeddings = model.encode_genes(prepared.gene_features, prepared.graph)
        for batch in indices.split(config.batch_size):
            prediction = _forward(model, prepared, batch, embeddings).float()
            targets = prepared.targets[batch]
            valid = torch.isfinite(targets)
            errors += torch.where(valid, prediction - targets, 0).square().sum(dim=0)
            counts += valid.sum(dim=0)
    weights = errors.new_tensor([1.0, config.gi_weight]) * (counts > 0)
    return float(((errors / counts.clamp_min(1)) * weights).sum() / weights.sum().clamp_min(1))


def fit_model(data, split, model_config=None, train_config=None, progress=None):
    """Fit on train, select an epoch on validation, and never inspect test labels.

    GNN training uses the graph induced by training genes. Inference applies
    the learned feature/aggregation functions to the complete external graph.
    A validation partition with at least one measured target is required.
    """
    model_config = ModelConfig() if model_config is None else model_config
    config = TrainConfig() if train_config is None else train_config
    start = time.monotonic()
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    if torch.device(config.device).type == "cuda":
        torch.cuda.manual_seed_all(config.seed)
        torch.backends.cuda.matmul.allow_tf32 = True
    preprocessing = Preprocessing.fit(data, split.train)
    context_dim = 0 if data.context_features is None else data.context_features.shape[1]
    model = build_model(
        data.gene_features.shape[1], len(data.gene_features), len(preprocessing.known_contexts),
        context_dim, torch.as_tensor(preprocessing.known_genes), model_config,
    ).to(config.device)
    prepared = _prepare(data, preprocessing, config.device, model.use_graph,
                        torch.as_tensor(preprocessing.known_genes, device=config.device))
    validation = replace(prepared, graph=normalized_graph(
        data.edge_index, len(data.gene_features), config.device
    )) if model.use_graph else prepared
    train_indices = torch.as_tensor(split.train, dtype=torch.long, device=config.device)
    val_indices = torch.as_tensor(split.val, dtype=torch.long, device=config.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate,
                                  weight_decay=config.weight_decay)
    amp = config.amp and torch.device(config.device).type == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    history, best_loss, best_epoch, stale = [], float("inf"), 0, 0
    best_state = None
    for epoch in range(1, config.max_epochs + 1):
        model.train()
        shuffled = train_indices[torch.randperm(len(train_indices), device=config.device)]
        train_loss = torch.zeros((), device=config.device)
        n_batches = 0
        for batch in shuffled.split(config.batch_size):
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type=torch.device(config.device).type, enabled=amp):
                prediction = _forward(model, prepared, batch)
                loss = _task_loss(prediction, prepared.targets[batch], config.gi_weight)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            nn.utils.clip_grad_norm_(model.parameters(), config.gradient_clip)
            scaler.step(optimizer)
            scaler.update()
            train_loss += loss.detach()
            n_batches += 1
        val_loss = _validation_loss(model, validation, val_indices, config)
        entry = {"epoch": epoch, "train_loss": float(train_loss / n_batches), "val_loss": val_loss}
        history.append(entry)
        if val_loss < best_loss - config.min_delta:
            best_loss, best_epoch, stale = val_loss, epoch, 0
            best_state = {name: tensor.detach().clone() for name, tensor in model.state_dict().items()}
        else:
            stale += 1
        if progress is not None and (epoch == 1 or epoch % config.log_every == 0):
            progress(entry)
        if stale >= config.patience:
            break
    model.load_state_dict(best_state)
    model.eval()
    return FitResult(model, preprocessing, model_config, config, split, history,
                     best_epoch, time.monotonic() - start, data.gi_target_name)


@torch.no_grad()
def predict(result, data, indices=None, mc_passes=1, seed=None):
    """Predict raw and standardized targets; MC dropout supplies heuristic spread.

    GI predictions are NaN if the model received no GI supervision. Unseen
    source scales use training-global mean/std, never target-dataset labels.
    MC spread is not a calibrated confidence interval.
    """
    model, config = result.model, result.train_config
    indices = np.arange(len(data.pairs)) if indices is None else np.asarray(indices)
    if seed is not None:
        torch.manual_seed(seed)
    prepared = _prepare(data, result.preprocessing, config.device, model.use_graph)
    rows = torch.as_tensor(indices, dtype=torch.long, device=config.device)
    stochastic = mc_passes > 1
    model.train(stochastic)
    # Keep contextual identity fixed during MC inference; only neural dropout varies.
    if model.context_encoder is not None:
        model.context_encoder.eval()
    predictions = []
    amp = config.amp and torch.device(config.device).type == "cuda"
    for _ in range(mc_passes):
        batches = []
        with torch.autocast(device_type=torch.device(config.device).type, enabled=amp):
            embeddings = model.encode_genes(prepared.gene_features, prepared.graph)
            for batch in rows.split(config.batch_size):
                batches.append(_forward(model, prepared, batch, embeddings).float())
        predictions.append(torch.cat(batches).cpu().numpy())
    model.eval()
    samples = np.stack(predictions)
    standardized, standardized_std = samples.mean(axis=0), samples.std(axis=0)
    active = result.preprocessing.labels.active
    standardized[:, ~active] = np.nan
    standardized_std[:, ~active] = np.nan
    sources = np.asarray(data.source_ids)[indices]
    raw = result.preprocessing.labels.inverse(standardized, sources)
    raw_std = standardized_std * result.preprocessing.labels.std[result.preprocessing.labels.rows(sources)]
    output = {
        "joint": raw[:, 0], "gi": raw[:, 1],
        "joint_std": raw_std[:, 0], "gi_std": raw_std[:, 1],
        "joint_standardized": standardized[:, 0], "gi_standardized": standardized[:, 1],
        "joint_std_standardized": standardized_std[:, 0], "gi_std_standardized": standardized_std[:, 1],
        "gi_target_name": result.gi_target_name,
        "source_scale_seen": np.isin(sources, result.preprocessing.labels.sources),
    }
    output[result.gi_target_name] = output["gi"]
    if result.gi_target_name == "additive_gi_proxy":
        output["gi_from_joint"] = np.where(output["source_scale_seen"],
            output["joint"] - data.singles()[indices].sum(axis=1), np.nan)
    return output


def evaluate_fit(result, data, indices=None, top_k=100, joint_direction=-1, gi_direction=-1):
    indices = result.split.test if indices is None else np.asarray(indices)
    prediction = predict(result, data, indices)
    truth = data.targets()[indices]
    standardized_truth = result.preprocessing.labels.transform(truth, np.asarray(data.source_ids)[indices])
    arguments = dict(top_k=top_k, joint_direction=joint_direction,
                     gi_direction=gi_direction, gi_target_name=result.gi_target_name,
                     source_scale_seen=prediction["source_scale_seen"])
    evaluation = {
        "raw": evaluate_predictions(truth[:, 0], prediction["joint"], truth[:, 1], prediction["gi"],
                                    np.asarray(data.source_ids)[indices], **arguments),
        "training_standardized": evaluate_predictions(
            standardized_truth[:, 0], prediction["joint_standardized"],
            standardized_truth[:, 1], prediction["gi_standardized"], **arguments),
        "unseen_source_count": int((~prediction["source_scale_seen"]).sum()),
    }
    if result.gi_target_name == "additive_gi_proxy":
        sources = np.asarray(data.source_ids)[indices]
        evaluation["gi_from_joint"] = {
            str(source): regression_metrics(
                truth[sources == source, 1], prediction["gi_from_joint"][sources == source],
                top_k=top_k, direction=gi_direction,
            ) for source in np.unique(sources)
        }
    return evaluation


def evaluate_additive_reference(data, indices, top_k=100, joint_direction=-1, gi_direction=-1):
    """Unfitted observed-single sum baseline, in the source's additive units.

    This reference is valid only for an observed-minus-additive GI definition.
    A separately fitted additive feature model remains available as 'additive'.
    """
    indices = np.asarray(indices)
    single_sum = data.singles()[indices].sum(axis=1)
    gi_prediction = np.where(np.isfinite(single_sum), 0.0, np.nan)
    truth = data.targets()[indices]
    return evaluate_predictions(
        truth[:, 0], single_sum, truth[:, 1], gi_prediction,
        np.asarray(data.source_ids)[indices], top_k=top_k,
        joint_direction=joint_direction, gi_direction=gi_direction,
        gi_target_name=data.gi_target_name,
    )


def benchmark(data, splits, variants=("additive", "mlp", "naiad_reimplementation", "gnn", "graphless", "contextless"),
              seeds=(11, 29, 47), model_config=None, train_config=None, progress=None):
    """Compare every variant on identical supplied partitions and initialization seeds."""
    base_model = ModelConfig() if model_config is None else model_config
    base_train = TrainConfig() if train_config is None else train_config
    records = []
    for split in splits:
        for seed in seeds:
            for variant in variants:
                result = fit_model(data, split, replace(base_model, variant=variant),
                                   replace(base_train, seed=seed), progress=progress)
                records.append({
                    "split": split.name, "split_seed": split.seed, "seed": seed, "model": variant,
                    "n_train": len(split.train), "n_val": len(split.val), "n_test": len(split.test),
                    "n_discarded": 0 if split.discarded is None else len(split.discarded),
                    "best_epoch": result.best_epoch, "elapsed_seconds": result.elapsed_seconds,
                    "inductive_gene_features": variant != "naiad_reimplementation",
                    "metrics": evaluate_fit(result, data),
                })
                if progress is not None:
                    progress(records[-1])
                del result
    return records


def active_learning(data, split, initial_size, acquisition_size, rounds,
                    strategies=("random", "effect", "uncertainty"), seed=11,
                    model_config=None, train_config=None, mc_passes=12,
                    task="joint", effect_direction=-1, progress=None):
    """Retrospective acquisition from measured training pairs only.

    Strategies share the same initial measured pool and fixed val/test sets.
    Neither validation nor test observations can enter the acquisition pool.
    Each round refits preprocessing using only the currently observed labels.
    Uncertainty is MC-dropout spread, not calibrated assay uncertainty.
    The optional balanced strategy uses GI: 0.7 effect percentile plus 0.3
    uncertainty percentile, with a 0.05 penalty per repeated gene in the same
    acquisition round. These fixed weights are a heuristic, not tuned values.
    Use task='gi' to compare all strategies on the same measured GI pool.
    """
    base_model = ModelConfig() if model_config is None else model_config
    base_train = TrainConfig() if train_config is None else train_config
    target_column = 0 if task == "joint" else 1
    pool = np.asarray(split.train)
    pool = pool[np.isfinite(data.targets()[pool, target_column])]
    initial = np.random.default_rng(seed).permutation(pool)[:initial_size]
    records = []
    for strategy in strategies:
        observed = initial.copy()
        rng = np.random.default_rng(seed)
        for step in range(rounds + 1):
            current = Split(observed.copy(), split.val, split.test, f"active_{strategy}", seed)
            fitted = fit_model(data, current, base_model, replace(base_train, seed=seed + step))
            record = {
                "strategy": strategy, "round": step, "seed": seed,
                "n_observed": len(observed), "observed_indices": observed.tolist(),
                "metrics": evaluate_fit(fitted, data), "elapsed_seconds": fitted.elapsed_seconds,
                "acquisition_task": "gi" if strategy == "balanced" else task,
                "acquired_indices": [],
            }
            if strategy == "balanced":
                record["acquisition_rule"] = {
                    "type": "fixed_heuristic", "effect_percentile_weight": 0.7,
                    "uncertainty_percentile_weight": 0.3, "gene_repeat_penalty": 0.05,
                }
            available = pool[~np.isin(pool, observed)]
            if step < rounds and len(available):
                n_acquire = min(acquisition_size, len(available))
                if strategy == "random":
                    selected = rng.permutation(available)[:n_acquire]
                else:
                    prediction = predict(fitted, data, available,
                                         mc_passes=mc_passes if strategy in {"uncertainty", "balanced"} else 1,
                                         seed=seed + step)
                    if strategy == "balanced":
                        effect_percentile = rankdata(
                            effect_direction * prediction["gi_standardized"], method="average"
                        ) / len(available)
                        uncertainty_percentile = rankdata(
                            prediction["gi_std_standardized"], method="average"
                        ) / len(available)
                        utility = 0.7 * effect_percentile + 0.3 * uncertainty_percentile
                        pair_genes = np.asarray(data.pairs)[available]
                        gene_counts = np.zeros(len(data.gene_features), dtype=np.int64)
                        selected = np.empty(n_acquire, dtype=available.dtype)
                        for i in range(n_acquire):
                            penalized = utility - 0.05 * gene_counts[pair_genes].sum(axis=1)
                            position = np.argmax(penalized)
                            selected[i] = available[position]
                            utility[position] = -np.inf
                            np.add.at(gene_counts, pair_genes[position], 1)
                    else:
                        score = (prediction[f"{task}_std_standardized"] if strategy == "uncertainty"
                                 else effect_direction * prediction[f"{task}_standardized"])
                        selected = available[np.argsort(-score, kind="stable")[:n_acquire]]
                record["acquired_indices"] = selected.tolist()
                observed = np.concatenate((observed, selected))
            records.append(record)
            if progress is not None:
                progress(record)
            del fitted
            if not len(available):
                break
    return records


def save_fit(result, path):
    """Save a local fitted model, including its train-only preprocessing."""
    torch.save(result, path)


def load_fit(path, device="cuda"):
    """Load a checkpoint created by save_fit; retain learned preprocessing."""
    result = torch.load(path, map_location=device, weights_only=False)
    result.model.to(device).eval()
    result.train_config = replace(result.train_config, device=device)
    return result
