"""Residual GI models with optional, independently gated annotation graphs.

The prediction is a single standardized additive GI residual. Joint effects
are reconstructed in source units by residual_training. Every graph uses the same
complete external annotation topology during training and prediction; this
is annotation-transductive, with held-out double-perturbation labels excluded.
"""

from dataclasses import dataclass
import torch
from torch import nn
from torch.nn import functional as F

from .models import ContextEncoder, invariant_pair, invariant_singles


@dataclass
class ResidualModelConfig:
    hidden_dim: int = 96
    context_dim: int = 24
    dropout: float = 0.15
    graph_names: tuple[str, ...] = ()
    gate_init: float = 0.1
    use_single_effects: bool = False


class GatedGeneEncoder(nn.Module):
    """A shared feature branch plus one optional hop per external graph."""

    def __init__(self, feature_dim, config):
        super().__init__()
        h = config.hidden_dim
        self.graph_names = config.graph_names
        self.dropout = config.dropout
        self.input = nn.Sequential(nn.Linear(feature_dim, h), nn.LayerNorm(h), nn.GELU())
        self.self_branch = nn.Sequential(
            nn.Linear(h, h), nn.GELU(), nn.Dropout(config.dropout), nn.Linear(h, h)
        )
        self.graph_projections = nn.ModuleDict({
            name: nn.Linear(h, h, bias=False) for name in config.graph_names
        })
        self.gate_logits = nn.ParameterDict({
            name: nn.Parameter(torch.logit(torch.tensor(config.gate_init)))
            for name in config.graph_names
        })
        self.norm = nn.LayerNorm(h)

    def forward(self, features, graphs):
        h0 = self.input(features)
        h = h0 + self.self_branch(h0)
        for name in self.graph_names:
            with torch.autocast(device_type=h0.device.type, enabled=False):
                neighbors = torch.sparse.mm(graphs[name], h0.float())
            message = F.gelu(self.graph_projections[name](neighbors))
            message = F.dropout(message, self.dropout, self.training)
            h = h + self.gate_logits[name].sigmoid() * message
        return self.norm(h)

    def gate_values(self):
        return {name: float(value.detach().sigmoid())
                for name, value in self.gate_logits.items()}


class ResidualGIModel(nn.Module):
    """Symmetric gene-pair model whose only learned outcome is additive GI."""

    def __init__(self, feature_dim, n_contexts, context_features_dim, config):
        super().__init__()
        self.config = config
        self.gene_encoder = GatedGeneEncoder(feature_dim, config)
        self.context_encoder = ContextEncoder(
            n_contexts, context_features_dim, config.context_dim, config.dropout
        )
        pair_dim = 3 * config.hidden_dim + config.context_dim
        if config.use_single_effects:
            pair_dim += 4
        self.gi_head = nn.Sequential(
            nn.Linear(pair_dim, 2 * config.hidden_dim), nn.GELU(), nn.Dropout(config.dropout),
            nn.Linear(2 * config.hidden_dim, config.hidden_dim), nn.GELU(),
            nn.Dropout(config.dropout), nn.Linear(config.hidden_dim, 1),
        )

    def encode_genes(self, gene_features, graphs):
        return self.gene_encoder(gene_features, graphs)

    def forward(self, gene_embeddings, pairs, context_ids, context_features,
                single_effects, single_missing):
        a, b = gene_embeddings[pairs[:, 0]], gene_embeddings[pairs[:, 1]]
        parts = [invariant_pair(a, b)]
        if self.config.use_single_effects:
            parts.append(invariant_singles(single_effects, single_missing))
        parts.append(self.context_encoder(context_ids, context_features))
        return self.gi_head(torch.cat(parts, dim=-1)).squeeze(-1)
