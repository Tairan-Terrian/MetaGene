"""CUDA-friendly invariant pair regressors with shared gene/context encoders.

Gene features and graph edges are external biological inputs. No gene-identity
embedding is used by the additive, MLP or graph models. The second output is
an explicitly supplied GI target, never an inferred author-reported score.
"""

from dataclasses import dataclass
import torch
from torch import nn
from torch.nn import functional as F


@dataclass
class ModelConfig:
    variant: str = "gnn"
    hidden_dim: int = 96
    graph_layers: int = 2
    context_dim: int = 24
    dropout: float = 0.15
    use_single_effects: bool = True
    naiad_embedding_dim: int = 16


def normalized_graph(edge_index, n_genes, device, included_genes=None):
    """Undirected GCN normalization, including one self edge per retained gene."""
    edges = torch.as_tensor(edge_index, dtype=torch.long, device=device)
    genes = torch.arange(n_genes, device=device)
    if included_genes is not None:
        keep = torch.zeros(n_genes, dtype=torch.bool, device=device)
        keep[included_genes] = True
        edges = edges[:, keep[edges].all(dim=0)]
        genes = genes[keep]
    edges = torch.cat((edges, edges.flip(0), torch.stack((genes, genes))), dim=1)
    graph = torch.sparse_coo_tensor(
        edges, torch.ones(edges.shape[1], device=device), (n_genes, n_genes)
    ).coalesce()
    edges = graph.indices()
    # Repeated input edges do not create stronger biological relationships.
    weights = torch.ones(edges.shape[1], device=device)
    degree = torch.zeros(n_genes, device=device).scatter_add_(0, edges[0], weights)
    inv_degree = degree.clamp_min(1).rsqrt()
    weights = inv_degree[edges[0]] * inv_degree[edges[1]]
    return torch.sparse_coo_tensor(edges, weights, (n_genes, n_genes)).coalesce()


class ContextEncoder(nn.Module):
    def __init__(self, n_contexts, context_features_dim, output_dim, dropout, linear_only=False):
        super().__init__()
        self.use_features = context_features_dim > 0
        self.dropout = dropout
        if self.use_features:
            layers = [nn.Linear(context_features_dim, output_dim)]
            if not linear_only:
                layers.append(nn.GELU())
            self.encoder = nn.Sequential(*layers)
        else:
            self.encoder = nn.Embedding(n_contexts + 1, output_dim)

    def forward(self, context_ids, context_features):
        if self.use_features:
            return self.encoder(context_features[context_ids])
        if self.training:
            # Train UNK explicitly so an unseen context has a learned fallback.
            unknown = torch.rand(context_ids.shape, device=context_ids.device) < self.dropout
            context_ids = context_ids.masked_fill(unknown, 0)
        return self.encoder(context_ids)


class GeneEncoder(nn.Module):
    def __init__(self, feature_dim, config):
        super().__init__()
        h = config.hidden_dim
        self.input = nn.Sequential(nn.Linear(feature_dim, h), nn.LayerNorm(h), nn.GELU())
        self.layers = nn.ModuleList(nn.Linear(h, h) for _ in range(config.graph_layers))
        self.norms = nn.ModuleList(nn.LayerNorm(h) for _ in range(config.graph_layers))
        self.dropout = config.dropout

    def forward(self, gene_features, graph):
        x = self.input(gene_features)
        for linear, norm in zip(self.layers, self.norms):
            if graph is None:
                messages = x
            else:
                # Sparse aggregation accumulates in FP32; dense layers use AMP.
                with torch.autocast(device_type=x.device.type, enabled=False):
                    messages = torch.sparse.mm(graph, x.float())
            update = F.dropout(F.gelu(linear(messages)), self.dropout, self.training)
            x = norm(x + update)
        return x


def invariant_pair(a, b):
    return torch.cat((a + b, torch.abs(a - b), a * b), dim=-1)


def invariant_singles(single_effects, single_missing, additive=False):
    observed_sum = single_effects.sum(dim=-1, keepdim=True)
    missing_count = single_missing.sum(dim=-1, keepdim=True)
    if additive:
        return torch.cat((observed_sum, missing_count), dim=-1)
    return torch.cat((observed_sum,
                      torch.abs(single_effects[:, :1] - single_effects[:, 1:]),
                      single_effects[:, :1] * single_effects[:, 1:], missing_count), dim=-1)


class PairRegressor(nn.Module):
    def __init__(self, feature_dim, n_contexts, context_features_dim, config):
        super().__init__()
        self.config = config
        self.variant = config.variant
        self.use_graph = self.variant in {"gnn", "contextless"}
        self.gene_encoder = (GeneEncoder(feature_dim, config)
                             if self.variant in {"gnn", "graphless", "contextless"} else None)
        self.context_encoder = None
        context_dim = 0
        if self.variant != "contextless":
            context_dim = config.context_dim
            self.context_encoder = ContextEncoder(
                n_contexts, context_features_dim, context_dim, config.dropout,
                linear_only=self.variant == "additive",
            )
        pair_dim = (config.hidden_dim if self.gene_encoder is not None else feature_dim)
        pair_dim *= 1 if self.variant == "additive" else 3
        if config.use_single_effects:
            pair_dim += 2 if self.variant == "additive" else 4
        pair_dim += context_dim
        if self.variant == "additive":
            self.head = nn.Linear(pair_dim, 2)
        else:
            self.head = nn.Sequential(
                nn.Linear(pair_dim, config.hidden_dim * 2), nn.GELU(), nn.Dropout(config.dropout),
                nn.Linear(config.hidden_dim * 2, config.hidden_dim), nn.GELU(),
                nn.Dropout(config.dropout), nn.Linear(config.hidden_dim, 2),
            )

    def encode_genes(self, gene_features, graph):
        if self.gene_encoder is None:
            return gene_features
        return self.gene_encoder(gene_features, graph if self.use_graph else None)

    def forward(self, gene_embeddings, pairs, context_ids, context_features,
                single_effects, single_missing):
        a, b = gene_embeddings[pairs[:, 0]], gene_embeddings[pairs[:, 1]]
        parts = [a + b if self.variant == "additive" else invariant_pair(a, b)]
        if self.config.use_single_effects:
            parts.append(invariant_singles(single_effects, single_missing, self.variant == "additive"))
        if self.context_encoder is not None:
            parts.append(self.context_encoder(context_ids, context_features))
        return self.head(torch.cat(parts, dim=-1))


class NaiadReimplementation(nn.Module):
    """Reference NAIAD MLP backbone with explicit evaluation adaptations.

    Reference: https://github.com/NeptuneBio/NAIAD, naiad/models.py,
    EmbedPhenoModel (MLP, both branches). The original 16-dimensional learned
    gene embeddings and 512-wide phenotype branch are retained. Adaptations:
    two outputs for the separately standardized tasks; averaging both single-
    effect orders for exact symmetry; unseen genes receive a fixed zero UNK.
    This ID-based comparator is not an inductive gene-feature model. It needs
    observed single effects and should be benchmarked within one context.
    """

    context_encoder = None

    def __init__(self, n_genes, known_genes, config):
        super().__init__()
        self.config = config
        self.variant = "naiad_reimplementation"
        self.use_graph = False
        d = config.naiad_embedding_dim
        p = 0.1
        lookup = torch.zeros(n_genes, dtype=torch.long)
        lookup[known_genes.cpu()] = torch.arange(1, len(known_genes) + 1)
        self.register_buffer("gene_lookup", lookup)
        self.embedding = nn.Embedding(len(known_genes) + 1, d, padding_idx=0)
        self.embedding_ffn = nn.Sequential(
            nn.Linear(d, d), nn.Dropout(p), nn.ReLU(), nn.Linear(d, d)
        )
        self.embedding_comb = nn.Sequential(
            nn.Linear(d, d), nn.Dropout(p), nn.ReLU(), nn.Linear(d, 2)
        )
        self.pheno_ffn = nn.Sequential(
            nn.Linear(2, 512, bias=False), nn.Dropout(p), nn.GELU(),
            nn.Linear(512, 512, bias=False), nn.Dropout(p), nn.GELU(), nn.Linear(512, 2),
        )

    def encode_genes(self, gene_features, graph):
        return self.embedding_ffn(self.embedding(self.gene_lookup))

    def forward(self, gene_embeddings, pairs, context_ids, context_features,
                single_effects, single_missing):
        embedded_pair = gene_embeddings[pairs].sum(dim=1)
        phenotypes = 0.5 * (self.pheno_ffn(single_effects)
                            + self.pheno_ffn(single_effects.flip(-1)))
        return self.embedding_comb(embedded_pair) + phenotypes


def build_model(feature_dim, n_genes, n_contexts, context_features_dim,
                known_genes, config):
    if config.variant == "naiad_reimplementation":
        return NaiadReimplementation(n_genes, known_genes, config)
    return PairRegressor(feature_dim, n_contexts, context_features_dim, config)
