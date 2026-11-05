"""Assemble measured perturbations and independent functional annotations."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CONTEXTS = ["horlbeck_k562", "horlbeck_jurkat", "ht29_colo1"]
CONTEXT_FEATURES = np.array([
    [1, 0, 0, 1, 0],  # myeloid; CRISPRi
    [0, 1, 0, 1, 0],  # T-cell; CRISPRi
    [0, 0, 1, 0, 1],  # colorectal; CRISPR knockout
], dtype=np.float32)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    data_dir = args.root / "data"
    out = data_dir / "processed"
    pairs = pd.read_csv(out / "pairs.csv")
    pairs = pairs[pairs.dataset.isin(CONTEXTS)].dropna(subset=["joint_effect", "single_a", "single_b"])
    pairs = pairs.sort_values(["dataset", "gene_a", "gene_b"]).reset_index(drop=True)
    genes = sorted(set(pairs.gene_a) | set(pairs.gene_b))
    gene_index = {gene: i for i, gene in enumerate(genes)}
    membership = pd.read_csv(data_dir / "annotations/pathway_membership.csv")
    selected = membership[membership.gene.isin(genes)]
    pathways = sorted(selected.pathway_id.unique())
    pathway_index = {path: i for i, path in enumerate(pathways)}
    features = np.zeros((len(genes), len(pathways) + 3), dtype=np.float32)
    features[selected.gene.map(gene_index), selected.pathway_id.map(pathway_index)] = 1
    annotations = pd.read_csv(data_dir / "annotations/gene_annotations.csv").set_index("gene")
    gene_table = annotations.reindex(genes).reset_index()
    gene_table["is_metabolic"] = gene_table.is_metabolic.fillna(0).astype(int)
    gene_table["is_epigenetic"] = gene_table.is_epigenetic.fillna(0).astype(int)
    features[:, -3] = gene_table.is_metabolic
    features[:, -2] = gene_table.is_epigenetic
    features[:, -1] = (features[:, :-3].sum(axis=1) > 0).astype(np.float32)
    edges = pd.read_csv(data_dir / "annotations/functional_edges.csv")
    edges = edges[edges.gene_a.isin(genes) & edges.gene_b.isin(genes)]
    forward = np.stack([edges.gene_a.map(gene_index), edges.gene_b.map(gene_index)])
    edge_index = np.concatenate([forward, forward[::-1]], axis=1).astype(np.int64)
    dataset_index = {name: i for i, name in enumerate(CONTEXTS)}
    pair_ids = pairs[["gene_a", "gene_b"]].apply(lambda column: column.map(gene_index)).to_numpy(dtype=np.int64)
    context_ids = pairs.dataset.map(dataset_index).to_numpy(dtype=np.int64)
    pairs["row_id"] = np.arange(len(pairs))
    pairs["is_metabolic_epigenetic"] = (
        ((features[pair_ids[:, 0], -3] == 1) & (features[pair_ids[:, 1], -2] == 1)) |
        ((features[pair_ids[:, 1], -3] == 1) & (features[pair_ids[:, 0], -2] == 1))
    )
    np.savez_compressed(
        out / "model_arrays.npz", gene_features=features, edge_index=edge_index,
        pairs=pair_ids, context_ids=context_ids, source_ids=context_ids,
        joint=pairs.joint_effect.to_numpy(dtype=np.float32),
        single_effects=pairs[["single_a", "single_b"]].to_numpy(dtype=np.float32),
        gi_proxy=pairs.additive_gi_proxy.to_numpy(dtype=np.float32),
        context_features=CONTEXT_FEATURES,
    )
    pairs.to_csv(out / "model_pairs.csv", index=False)
    gene_table.to_csv(out / "model_genes.csv", index=False)
    pd.DataFrame({"pathway_id": pathways}).to_csv(out / "model_pathways.csv", index=False)
    summary = {
        "contexts": CONTEXTS,
        "context_features": ["myeloid", "T_cell", "colorectal", "CRISPRi", "CRISPR_KO"],
        "secondary_target": "additive_gi_proxy = joint_effect - single_a - single_b",
        "secondary_target_scope": "Source-scale additive reference; genetic-interaction prioritization proxy",
        "gene_count": len(genes), "pathway_feature_count": len(pathways),
        "directed_functional_edges": edge_index.shape[1],
        "labelled_pair_records": len(pairs),
        "pairs_per_context": pairs.groupby("dataset").size().to_dict(),
        "metabolic_epigenetic_pairs_per_context": pairs[pairs.is_metabolic_epigenetic].groupby("dataset").size().to_dict(),
        "gene_features_without_reactome": int((features[:, -1] == 0).sum()),
        "norman_role": "CRISPRa data prepared separately; primary inhibition/knockout model uses three contexts above",
    }
    (out / "model_data_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
