"""Independent PPI inputs and observed-single evidence for the current analysis."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    root = args.root
    out = root / "data/evidence"
    sources = out / "sources"
    genes = pd.read_csv(root / "data/processed/model_genes.csv")
    index = {g: i for i, g in enumerate(genes.gene)}
    entrez = {str(int(r.entrez_id)): r.gene for r in genes.dropna(subset=["entrez_id"]).itertuples()}
    mapping = pd.read_csv(root / "data/annotations/symbol_mapping.csv", keep_default_na=False).set_index("source_symbol").gene.to_dict()
    graphs = {"reactome": np.load(root / "data/processed/model_arrays.npz")["edge_index"]}
    graph_counts = {}
    for name, file in [("bioplex", "BioPlex_HCT116_Network_5.5K_Dec_2019.tsv"),
                       ("bioplex_293t", "BioPlex_293T_Network_10K_Dec_2019.tsv")]:
        raw = pd.read_csv(sources / file, sep="\t", dtype={"GeneA": str, "GeneB": str})
        raw["gene_a"] = raw.GeneA.map(entrez)
        raw["gene_b"] = raw.GeneB.map(entrez)
        retained = raw.dropna(subset=["gene_a", "gene_b"]).copy()
        retained = retained[retained.gene_a.ne(retained.gene_b)]
        retained[["gene_a", "gene_b"]] = np.sort(retained[["gene_a", "gene_b"]].to_numpy(), axis=1)
        retained = retained.groupby(["gene_a", "gene_b"], as_index=False).pInt.max()
        retained.to_csv(out / f"{name}_edges.csv", index=False)
        graphs[name] = retained[["gene_a", "gene_b"]].apply(lambda col: col.map(index)).to_numpy(np.int64).T
        covered = set(retained.gene_a) | set(retained.gene_b)
        graph_counts[name] = {"source_edges": len(raw), "retained_undirected_edges": len(retained),
                              "covered_model_genes": len(covered)}
    np.savez_compressed(out / "graphs.npz", **graphs)

    # R preserves the deposited names and dimensions of these public RDS objects.
    r_program = '''args <- commandArgs(trailingOnly=TRUE)
    s <- readRDS(file.path(args[1], "SCORE_median_lfc.rds"))
    write.csv(data.frame(gene=names(s), score_ht29_lfc=as.numeric(s)), file.path(args[2], "score_ht29_singles.csv"), row.names=FALSE)
    s <- readRDS(file.path(args[1], "colo1_plm_scale.rds"))
    write.csv(data.frame(guide_pair_id=rownames(s), s, check.names=FALSE), file.path(args[2], "colo1_scaled_guide_effects.csv"), row.names=FALSE)
    '''
    subprocess.run(["Rscript", "-e", r_program, str(sources.resolve()), str(out.resolve())], check=True)
    scaled = pd.read_csv(out / "colo1_scaled_guide_effects.csv").set_index("guide_pair_id")
    scaled["scaled_effect"] = scaled.mean(axis=1)
    annotation = pd.read_csv(root / "data/pairs_sources/HT29_annotations.csv", keep_default_na=False)
    annotation["gene_a"] = annotation.Gene1.map(mapping).fillna(annotation.Gene1)
    annotation["gene_b"] = annotation.Gene2.map(mapping).fillna(annotation.Gene2)
    table = annotation.merge(scaled[["scaled_effect"]], left_on="ID", right_index=True)
    controls = {"ADAD1", "CYLC2", "KLK12"}
    gi = table.MyNote.eq("GIControlsSingletons")
    anchor = table[(table.MyNote.eq("AnchorSingletons") | gi) & table.gene_a.isin(controls)]
    library = table[(table.MyNote.eq("LibrarySingletons") | gi) & table.gene_b.isin(controls)]
    first = anchor.groupby("gene_b").scaled_effect.agg(["median", "count"]).rename_axis("gene")
    second = library.groupby("gene_a").scaled_effect.agg(["median", "count"]).rename_axis("gene")
    single = pd.concat([first, second])
    single = single.loc[~single.index.duplicated(keep="first")].rename(columns={
        "median": "colo1_single_control_scaled", "count": "n_guide_pairs_in_single_estimate"})
    single["strong_single_depletion_scaled"] = single.colo1_single_control_scaled.le(-0.5)
    score = pd.read_csv(out / "score_ht29_singles.csv", keep_default_na=False)
    score["gene"] = score.gene.map(mapping).fillna(score.gene)
    score = score.groupby("gene").score_ht29_lfc.median().to_frame()
    score["score_ht29_relative_rank"] = score.score_ht29_lfc.rank(pct=True)
    evidence = single.join(score, how="outer")
    evidence.to_csv(out / "single_gene_evidence.csv")
    control_values = table[table.Note.isin(["NegativeControls", "PositiveControls"])].groupby("Note").scaled_effect.agg(["median", "count"])
    meta = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "graphs": graph_counts,
        "bioplex_publication": "https://doi.org/10.1016/j.cell.2021.04.011",
        "bioplex_source_page": "https://bioplex.hms.harvard.edu/interactions.php",
        "graph_definition": "Released, filtered AP-MS protein associations. HCT116 is the prespecified current PPI graph; 293T is retained as evidence. Model uses binary edges plus self loops and symmetric degree normalization. AP-MS co-complex association does not establish direct binding.",
        "graph_mapping": "Entrez gene IDs to the existing HGNC model gene table; all graph sources independent of combinatorial phenotype labels.",
        "graph_learning_scope": "Complete external annotation graph is fixed during training and inference; held-out pair/gene phenotype labels remain excluded. Annotation-transductive, label-inductive evaluation.",
        "single_source": "https://github.com/ibarrioh/DualGuide_COLO1/tree/main/Analysis/ENCORE/input",
        "score_source_definition": "Authors' deposited SCORE_median_lfc: independent HT29 gene-level median guide LFC downloaded by the original study from Project SCORE; original release date and exact control scaling not specified in its readme. Retain raw LFC and within-file relative rank.",
        "scaled_source_definition": "Author-deposited colo1_plm_scale guide-pair values, negative and positive controls scaled to 0 and -1. Mean available biological replicate values per guide pair, then median across gene-control guide pairs; anchor singleton position preferred when both roles are available.",
        "strong_single_rule": "Control-scaled singleton median <= -0.5: an explicit midpoint depletion filter, not a formal essentiality probability or proof of tumor selectivity. Sensitivity to -0.3/-0.7 is reported in candidate outputs.",
        "single_evidence_counts": {"control_scaled_genes": int(single.colo1_single_control_scaled.notna().sum()),
                                   "independent_score_genes": len(score)},
        "deposited_control_summaries": control_values.reset_index().to_dict("records"),
        "model_target": "Original joint and additive residual labels are retained for directly comparable benchmarking; control-scaled singles are an independent candidate triage layer.",
    }
    (out / "feature_manifest.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
