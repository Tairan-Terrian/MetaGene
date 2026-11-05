#!/usr/bin/env python3
"""Prepare independent human gene categories and sparse biological features.

Run from any directory. Cached source files make the preparation reproducible;
--download explicitly refreshes them. No interaction-screen labels are used.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
import json
from pathlib import Path
import urllib.request
import zipfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "ReactomePathways.gmt.zip": "https://reactome.org/download/current/ReactomePathways.gmt.zip",
    "reactome_release.txt": "https://reactome.org/ContentService/data/database/version",
    "EpiGenes_main_v2.1.csv": "https://epifactors.autosome.org/public_data/v2.1/EpiGenes_main.csv",
    "hgnc_complete_set.txt": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
}


def symbol_mapping(hgnc):
    """Only resolve unique aliases; approved symbols take priority over aliases."""
    candidates = defaultdict(set)
    for row in hgnc.itertuples(index=False):
        for field in (row.alias_symbol, row.prev_symbol):
            for alias in field.split("|"):
                if alias:
                    candidates[alias].add(row.symbol)
    mapping = {alias: next(iter(genes)) for alias, genes in candidates.items() if len(genes) == 1}
    mapping.update({gene: gene for gene in hgnc.symbol})
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data/annotations")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--max-edge-pathway-size", type=int, default=50)
    args = parser.parse_args()
    out = args.out_dir
    source = out / "sources"
    source.mkdir(parents=True, exist_ok=True)
    source_records = []
    for name, url in SOURCES.items():
        path = source / name
        if args.download or not path.exists():
            urllib.request.urlretrieve(url, path)
        source_records.append({
            "file": str(path.relative_to(out)), "url": url,
            "file_saved_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
            "bytes": path.stat().st_size,
        })
    hgnc = pd.read_csv(source / "hgnc_complete_set.txt", sep="\t", dtype=str).fillna("")
    hgnc = hgnc[hgnc.status == "Approved"]
    mapping = symbol_mapping(hgnc)
    pd.DataFrame(sorted(mapping.items()), columns=["source_symbol", "gene"]).to_csv(out / "symbol_mapping.csv", index=False)
    epifactors = pd.read_csv(source / "EpiGenes_main_v2.1.csv").fillna("")
    hgnc_id_map = dict(zip(hgnc.hgnc_id, hgnc.symbol))
    epifactors["gene"] = ("HGNC:" + epifactors.HGNC_ID.astype(str)).map(hgnc_id_map)
    epifactors["gene"] = epifactors.gene.fillna(epifactors.HGNC_symbol.map(mapping))
    epi_unmapped = epifactors[epifactors.gene.isna()]
    epi_unmapped.to_csv(out / "epifactors_unmapped.csv", index=False)
    epifactors = epifactors.dropna(subset=["gene"])
    epifactors.to_csv(out / "epifactors_curated.csv", index=False)
    epi = set(epifactors.gene)

    rows = []
    pathways = []
    metabolic = set()
    unresolved = []
    with zipfile.ZipFile(source / "ReactomePathways.gmt.zip") as archive:
        lines = archive.read("ReactomePathways.gmt").decode().splitlines()
    for line in lines:
        name, pathway_id, *symbols = line.split("\t")
        genes = sorted({mapping[symbol] for symbol in symbols if symbol in mapping})
        unresolved.extend((pathway_id, symbol) for symbol in symbols if symbol not in mapping)
        pathways.append((pathway_id, name, genes))
        rows.extend((gene, pathway_id, name) for gene in genes)
        if pathway_id == "R-HSA-1430728":
            metabolic.update(genes)
    membership = pd.DataFrame(rows, columns=["gene", "pathway_id", "pathway_name"])
    membership.to_csv(out / "pathway_membership.csv", index=False)
    pd.DataFrame(unresolved, columns=["pathway_id", "unmapped_symbol"]).to_csv(out / "reactome_unmapped.csv", index=False)
    pd.DataFrame([(pid, name, len(genes)) for pid, name, genes in pathways],
                 columns=["pathway_id", "pathway_name", "n_genes"]).to_csv(out / "pathways.csv", index=False)

    # Edges encode shared curated function, not measured physical interactions.
    # Compact pathways bound enumeration; no genome-wide dense pair matrix exists.
    edges = defaultdict(lambda: [0, 0.0])
    for _, _, genes in pathways:
        if 2 <= len(genes) <= args.max_edge_pathway_size:
            for pair in combinations(genes, 2):
                edges[pair][0] += 1
                edges[pair][1] += 1.0 / (len(genes) - 1)
    pd.DataFrame([(a, b, count, weight) for (a, b), (count, weight) in sorted(edges.items())],
                 columns=["gene_a", "gene_b", "shared_pathways", "weight"]).to_csv(out / "functional_edges.csv", index=False)
    annotation = hgnc[["symbol", "hgnc_id", "entrez_id", "ensembl_gene_id", "locus_group"]].rename(columns={"symbol": "gene"})
    annotation["is_metabolic"] = annotation.gene.isin(metabolic).astype(int)
    annotation["is_epigenetic"] = annotation.gene.isin(epi).astype(int)
    annotation["n_reactome_pathways"] = annotation.gene.map(membership.groupby("gene").size()).fillna(0).astype(int)
    annotation["metabolic_source"] = annotation.is_metabolic.map({1: "Reactome:R-HSA-1430728", 0: ""})
    annotation["epigenetic_source"] = annotation.is_epigenetic.map({1: "EpiFactors:2.1:proteins", 0: ""})
    annotation.to_csv(out / "gene_annotations.csv", index=False)
    manifest = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "sources": source_records,
        "reactome_release": (source / "reactome_release.txt").read_text().strip(),
        "epifactors_version": "2.1, September 10 2024",
        "epifactors_citation": "https://doi.org/10.1093/nar/gkac989",
        "epifactors_description": "https://epifactors.autosome.org/description",
        "hgnc_version": "Complete approved gene snapshot; source file modification time above",
        "metabolic_definition": "All HGNC-resolved genes in Reactome Metabolism R-HSA-1430728; broad metabolic pathway membership, not an enzyme-only label.",
        "epigenetic_definition": "HGNC-resolved proteins in curated EpiFactors 2.1 main protein table; includes cofactors and RNA modification factors; not limited to chromatin enzymes.",
        "graph_definition": f"Undirected gene pairs sharing a Reactome pathway of 2-{args.max_edge_pathway_size} human HGNC-resolved genes; weight=sum(1/(pathway_size-1)). These are functional co-membership edges, not physical PPI evidence. Stored once with gene_a < gene_b; symmetrize for message passing.",
        "feature_definition": "Binary sparse gene-pathway membership across full Reactome vocabulary. Features are available for unseen screen genes independently of labels. Unknown genes remain unannotated, not inferred from their name.",
        "symbol_resolution": "EpiFactors proteins use their HGNC IDs first. For symbol inputs, approved HGNC symbols take precedence; previous and alias symbols are resolved only if unique; unresolved entries exported separately.",
        "counts": {"hgnc_genes": len(annotation), "pathways": len(pathways), "pathway_memberships": len(membership), "functional_edges": len(edges), "metabolic_genes": len(metabolic), "epigenetic_genes": len(epi), "overlap": len(metabolic & epi), "epifactors_unmapped": len(epi_unmapped)},
    }
    (out / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
