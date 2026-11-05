#!/usr/bin/env python3
"""Reanalyse independent HT29/HTOXAR3 expression evidence from GSE76092.

The supplied array matrix already has normexp background correction, log2
transformation, quantile normalization, and gene-identifier collapse. It is
neither a count matrix nor an absolute expression/detection measurement.
"""
import argparse
from datetime import datetime, timezone
import gzip
import io
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def prepare(args):
    args.out_dir.mkdir(parents=True, exist_ok=True)
    symbols = pd.read_csv(args.annotation_dir / "symbol_mapping.csv", keep_default_na=False, dtype=str)
    symbol_map = dict(zip(symbols.source_symbol, symbols.gene))
    approved = set(pd.read_csv(args.annotation_dir / "gene_annotations.csv").gene)
    platform_text = (args.source_dir / "GPL21253.txt").read_text()
    platform = pd.read_csv(io.StringIO(platform_text.split("!platform_table_begin\n")[1]
                                      .split("!platform_table_end")[0]), sep="\t", dtype=str,
                           keep_default_na=False)
    platform = platform.rename(columns={"ID": "ID_REF", "GeneName": "source_symbol"})
    platform["gene"] = platform.source_symbol.map(symbol_map)
    platform["included"] = platform.gene.isin(approved)
    platform[["ID_REF", "source_symbol", "gene", "included"]].to_csv(
        args.out_dir / "GSE76092_symbol_mapping.csv", index=False)
    with gzip.open(args.source_dir / "GSE76092_series_matrix.txt.gz", "rt") as stream:
        matrix_text = stream.read()
    matrix = pd.read_csv(io.StringIO(matrix_text.split("!series_matrix_table_begin\n")[1]
                                    .split("!series_matrix_table_end")[0]), sep="\t",
                         dtype={"ID_REF": str})
    mapped = matrix.merge(platform.loc[platform.included, ["ID_REF", "gene"]], on="ID_REF")
    expression = mapped.drop(columns="ID_REF").groupby("gene").median(numeric_only=True)
    expression.to_csv(args.out_dir / "GSE76092_log2_intensity.csv")

    sample_rows = []
    for block in (args.source_dir / "GSE76092_samples.txt").read_text().split("^SAMPLE = ")[1:]:
        sample = block.splitlines()[0].strip()
        title = re.search(r"!Sample_title = (.*)", block).group(1).strip()
        cell_line = title.split("_")[0]
        treatment = ("untreated" if "untreated" in title else
                     "OXA_CUR" if "plus_curcumin" in title else "OXA")
        sample_rows.append({
            "sample": sample, "title": title, "cell_line": cell_line,
            "resistance_state": "parental" if cell_line == "HT29" else "acquired_resistant",
            "treatment": treatment, "group": f"{cell_line}_{treatment}",
            "biological_replicate": int(title.rsplit("repl", 1)[1]),
            "exposure_hours": 0 if treatment == "untreated" else 24,
            "oxaliplatin_uM": 0 if treatment == "untreated" else 10 if cell_line == "HT29" else 30,
            "curcumin_uM": 10 if treatment == "OXA_CUR" else 0,
            "primary_analysis": treatment == "untreated",
            "secondary_analysis": treatment in {"untreated", "OXA"},
        })
    metadata = pd.DataFrame(sample_rows).set_index("sample").loc[expression.columns]
    metadata.to_csv(args.out_dir / "GSE76092_sample_mapping.csv")
    expression.corr().to_csv(args.out_dir / "GSE76092_sample_correlations.csv")
    centered = expression.T - expression.T.mean(axis=0)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    pc = pd.DataFrame(u[:, :3] * singular[:3], index=expression.columns,
                      columns=["PC1", "PC2", "PC3"])
    metadata.join(pc).to_csv(args.out_dir / "GSE76092_sample_pca.csv")

    intensity = pd.DataFrame(index=expression.index)
    for group, samples in metadata.groupby("group"):
        block = expression[samples.index]
        intensity[f"{group}_mean_log2_intensity"] = block.mean(axis=1)
        intensity[f"{group}_sd_log2_intensity"] = block.std(axis=1)
        intensity[f"{group}_relative_intensity_percentile"] = block.mean(axis=1).rank(pct=True)
    intensity["interpretation"] = "relative array intensity; not absolute expression or a detection call"
    intensity.to_csv(args.out_dir / "GSE76092_gene_intensity.csv")

    pathways = pd.read_csv(args.annotation_dir / "pathway_membership.csv")
    pathways = pathways[pathways.gene.isin(expression.index)].drop_duplicates(["gene", "pathway_id"])
    sizes = pathways.groupby("pathway_id").gene.nunique()
    eligible = sizes[(sizes >= 10) & (sizes <= 500)].index
    pathways[pathways.pathway_id.isin(eligible)].to_csv(args.out_dir / "reactome_tested_membership.csv", index=False)
    source_info = {
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "study": "GSE76092", "platform": "GPL21253", "publication_pmid": "27091625",
        "urls": {
            "series": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092",
            "platform": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21253",
            "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76092/matrix/GSE76092_series_matrix.txt.gz",
            "sample_metadata": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76092&targ=gsm&form=text&view=brief",
        },
        "design": "2 cell states x 3 treatments x 3 independent biological replicates = 18 samples",
        "primary_contrast": "untreated HTOXAR3 minus untreated HT29; independent 3 versus 3 fit",
        "secondary_contrasts": ["HT29 OXA minus untreated", "HTOXAR3 OXA minus untreated",
                                "(HTOXAR3 OXA - untreated) - (HT29 OXA - untreated)"],
        "secondary_fit": "four groups, 12 samples; state by drug response interaction",
        "dose_caveat": "HT29 10 uM and HTOXAR3 30 uM oxaliplatin for 24 h, each at its own IC50; interaction is dose-confounded and is not a same-dose sensitization effect",
        "curcumin_caveat": "OXA+curcumin samples are retained in exported matrices and QC but excluded from primary and secondary inference; no curcumin-only arm exists",
        "primary_interpretation": "association with acquired resistance; not a causal sensitization or gene-knockout effect",
        "supplied_preprocessing": "normexp background correction; log2; quantile normalization; common gene-identifier collapse; control probes excluded by authors",
        "additional_preprocessing": "HGNC symbol canonicalization; retain approved symbols; median collapse only where identifiers map to the same current symbol; no re-normalization",
        "intensity_interpretation": "within-array relative intensity percentile, not absolute expression or measured present/absent calls",
        "n_input_identifiers": len(matrix), "n_mapped_identifiers": len(mapped),
        "n_approved_genes": len(expression), "n_samples": len(metadata),
        "pca_variance_explained": (singular[:3] ** 2 / np.sum(singular ** 2)).tolist(),
        "n_tested_reactome_pathways": len(eligible),
        "pathway_universe": "all HGNC-mapped genes tested in the corresponding expression fit; pathway size 10-500 within this universe",
        "differential_method": "limma lmFit + contrasts.fit + eBayes(trend=TRUE, robust=TRUE)",
        "multiple_testing": "BH and BY within each contrast across all tested genes; original publication used BY",
        "pathway_methods": "CAMERA competitive gene-set test (inter.gene.cor=0.01) and directional hypergeometric ORA (BH<0.05, |log2FC|>=1); pathway BH correction within contrast and direction",
        "independent_evidence_role": "disease prioritization only; no pairwise model training, tuning, or outcome relabeling",
        "inference_status": "prepared_only",
    }
    (args.out_dir / "source_manifest.json").write_text(json.dumps(source_info, indent=2) + "\n")
    print(json.dumps({k: source_info[k] for k in ["n_input_identifiers", "n_mapped_identifiers", "n_approved_genes", "n_samples", "n_tested_reactome_pathways"]}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=ROOT / "data/evidence/sources")
    parser.add_argument("--annotation-dir", type=Path, default=ROOT / "data/annotations")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "results/disease")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--r-lib", default="")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    prepare(args)
    if not args.prepare_only:
        subprocess.run([args.rscript, str(Path(__file__).with_name("fit_resistance.R")),
                        str(args.out_dir.resolve()), args.r_lib], check=True)
        manifest_path = args.out_dir / "source_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["inference_status"] = "limma_completed"
        manifest["contrast_summary"] = pd.read_csv(args.out_dir / "contrast_summary.csv").to_dict("records")
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
