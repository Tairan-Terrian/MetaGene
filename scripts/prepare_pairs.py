"""Prepare measured, distinct-gene combinatorial perturbation labels.

Horlbeck/Norman follow the public NAIAD data-preparation conventions.
HT29 follows the deposited gene-level LFC and singleton aggregation in
DualGuide_COLO1/Analysis/ENCORE/scripts/Rscript_final.R, lines 941-1017.
See data/pairs_sources/source_definitions.json for the label definitions.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
HT29_CONTROL_GENES = {"ADAD1", "CYLC2", "KLK12"}


def convert_ht29_rds(raw_dir, source_dir):
    """Export the author's numeric dataframe and character matrix losslessly."""
    import rdata

    lfc = rdata.read_rds(raw_dir / "HT29_gene_lfc.rds")
    lfc.to_csv(source_dir / "HT29_gene_lfc.csv", index_label="source_key")
    parsed = rdata.parser.parse_file(raw_dir / "HT29_annotations.rds")
    attributes = rdata.conversion.convert(parsed.object.attributes)
    annotation = pd.DataFrame(
        rdata.conversion.convert(parsed).reshape(tuple(attributes["dim"]), order="F"),
        columns=attributes["dimnames"][1],
    )
    annotation.to_csv(source_dir / "HT29_annotations.csv", index=False)


def normalize_symbols(values, mapping):
    return values.map(mapping).fillna(values)


def canonicalize(frame):
    frame = frame.copy()
    frame[["gene_a", "gene_b"]] = np.sort(frame[["gene_a", "gene_b"]].to_numpy(), axis=1)
    return frame


def join_labels(pairs, singles, dataset, perturbation, label_type, author_gi=None):
    pairs = canonicalize(pairs)
    pairs = pairs.groupby(["gene_a", "gene_b"], as_index=False).agg(
        joint_effect=("joint_effect", "mean"),
        n_source_records=("n_source_records", "sum"),
        source_pair_ids=("source_pair_ids", lambda x: "|".join(sorted(set(x)))),
    )
    singles = singles.groupby("gene", as_index=False).single_effect.mean()
    single_map = singles.set_index("gene").single_effect
    pairs["single_a"] = pairs.gene_a.map(single_map)
    pairs["single_b"] = pairs.gene_b.map(single_map)
    pairs["exclusion_reason"] = ""
    pairs.loc[pairs.gene_a == pairs.gene_b, "exclusion_reason"] = "same_gene_dual_guide"
    pairs.loc[pairs.joint_effect.isna(), "exclusion_reason"] = "joint_effect_unmeasured"
    missing_single = pairs[["single_a", "single_b"]].isna().any(axis=1)
    pairs.loc[missing_single & pairs.exclusion_reason.eq(""), "exclusion_reason"] = "single_effect_unmeasured"
    excluded = pairs.loc[pairs.exclusion_reason.ne("")].copy()
    pairs = pairs.loc[pairs.exclusion_reason.eq("")].drop(columns="exclusion_reason")
    pairs["author_gi"] = np.nan
    if author_gi is not None:
        author_gi = canonicalize(author_gi).groupby(["gene_a", "gene_b"]).author_gi.mean()
        keys = pd.MultiIndex.from_frame(pairs[["gene_a", "gene_b"]])
        pairs["author_gi"] = author_gi.reindex(keys).to_numpy()
    pairs["additive_gi_proxy"] = pairs.joint_effect - pairs.single_a - pairs.single_b
    pairs["hsa_gi_proxy"] = pairs.joint_effect - pairs[["single_a", "single_b"]].min(axis=1)
    pairs["dataset"] = dataset
    pairs["perturbation_type"] = perturbation
    pairs["label_type"] = label_type
    pairs["label_direction"] = "lower_is_less_growth"
    singles = singles.dropna(subset=["single_effect"])
    singles["dataset"] = dataset
    singles["perturbation_type"] = perturbation
    excluded["dataset"] = dataset
    return pairs, singles, excluded


def prepare_naiad(raw, dataset, perturbation, mapping, author_gi_raw=None):
    raw = raw.copy()
    raw["source_pair_ids"] = raw.id1 + "+" + raw.id2
    raw["gene_a"] = normalize_symbols(raw.id1, mapping)
    raw["gene_b"] = normalize_symbols(raw.id2, mapping)
    raw = canonicalize(raw)
    has_control = raw.gene_a.eq("negative") | raw.gene_b.eq("negative")
    single_rows = raw.loc[has_control & raw.gene_a.ne(raw.gene_b)].copy()
    single_rows["gene"] = np.where(single_rows.gene_a.eq("negative"), single_rows.gene_b, single_rows.gene_a)
    singles = single_rows.groupby("gene", as_index=False).score.mean().rename(columns={"score": "single_effect"})
    pair_rows = raw.loc[~has_control].rename(columns={"score": "joint_effect"})
    pair_rows["n_source_records"] = 1
    author_gi = None
    if author_gi_raw is not None:
        author_gi = pd.DataFrame({
            "gene_a": normalize_symbols(author_gi_raw.id1, mapping),
            "gene_b": normalize_symbols(author_gi_raw.id2, mapping),
            "author_gi": author_gi_raw.score,
        })
    pairs, singles, excluded = join_labels(pair_rows, singles, dataset, perturbation, "growth_gamma", author_gi)
    singles["single_definition"] = "mean_gene_plus_negative_control_gamma"
    summary = {
        "input_rows": len(raw),
        "input_single_control_rows": len(single_rows),
        "input_negative_negative_rows": int((raw.gene_a.eq("negative") & raw.gene_b.eq("negative")).sum()),
        "input_same_gene_dual_guide_rows": int((~has_control & raw.gene_a.eq(raw.gene_b)).sum()),
        "input_missing_scores": int(raw.score.isna().sum()),
        "aggregation": "Mean measured scores across guide pairs and swapped orientations after symbol normalization. Horlbeck CSVs are already gene-pair aggregated; Norman gamma contains guide-pair records.",
    }
    return pairs, singles, excluded, summary


def prepare_ht29(source_dir, mapping):
    lfc = pd.read_csv(source_dir / "HT29_gene_lfc.csv", index_col=0)
    annotation = pd.read_csv(source_dir / "HT29_annotations.csv")
    annotation["source_key"] = annotation.Note + "_" + annotation.Gene1 + "_" + annotation.Gene2
    design = annotation.drop_duplicates("source_key").set_index("source_key")
    table = lfc.join(design)
    table["joint_effect"] = lfc.mean(axis=1)
    table["gene_a"] = normalize_symbols(table.Gene1, mapping)
    table["gene_b"] = normalize_symbols(table.Gene2, mapping)
    table["source_pair_ids"] = table.Gene1 + "+" + table.Gene2
    control_a = table.Gene1.isin(HT29_CONTROL_GENES)
    control_b = table.Gene2.isin(HT29_CONTROL_GENES)
    gi_single = table.MyNote.eq("GIControlsSingletons")
    anchor = table.loc[table.MyNote.eq("AnchorSingletons") | (gi_single & control_a)]
    library = table.loc[table.MyNote.eq("LibrarySingletons") | (gi_single & control_b)]
    anchor_single = anchor.groupby("gene_b").joint_effect.median()
    library_single = library.groupby("gene_a").joint_effect.median()
    # The authors concatenate anchor values before library values; R named
    # indexing selects the first match for genes present in both positions.
    single = pd.concat([anchor_single, library_single])
    single = single.loc[~single.index.duplicated(keep="first")]
    singles = single.rename_axis("gene").reset_index(name="single_effect")
    pair_rows = table.loc[table.MyNote.isin([
        "AnchorCombinations", "LibraryCombinations", "GIControlsCombinations",
    ])].copy()
    pair_rows["n_source_records"] = 1
    pairs, singles, excluded = join_labels(pair_rows, singles, "ht29_colo1", "CRISPRko", "log2_fold_change_vs_plasmid")
    singles["single_definition"] = "median_control_pair_replicate_means_anchor_position_preferred"
    support = table.loc[pair_rows.index, ["gene_a", "gene_b"]].copy()
    support["n_biological_replicates"] = lfc.loc[support.index].notna().sum(axis=1)
    support["n_guide_pairs_in_design"] = annotation.groupby("source_key").size().reindex(support.index)
    support = canonicalize(support).groupby(["gene_a", "gene_b"], as_index=False).agg(
        n_biological_replicates_min=("n_biological_replicates", "min"),
        n_biological_replicates_max=("n_biological_replicates", "max"),
        n_guide_pairs_in_design=("n_guide_pairs_in_design", "sum"),
    )
    pairs = pairs.merge(support, on=["gene_a", "gene_b"], how="left")
    summary = {
        "input_gene_level_rows": len(lfc),
        "input_guide_pair_annotation_rows": len(annotation),
        "input_biological_replicates": list(lfc.columns),
        "input_combination_designs": len(pair_rows),
        "input_combination_designs_without_any_measured_replicate": int(pair_rows.joint_effect.isna().sum()),
        "input_same_gene_combination_designs": int(pair_rows.gene_a.eq(pair_rows.gene_b).sum()),
        "aggregation": "Use deposited gene-level LFC; mean available biological replicates per orientation, then mean observed orientations of each unordered gene pair. Singleton controls: median of replicate means across ADAD1/CYLC2/KLK12 control arms; anchor position preferred when both positions exist, matching author R code.",
        "author_gi_status": "Unavailable: final article uses GEMINI, but deposited main and Release2 analysis supply additive/HSA calculations rather than GEMINI result tables. author_gi is intentionally empty.",
    }
    return pairs, singles, excluded, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=PROJECT / "data/raw")
    parser.add_argument("--source-dir", type=Path, default=PROJECT / "data/pairs_sources")
    parser.add_argument("--output-dir", type=Path, default=PROJECT / "data/processed")
    parser.add_argument("--gene-map", type=Path, default=PROJECT / "data/annotations/symbol_mapping.csv")
    parser.add_argument("--convert-rds", action="store_true", help="Regenerate HT29 CSV exports from author RDS inputs; requires rdata.")
    args = parser.parse_args()
    if args.convert_rds:
        convert_ht29_rds(args.raw_dir, args.source_dir)
    gene_map = pd.read_csv(args.gene_map, keep_default_na=False)
    mapping = gene_map.set_index("source_symbol").gene.to_dict()
    results = []
    for cell in ["k562", "jurkat"]:
        raw = pd.read_csv(args.raw_dir / f"horlbeck_{cell}_gamma.csv")
        results.append(prepare_naiad(raw, f"horlbeck_{cell}", "CRISPRi", mapping))
    results.append(prepare_naiad(
        pd.read_csv(args.source_dir / "norman_gamma.csv"), "norman_k562", "CRISPRa", mapping,
        pd.read_csv(args.raw_dir / "norman_gi.csv"),
    ))
    results.append(prepare_ht29(args.source_dir, mapping))
    pairs = pd.concat([result[0] for result in results], ignore_index=True)
    singles = pd.concat([result[1] for result in results], ignore_index=True)
    excluded = pd.concat([result[2] for result in results], ignore_index=True)
    columns = ["dataset", "gene_a", "gene_b", "joint_effect", "single_a", "single_b", "author_gi", "additive_gi_proxy", "hsa_gi_proxy", "label_type", "perturbation_type", "label_direction", "n_source_records", "source_pair_ids", "n_biological_replicates_min", "n_biological_replicates_max", "n_guide_pairs_in_design"]
    pairs = pairs[columns].sort_values(["dataset", "gene_a", "gene_b"])
    singles = singles[["dataset", "gene", "single_effect", "perturbation_type", "single_definition"]].sort_values(["dataset", "gene"])
    summary = {}
    for pair_table, single_table, excluded_table, metadata in results:
        dataset = pair_table.dataset.iloc[0]
        metadata.update({
            "measured_distinct_gene_pairs": len(pair_table),
            "measured_single_genes": len(single_table),
            "pair_genes": len(set(pair_table.gene_a) | set(pair_table.gene_b)),
            "exclusions": excluded_table.exclusion_reason.value_counts().to_dict(),
            "author_gi_available_pairs": int(pair_table.author_gi.notna().sum()),
            "joint_effect_min": float(pair_table.joint_effect.min()),
            "joint_effect_max": float(pair_table.joint_effect.max()),
            "joint_effect_mean": float(pair_table.joint_effect.mean()),
            "perturbation_type": pair_table.perturbation_type.iloc[0],
            "label_type": pair_table.label_type.iloc[0],
        })
        summary[dataset] = metadata
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.output_dir / "pairs.csv", index=False)
    singles.to_csv(args.output_dir / "singles.csv", index=False)
    excluded.to_csv(args.output_dir / "excluded_pair_designs.csv", index=False)
    (args.output_dir / "dataset_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    quality = {
        "total_measured_distinct_pairs": len(pairs),
        "duplicate_dataset_gene_pairs": int(pairs.duplicated(["dataset", "gene_a", "gene_b"]).sum()),
        "nonfinite_primary_labels": int((~np.isfinite(pairs[["joint_effect", "single_a", "single_b"]])).sum().sum()),
        "noncanonical_or_self_rows": int(pairs.gene_a.ge(pairs.gene_b).sum()),
        "single_rows": len(singles),
        "duplicate_single_genes_per_dataset": int(singles.duplicated(["dataset", "gene"]).sum()),
        "missing_as_negative_labels_added": 0,
    }
    (args.output_dir / "preparation_quality.json").write_text(json.dumps(quality, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
