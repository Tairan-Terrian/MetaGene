#!/usr/bin/env python3
"""Prepare independent In4mer and Dede cohorts for frozen external evaluation.

No external pair or single-gene phenotype belongs in fitting or model selection.
Non-hits are measured pairs that do not satisfy the source effect-size rule;
they are not proven biological nulls. The sources do not supply pairwise FDR.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN4MER_COLUMNS = {"k562_T8": "K562", "A549_T21": "A549",
                  "MeljusoCP1905_T21": "MELJUSO", "A375CP1906_T21": "A375"}


def in4mer_pairs(source, symbol_map, approved):
    table = pd.read_csv(source / "In4mer_Supp_table7.txt", sep="\t", keep_default_na=False)
    table[list(IN4MER_COLUMNS)] = table[list(IN4MER_COLUMNS)].apply(pd.to_numeric, errors="coerce")
    table = table.set_index("Gene")
    pair_ids = [x for x in table.index if x.count("_") == 1 and "," not in x]
    single_ids = [x for x in table.index if "_" not in x and "," not in x]
    records, singles = [], []
    for column, cell_line in IN4MER_COLUMNS.items():
        cohort = f"in4mer2024_{cell_line.lower()}"
        for g in single_ids:
            mapped = symbol_map.get(g, g)
            value = table.loc[g, column]
            if mapped in approved and pd.notna(value):
                singles.append({"dataset": cohort, "study": "In4mer2024", "cell_line": cell_line,
                                "gene": mapped, "source_symbol": g, "single_effect": value,
                                "source": "Supplementary table 7", "effect_units": "log2 fold change"})
        for source_pair in pair_ids:
            ga, gb = source_pair.split("_")
            a, b = symbol_map.get(ga, ga), symbol_map.get(gb, gb)
            if a not in approved or b not in approved or ga not in table.index or gb not in table.index:
                continue
            observed, sa, sb = table.loc[source_pair, column], table.loc[ga, column], table.loc[gb, column]
            if pd.isna([observed, sa, sb]).any():
                continue
            if b < a:
                a, b, sa, sb = b, a, sb, sa
            expected = sa + sb
            gi = observed - expected
            records.append({"dataset": cohort, "study": "In4mer2024", "cell_line": cell_line,
                            "gene_a": a, "gene_b": b, "joint_effect": observed,
                            "single_a": sa, "single_b": sb, "expected_additive": expected,
                            "author_gi": gi, "cohen_d": np.nan, "source_pair_id": source_pair,
                            "source_table": "Supplementary table 7", "source_row": table.index.get_loc(source_pair)+2,
                            "study_hit": gi < -1, "hit_definition": "dLFC < -1 (source effect-size threshold; no FDR)",
                            "source_doi": "10.1038/s41467-024-47795-3", "perturbation_type": "enCas12a_KO"})
    return records, singles


def dede_pairs(source, symbol_map, approved):
    records, singles = [], []
    for cell_line in ["A549", "HT29", "OVCAR8"]:
        sheet = f"Dede_{cell_line}"
        table = pd.read_excel(source / "In4mer_Supp_table1.xlsx", sheet_name=sheet, keep_default_na=False)
        cohort = f"dede2020_{cell_line.lower()}"
        for index, row in table.iterrows():
            source_pair = row.iloc[0]
            ga, gb = source_pair.split(":")
            a, b = symbol_map.get(ga, ga), symbol_map.get(gb, gb)
            if a not in approved or b not in approved:
                continue
            sa, sb = row.gene1_mean, row.gene2_mean
            for mapped, original, value in [(a, ga, sa), (b, gb, sb)]:
                singles.append({"dataset": cohort, "study": "Dede2020_reanalysis2024", "cell_line": cell_line,
                                "gene": mapped, "source_symbol": original, "single_effect": value,
                                "source": f"Supplementary table 1 / {sheet}", "effect_units": "log2 fold change"})
            if b < a:
                a, b, sa, sb = b, a, sb, sa
            records.append({"dataset": cohort, "study": "Dede2020_reanalysis2024", "cell_line": cell_line,
                            "gene_a": a, "gene_b": b, "joint_effect": row.obs_mean,
                            "single_a": sa, "single_b": sb, "expected_additive": row.exp_mean,
                            "author_gi": row.deltaMean, "cohen_d": row.CohenD, "source_pair_id": source_pair,
                            "source_table": f"Supplementary table 1 / {sheet}", "source_row": index+2,
                            "study_hit": row.deltaMean < -1 and row.CohenD > 0.8,
                            "hit_definition": "dLFC < -1 and Cohen D > 0.8 (2024 source reanalysis; no FDR)",
                            "source_doi": "10.1186/s13059-020-02173-2", "perturbation_type": "enCas12a_KO"})
    return records, singles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=ROOT / "data/evidence/sources")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data/evidence/external")
    parser.add_argument("--annotation-dir", type=Path, default=ROOT / "data/annotations")
    parser.add_argument("--model-dir", type=Path, default=ROOT / "data/processed")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    symbols = pd.read_csv(args.annotation_dir / "symbol_mapping.csv", keep_default_na=False, dtype=str)
    symbol_map = dict(zip(symbols.source_symbol, symbols.gene))
    approved = set(pd.read_csv(args.annotation_dir / "gene_annotations.csv", keep_default_na=False).gene)
    model_genes = set(pd.read_csv(args.model_dir / "model_genes.csv", keep_default_na=False).gene)
    train = pd.read_csv(args.model_dir / "model_pairs.csv", usecols=["gene_a", "gene_b"],
                        keep_default_na=False, dtype=str)
    seen = set(train.apply(lambda x: "|".join(sorted([x.gene_a, x.gene_b])), axis=1))
    records, singles = in4mer_pairs(args.source_dir, symbol_map, approved)
    more_records, more_singles = dede_pairs(args.source_dir, symbol_map, approved)
    pairs = pd.DataFrame(records + more_records)
    pairs["pair_key"] = pairs.gene_a + "|" + pairs.gene_b
    pairs["sl_score"] = -pairs.author_gi
    pairs["label_direction"] = "negative author_gi means stronger synthetic/synergistic growth defect"
    pairs["label_type"] = "delta_log2_fold_change"
    pairs["pvalue"] = np.nan
    pairs["padj"] = np.nan
    pairs["is_measured_nonhit"] = ~pairs.study_hit
    pairs["both_genes_in_model_domain"] = pairs.gene_a.isin(model_genes) & pairs.gene_b.isin(model_genes)
    pairs["pair_seen_in_any_training_study"] = pairs.pair_key.isin(seen)
    pairs["independent_study"] = True
    pairs["oxaliplatin_experiment"] = False
    pairs.to_csv(args.out_dir / "external_pairs.csv", index=False)
    pairs[pairs.both_genes_in_model_domain].to_csv(args.out_dir / "external_pairs_model_domain.csv", index=False)
    single_table = pd.DataFrame(singles + more_singles)
    single_table = single_table.groupby(["dataset", "study", "cell_line", "gene", "source", "effect_units"],
                                       as_index=False).agg(single_effect=("single_effect", "median"),
                                                           source_symbol=("source_symbol", lambda x: ";".join(sorted(set(x)))))
    single_table.to_csv(args.out_dir / "external_singles.csv", index=False)
    counts = []
    for name, group in pairs.groupby("dataset"):
        covered = group[group.both_genes_in_model_domain]
        unseen = covered[~covered.pair_seen_in_any_training_study]
        counts.append({"dataset": name, "cell_line": group.cell_line.iloc[0],
                       "n_full_measured_pairs": len(group), "n_full_hits": int(group.study_hit.sum()),
                       "n_model_domain_pairs": len(covered), "n_model_domain_hits": int(covered.study_hit.sum()),
                       "n_unseen_model_domain_pairs": len(unseen), "n_unseen_model_domain_hits": int(unseen.study_hit.sum()),
                       "n_full_measured_nonhits": int(group.is_measured_nonhit.sum())})
    summary = pd.DataFrame(counts)
    summary.to_csv(args.out_dir / "coverage_summary.csv", index=False)
    manifest = {
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "experimental_studies": {
            "In4mer2024": {"doi": "10.1038/s41467-024-47795-3", "data_doi": "10.6084/m9.figshare.24243832.v1",
                            "table": "Supp_table7.txt", "download_url": "https://ndownloader.figshare.com/files/44435639",
                            "source_code_url": "https://ndownloader.figshare.com/files/44435618",
                            "gi": "observed pair log2FC minus sum of two measured singleton log2FC",
                            "hit": "dLFC < -1; author-defined effect-size rule; no P value or FDR claim",
                            "platforms": "K562/A549 prototype library; MELJUSO/A375 Inzolia library"},
            "Dede2020_reanalysis2024": {"original_doi": "10.1186/s13059-020-02173-2",
                            "reanalysis_doi": "10.1038/s41467-024-47795-3", "data_doi": "10.6084/m9.figshare.24243832.v1",
                            "table": "Supp_table1.xlsx / Dede_A549, Dede_HT29, Dede_OVCAR8",
                            "download_url": "https://ndownloader.figshare.com/files/44435633",
                            "gi": "published deltaMean retained exactly, not recomputed from rounded singleton columns",
                            "hit": "deltaMean < -1 and CohenD > 0.8; threshold from 2024 reanalysis",
                            "source_rows_per_cell": 403,
                            "mapped_rows_per_cell": 402,
                            "unmapped_source_pair_excluded": "MPP2:MPP6 (MPP6 absent from project's approved HGNC symbol set)"}},
        "preprocessing": "HGNC symbol mapping; canonical unordered pairs; pairs only (triples/quads and ambiguous multigene names excluded); retain only finite pair and both-singleton phenotypes; no selection by hit status",
        "evaluation_only": True, "training_or_tuning_use": False,
        "single_gene_use": "Frozen external scoring may use independently measured source singleton phenotypes; explicitly report this as conditional on observed singles. Do not fit or tune on external pairs.",
        "nonhit_definition": "measured pair not passing its source hit threshold; not a proven biological null; missing pairs are excluded, never encoded as zero/non-hit",
        "applicability": "Dede HT29 is colorectal pairwise evidence under basal growth. Other cohorts assess study/platform/cell-context transfer. No cohort contains oxaliplatin or an acquired-resistant cell state.",
        "overlap_rule": "Independent study can contain training-seen gene pairs. Report all-covered and pair-unseen subsets separately; no assertion of pair-level novelty from study independence alone.",
        "model_domain_size": len(model_genes), "coverage": counts,
        "newer_candidate_sources": {
            "COLO1_2025": {"doi": "10.1038/s41467-025-67256-9", "role": "already a training study; not independent external evaluation"},
            "Functional_modules": {"doi": "10.64898/2025.12.18.695201", "latest_checked_version": "v2 (2026-09-22)",
                       "title": "Functional modules predict genetic interactions, revealing evolutionary conservation and rewiring",
                       "data_repository": "https://github.com/hart-lab/Functional_modules_predict_genetic_interactions",
                       "availability": "Repository and public contents endpoint returned 404. No full measured gene-pair score table obtained; excluded from computational external validation.",
                       "v2_reported_hit_rule": "Zgi FDR < 0.25, normalized effect size > 0.15, and double knockout Zfc < -2",
                       "literature_role": "glycosylation network and 29-gene in vitro/3D/xenograft validation evidence; cannot calculate recall or FDR from a literature-only positive list"}},
    }
    (args.out_dir / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(f"Saved {len(pairs)} measured pairs and {len(single_table)} singleton records")


if __name__ == "__main__":
    main()
