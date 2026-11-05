#!/usr/bin/env python3
"""Prepare CRC oxaliplatin evidence without treating expression as causality.

Counts use PyDESeq2 by default. --de-method exploratory exports normalized
effect sizes only; it never assigns p-values to FPKM, arrays, or the screen
with one PBS control. Cell lines and MD/HD contrasts remain separate.
"""
import argparse
from datetime import datetime, timezone
import gzip
import importlib.metadata
import io
import json
from pathlib import Path
import re
import tarfile
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_URLS = {
    "GSE299428_counts.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE299nnn/GSE299428/suppl/GSE299428_20190218_HCT116_Ox_R_Consolidated_Counts.xlsx",
    "GSE289766_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289766/suppl/GSE289766_RAW.tar",
    "GSE254136_count_table.count.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254136/suppl/GSE254136_count_table.count.txt.gz",
}
META_URLS = {
    "GSE299428_samples.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299428&targ=gsm&form=text&view=brief",
    "GSE42387_samples.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42387&targ=gsm&form=text&view=brief",
    "GPL16297.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16297&targ=self&form=text&view=full",
    "GSE42387_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42387/matrix/GSE42387_series_matrix.txt.gz",
    "GSE254136_metadata.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136&targ=self&form=text&view=brief",
    "GSE254136_samples.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254136&targ=gsm&form=text&view=brief",
}


def soft_samples(path):
    samples = []
    for block in path.read_text().split("^SAMPLE = ")[1:]:
        accession = block.splitlines()[0].strip()
        title = re.search(r"!Sample_title = (.*)", block).group(1).strip()
        descriptions = re.findall(r"!Sample_description = (.*)", block)
        samples.append({"accession": accession, "title": title, "descriptions": descriptions, "block": block})
    return samples


def canonical_matrix(frame, symbol_map, aggregation="sum"):
    frame.index = [symbol_map.get(str(g), str(g)) for g in frame.index]
    return frame.groupby(level=0).agg(aggregation)


def evidence_frame(matrix, controls, treated, study, cell_line, contrast,
                   assay, method, log_scale=False, pseudocount=0.1):
    control_mean = matrix[controls].mean(axis=1)
    treated_mean = matrix[treated].mean(axis=1)
    effect = (treated_mean - control_mean if log_scale else
              np.log2((treated_mean + pseudocount) / (control_mean + pseudocount)))
    return pd.DataFrame({
        "gene": matrix.index, "study": study, "cell_line": cell_line,
        "contrast": contrast, "assay": assay, "method": method,
        "log2fc": effect.values, "control_mean": control_mean.values,
        "treated_mean": treated_mean.values, "n_control": len(controls),
        "n_treated": len(treated), "pvalue": np.nan, "padj": np.nan,
        "adequate_expression": True,
    }).set_index("gene")


def count_contrasts(args, symbols, sample_rows):
    raw = pd.read_excel(args.source_dir / "GSE299428_counts.xlsx").set_index("GENE_ID")
    counts = canonical_matrix(raw, symbols)
    metadata = []
    for sample in soft_samples(args.metadata_dir / "GSE299428_samples.txt"):
        cell_line, group = sample["title"].split(",")[0].split()
        library = next(x.split(": ")[1] for x in sample["descriptions"] if x.startswith("Library name:"))
        record = {"study": "GSE299428", "sample": library, "accession": sample["accession"],
                  "cell_line": cell_line, "condition": group, "title": sample["title"]}
        sample_rows.append(record)
        metadata.append(record)
    meta = pd.DataFrame(metadata).set_index("sample")
    outputs = []
    for cell_line in ("HCT116", "SW480"):
        for condition in ("MD", "HD"):
            submeta = meta[(meta.cell_line == cell_line) & meta.condition.isin(["PAR", condition])]
            selected = counts[submeta.index]
            controls = submeta.index[submeta.condition == "PAR"].tolist()
            treated = submeta.index[submeta.condition == condition].tolist()
            expressed = (selected >= 10).sum(axis=1) >= 3
            if args.de_method == "pydeseq2":
                from pydeseq2.dds import DeseqDataSet
                from pydeseq2.ds import DeseqStats
                from pydeseq2.default_inference import DefaultInference
                inference = DefaultInference(n_cpus=args.n_cpus)
                dds = DeseqDataSet(counts=selected.loc[expressed].T.astype(int),
                                  metadata=submeta[["condition"]], design="~condition",
                                  refit_cooks=True, inference=inference, quiet=True)
                dds.deseq2()
                stats = DeseqStats(dds, contrast=["condition", condition, "PAR"],
                                   inference=inference, quiet=True)
                stats.summary()
                normalized = pd.DataFrame(dds.layers["normed_counts"].T,
                                          index=dds.var_names, columns=dds.obs_names)
                result = evidence_frame(normalized, controls, treated, "GSE299428", cell_line,
                                        f"{condition}_vs_PAR", "RNAseq_counts", "PyDESeq2", pseudocount=1)
                result[["log2fc", "pvalue", "padj"]] = stats.results_df[["log2FoldChange", "pvalue", "padj"]].rename(columns={"log2FoldChange": "log2fc"})
            else:
                positive = selected[(selected > 0).all(axis=1)]
                geometric_mean = np.exp(np.log(positive).mean(axis=1))
                factors = positive.div(geometric_mean, axis=0).median(axis=0)
                factors /= np.exp(np.log(factors).mean())
                normalized = selected.div(factors, axis=1)
                result = evidence_frame(normalized, controls, treated, "GSE299428", cell_line,
                                        f"{condition}_vs_PAR", "RNAseq_counts",
                                        "exploratory_median_ratio_effect", pseudocount=1)
                result["adequate_expression"] = expressed
                result.loc[~expressed, "log2fc"] = np.nan
            outputs.append(result.reset_index())
            print(f"GSE299428 {cell_line} {condition}: {len(result)} genes, {args.de_method}", flush=True)
    return outputs


def fpkm_contrasts(args, symbols, sample_rows):
    series = []
    meta = []
    with tarfile.open(args.source_dir / "GSE289766_RAW.tar") as archive:
        for entry in archive.getmembers():
            accession, cell_line, group, replicate, _ = entry.name.split("_", 4)
            record = {"study": "GSE289766", "sample": entry.name.removesuffix(".txt.gz"),
                      "accession": accession, "cell_line": cell_line, "condition": group, "title": entry.name}
            sample_rows.append(record)
            if group not in ("WT", "OX"):
                continue
            frame = pd.read_csv(io.BytesIO(gzip.decompress(archive.extractfile(entry).read())), sep="\t", index_col=0)
            frame.columns = [record["sample"]]
            series.append(canonical_matrix(frame, symbols))
            meta.append(record)
    matrix = pd.concat(series, axis=1)
    metadata = pd.DataFrame(meta).set_index("sample")
    outputs = []
    for cell_line in ("HCT116", "SW480", "DLD1"):
        selected = metadata[metadata.cell_line == cell_line]
        controls = selected.index[selected.condition == "WT"].tolist()
        treated = selected.index[selected.condition == "OX"].tolist()
        result = evidence_frame(matrix, controls, treated, "GSE289766", cell_line, "OX_vs_WT",
                                "RNAseq_FPKM", "exploratory_log2_ratio_mean_FPKM", pseudocount=0.1)
        result["adequate_expression"] = result[["control_mean", "treated_mean"]].max(axis=1) >= 1
        result.loc[~result.adequate_expression, "log2fc"] = np.nan
        outputs.append(result.reset_index())
    return outputs


def array_contrasts(args, symbols, sample_rows):
    platform_text = (args.metadata_dir / "GPL16297.txt").read_text()
    platform_table = platform_text.split("!platform_table_begin\n")[1].split("!platform_table_end")[0]
    platform = pd.read_csv(io.StringIO(platform_table), sep="\t", dtype=str)
    platform = platform.dropna(subset=["Gene Symbol"])
    # A probe assigned to multiple genes cannot identify one gene unambiguously.
    platform = platform[~platform["Gene Symbol"].str.contains(r"///|;|\|", regex=True)]
    platform["gene"] = platform["Gene Symbol"].map(lambda g: symbols.get(g, g))
    with gzip.open(args.metadata_dir / "GSE42387_series_matrix.txt.gz", "rt") as stream:
        text = stream.read()
    table = text.split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    matrix = pd.read_csv(io.StringIO(table), sep="\t").rename(columns={"ID_REF": "ID"})
    n_input = len(matrix)
    matrix = matrix.merge(platform[["ID", "gene"]], on="ID", how="inner").drop(columns="ID")
    n_mapped = len(matrix)
    matrix = matrix.groupby("gene").median(numeric_only=True)
    metadata = []
    for sample in soft_samples(args.metadata_dir / "GSE42387_samples.txt"):
        base = sample["title"].split("_rep")[0]
        cell_line, condition = base.split("-", 1) if "-" in base else (base, "PAR")
        record = {"study": "GSE42387", "sample": sample["accession"], "accession": sample["accession"],
                  "cell_line": cell_line, "condition": condition, "title": sample["title"],
                  "original_array_column": sample["descriptions"][0]}
        sample_rows.append(record)
        metadata.append(record)
    meta = pd.DataFrame(metadata).set_index("sample")
    outputs = []
    for cell_line in ("HCT116", "HT29", "LoVo"):
        selected = meta[meta.cell_line == cell_line]
        controls = selected.index[selected.condition == "PAR"].tolist()
        treated = selected.index[selected.condition == "OxPt"].tolist()
        result = evidence_frame(matrix, controls, treated, "GSE42387", cell_line, "OX_vs_PAR",
                                "microarray_log2", "exploratory_difference_mean_log2_intensity", log_scale=True)
        outputs.append(result.reset_index())
    return outputs, {"collapsed_rows_input": n_input, "rows_mapped_to_gene": n_mapped, "genes_after_median_collapse": len(matrix)}


def crispr_effects(args, symbols, sample_rows):
    table = pd.read_csv(args.source_dir / "GSE254136_count_table.count.txt.gz", sep="\t")
    table["gene"] = table.Gene.map(lambda g: symbols.get(g, g))
    columns = ["PBS", "Replicate1", "Replicate2", "Replicate3"]
    normalized = table[columns].div(table[columns].sum(axis=0), axis=1) * 1_000_000
    effects = np.log2(normalized[["Replicate1", "Replicate2", "Replicate3"]].add(1).div(normalized.PBS + 1, axis=0))
    effects.columns = ["replicate1_log2fc", "replicate2_log2fc", "replicate3_log2fc"]
    guides = pd.concat([table[["sgRNA", "Gene", "gene", "PBS"]], effects], axis=1)
    guides["mean_log2fc"] = effects.mean(axis=1)
    guides["baseline_usable"] = guides.PBS >= args.min_control_guide_count
    guides["is_hgnc_gene"] = guides.gene.isin(set(symbols.values()))
    guides.to_csv(args.out_dir / "GSE254136_guide_effects.csv", index=False)
    usable = guides[guides.baseline_usable & guides.is_hgnc_gene]
    aggregated = usable.groupby("gene")[["replicate1_log2fc", "replicate2_log2fc", "replicate3_log2fc"]].median()
    aggregated["screen_log2fc"] = usable.groupby("gene").mean_log2fc.median()
    aggregated["n_usable_guides"] = usable.groupby("gene").size()
    aggregated["n_total_guides"] = guides.groupby("gene").size()
    aggregated["guide_depletion_fraction"] = usable.groupby("gene").mean_log2fc.apply(lambda x: (x < 0).mean())
    aggregated["treated_replicate_depletion_fraction"] = (aggregated[["replicate1_log2fc", "replicate2_log2fc", "replicate3_log2fc"]] < 0).mean(axis=1)
    aggregated["screen_pvalue"] = np.nan
    aggregated["screen_padj"] = np.nan
    aggregated["screen_interpretation"] = "single-control exploratory KO drug sensitization; negative means depletion under oxaliplatin"
    aggregated.to_csv(args.out_dir / "GSE254136_gene_effects.csv", index_label="gene")
    for sample in soft_samples(args.metadata_dir / "GSE254136_samples.txt"):
        condition = "PBS" if "control" in sample["title"] else "OX"
        sample_rows.append({"study": "GSE254136", "sample": "PBS" if condition == "PBS" else "Replicate" + sample["title"][-1],
                            "accession": sample["accession"], "cell_line": "HCT116", "condition": condition, "title": sample["title"]})
    return aggregated, {
        "sgRNAs": len(guides),
        "hgnc_target_genes_in_input": int(guides[guides.is_hgnc_gene].gene.nunique()),
        "non_hgnc_control_labels": int(guides[~guides.is_hgnc_gene].gene.nunique()),
        "hgnc_target_genes_with_usable_guides": len(aggregated),
        "guides_with_usable_PBS_count": int(guides.baseline_usable.sum()),
        "usable_gene_targeting_guides": len(usable),
    }


def consolidate(long, screen, threshold):
    long["evidence_key"] = long.study + "_" + long.cell_line + "_" + long.contrast
    wide = long.pivot(index="gene", columns="evidence_key", values="log2fc").add_suffix("_log2fc")
    pvalues = long[long.study == "GSE299428"].pivot(index="gene", columns="evidence_key", values="padj").add_suffix("_padj")
    wide = wide.join(pvalues)
    hct116_columns = ["GSE299428_HCT116_MD_vs_PAR_log2fc", "GSE289766_HCT116_OX_vs_WT_log2fc", "GSE42387_HCT116_OX_vs_PAR_log2fc"]
    hct = wide[hct116_columns]
    # Threshold avoids calling arbitrarily tiny fluctuations a reproducible direction.
    positive = (hct >= threshold).sum(axis=1)
    negative = (hct <= -threshold).sum(axis=1)
    measured = hct.notna().sum(axis=1)
    directional = positive + negative
    wide["hct116_studies_measured"] = measured
    wide["hct116_up_studies"] = positive
    wide["hct116_down_studies"] = negative
    wide["hct116_directional_studies"] = directional
    wide["hct116_direction_concordance"] = pd.concat([positive, negative], axis=1).max(axis=1).div(directional.where(directional >= 2))
    wide["hct116_directional_support_fraction"] = pd.concat([positive, negative], axis=1).max(axis=1).div(measured.where(measured >= 2))
    wide["hct116_all_sign_concordance"] = pd.concat([(hct > 0).sum(axis=1), (hct < 0).sum(axis=1)], axis=1).max(axis=1).div(measured.where(measured >= 2))
    wide["hct116_direction"] = np.select([positive > negative, negative > positive], ["up_in_resistance", "down_in_resistance"], default="mixed_or_small")
    wide["hct116_consistent_three_studies"] = (measured == 3) & ((positive == 3) | (negative == 3))
    wide["expression_evidence_role"] = "resistance-associated expression; not causal target evidence"
    screen = screen.rename(columns={c: "GSE254136_" + c for c in screen.columns})
    return wide.join(screen, how="outer").reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=ROOT / "data/raw")
    parser.add_argument("--metadata-dir", type=Path, default=ROOT / "data/disease/sources")
    parser.add_argument("--annotation-dir", type=Path, default=ROOT / "data/annotations")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "results/disease")
    parser.add_argument("--de-method", choices=["pydeseq2", "exploratory"], default="pydeseq2")
    parser.add_argument("--n-cpus", type=int, default=4)
    parser.add_argument("--direction-threshold", type=float, default=0.5)
    parser.add_argument("--min-control-guide-count", type=int, default=30)
    parser.add_argument("--download", action="store_true", help="Refresh cached public sources")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for directory, urls in ((args.source_dir, RAW_URLS), (args.metadata_dir, META_URLS)):
        directory.mkdir(parents=True, exist_ok=True)
        for name, url in urls.items():
            path = directory / name
            if args.download or not path.exists():
                urllib.request.urlretrieve(url, path)
            records.append({"file": name, "url": url, "bytes": path.stat().st_size,
                            "file_saved_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()})
    mapping = pd.read_csv(args.annotation_dir / "symbol_mapping.csv", keep_default_na=False)
    symbols = dict(zip(mapping.source_symbol, mapping.gene))
    samples = []
    contrasts = count_contrasts(args, symbols, samples)
    contrasts.extend(fpkm_contrasts(args, symbols, samples))
    array, array_mapping = array_contrasts(args, symbols, samples)
    contrasts.extend(array)
    screen, screen_mapping = crispr_effects(args, symbols, samples)
    long = pd.concat(contrasts, ignore_index=True)
    long.to_csv(args.out_dir / "expression_contrasts.csv.gz", index=False)
    evidence = consolidate(long, screen, args.direction_threshold)
    evidence.to_csv(args.out_dir / "disease_evidence.csv", index=False)
    pd.DataFrame(samples).to_csv(args.out_dir / "sample_mapping.csv", index=False)
    manifest = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "sources": records,
        "count_method": args.de_method,
        "pydeseq2_version": importlib.metadata.version("pydeseq2") if args.de_method == "pydeseq2" else None,
        "count_design": "Each cell line and dose analyzed separately: ~condition, MD vs PAR and HD vs PAR, 3 vs 3. Genes with >=10 counts in >=3 samples retained for inferential DE. HD is excluded from primary cross-study direction concordance.",
        "exploratory_count_effect": "Median-of-ratios library normalization from genes positive in all 6 contrast samples; log2((mean_treated+1)/(mean_control+1)); no p-values. Low-expression effects missing using the same >=10 counts in >=3 samples rule.",
        "fpkm_method": "GSE289766: log2((mean OX FPKM+0.1)/(mean WT FPKM+0.1)); effects missing if both means <1 FPKM; 3 vs 3 separately for HCT116, SW480, DLD1; no count model and no p-values.",
        "array_method": "GSE42387 official series matrix contains quantile-normalized log2 signals with systematic identifiers matching GPL16297; unambiguous platform symbols are resolved to HGNC and multiple rows per gene collapsed by median; mean OX minus mean parental, 3 vs 3 separately for HCT116/HT29/LoVo; no additional log transform or p-values.",
        "array_mapping": array_mapping,
        "array_source_choice": "The pre-existing 41000-row probe_level_data supplementary file uses Agilent probe IDs and must NOT be joined directly to the 32705-row collapsed-systematic-ID GPL16297 table. The matching GEO series matrix is used instead.",
        "screen_method": f"GSE254136: total-library CPM per guide/sample, pseudocount 1 CPM; guide log2(treatment/PBS) in each treated replicate, mean across 3 treated replicates, median of guide mean effects per gene. Guides require PBS >= {args.min_control_guide_count} reads. Non-HGNC labels (HGLibA controls) remain in the guide table but are excluded from gene-level evidence. One PBS control prevents estimating control biological variation: no inferential p-values/FDR or causal claim from this effect alone.",
        "screen_direction": "Negative screen log2FC indicates guide depletion under oxaliplatin (loss of gene may sensitize cells); positive indicates enrichment. This is not a direct gene-expression effect.",
        "screen_mapping": screen_mapping,
        "screen_metadata_note": "GEO describes human HCT116 and Homo sapiens but includes Assembly: mm10 in sample processing text; supplied gene symbols/counts are used without genome realignment. This deposited metadata inconsistency is preserved here.",
        "concordance": f"Primary cross-study HCT116 only: GSE299428 MD, GSE289766 OX, GSE42387 OX. Directional study requires abs(log2FC)>={args.direction_threshold}; concordance=max(up,down)/(up+down), missing unless >=2 directional studies. Directional_support_fraction uses all measured studies in the denominator; all_sign_concordance uses all nonzero signs without the effect threshold. Three-study consistency requires all 3 effects exceeding threshold in one direction. Other cell lines and HD are retained in separate columns, never pooled.",
        "interpretation": "Disease expression evidence is an independent prioritization/validation layer, not pairwise synthetic-lethality labels and not proof of causal targets. Single-gene CRISPR evidence also does not establish double-perturbation synergy.",
        "n_genes": len(evidence), "n_expression_contrasts": long[["study", "cell_line", "contrast"]].drop_duplicates().shape[0],
        "n_screen_genes": len(screen), "n_samples": len(samples),
    }
    (args.out_dir / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({key: manifest[key] for key in ("count_method", "n_genes", "n_expression_contrasts", "n_screen_genes", "n_samples", "array_mapping")}, indent=2))


if __name__ == "__main__":
    main()
