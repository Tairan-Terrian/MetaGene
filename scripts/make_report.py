#!/usr/bin/env python3
"""Generate the current Chinese research report and verified scientific figures."""
import argparse
from datetime import datetime
from zoneinfo import ZoneInfo
import json
from pathlib import Path
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from review_delivery import review_evidence

ROOT = Path(__file__).resolve().parents[1]
BLUE = "#0F4D92"
TEAL = "#42949E"
RED = "#B64342"
GRAY = "#767676"
LIGHT = "#CFCECE"
MODEL_NAMES = {
    "linear": "Linear", "naiad": "NAIAD reimpl.", "twohead_gnn": "two-head GNN",
    "twohead_graphless": "Two-head / no graph", "xgb": "XGBoost",
    "target_graphless": "Target / no graph", "target_bioplex": "Target / BioPlex",
    "pooled_graphless": "Pooled / no graph", "pooled_reactome": "Pooled / Reactome",
    "pooled_bioplex": "Pooled / BioPlex", "pooled_dualgraph": "Pooled / dual graph",
    "target_string": "Target / STRING", "pooled_string": "Pooled / STRING", "pooled_multigraph": "Pooled / multi-graph",
}
MODEL_ORDER = list(MODEL_NAMES)
PROTOCOL_NAMES = {
    "pair_holdout": "Pair holdout", "gene_holdout": "Gene holdout",
    "anchor_holdout": "Anchor holdout", "transfer_ht29": "Transfer to HT29",
}
STRATEGY_NAMES = {
    "random": "Random", "effect": "Predicted effect", "uncertainty": "Uncertainty",
    "expected_improvement": "Expected improvement",
}
STRATEGY_COLORS = {"random": GRAY, "effect": RED, "uncertainty": TEAL,
                   "expected_improvement": BLUE}


def read_json(path):
    return json.loads(path.read_text())


def number(value, digits=3):
    return "不可估计" if pd.isna(value) else f"{value:.{digits}f}"


def measurement(value, digits=3):
    return "未测" if pd.isna(value) else f"{value:.{digits}f}"


def mean_sd(values, digits=3):
    x = pd.Series(values).dropna()
    if len(x) == 0:
        return "不可估计"
    if len(x) == 1:
        return f"{x.iloc[0]:.{digits}f}（单次运行）"
    return f"{x.mean():.{digits}f} ± {x.std(ddof=1):.{digits}f}"


def count_range(values):
    x = pd.Series(values).dropna().astype(int)
    low, high = x.min(), x.max()
    return f"{low:,}" if low == high else f"{low:,}–{high:,}"


def table(headers, rows):
    def cell(x):
        return str(x).replace("|", "\\|").replace("\n", "<br>")
    lines = ["| " + " | ".join(map(cell, headers)) + " |",
             "| " + " | ".join(["---"] * len(headers)) + " |"]
    return "\n".join(lines + ["| " + " | ".join(map(cell, r)) + " |" for r in rows])


def ordered_models(models):
    return sorted(models, key=lambda m: (MODEL_ORDER.index(m), m))


def validate_report_inputs(config,metrics,curves,metadata):
    """A final report must not silently use a partial run or stale split design."""
    expected={(protocol,model,seed)
              for protocol,models in config["protocols"].items()
              for model in models
              for seed in config["primary_seeds"] if protocol in {"pair_holdout","gene_holdout"}}
    expected|={(protocol,model,seed)
               for protocol,models in config["protocols"].items()
               for model in models
               for seed in config["secondary_seeds"] if protocol not in {"pair_holdout","gene_holdout"}}
    observed=set(metrics[["protocol","model","seed"]].itertuples(index=False,name=None))
    if expected!=observed:
        raise ValueError(f"Incomplete benchmark registry: missing={sorted(expected-observed)}, extra={sorted(observed-expected)}")
    keys=["protocol","model","seed","dataset","target","subset"]
    if metrics.duplicated(keys).any():
        raise ValueError("Duplicated benchmark metric rows")
    expected_al={(strategy,seed,r)
                 for strategy in config["active_learning"]["strategies"]
                 for seed in config["primary_seeds"]
                 for r in range(config["active_learning"]["query_rounds"]+1)}
    observed_al=set(curves[["strategy","seed","round"]].itertuples(index=False,name=None))
    if expected_al!=observed_al or curves.duplicated(["strategy","seed","round"]).any():
        raise ValueError("Incomplete or duplicated active-learning curves")
    if not curves.labels_observed_plus_validation.eq(curves.n_observed+curves.val_size).all():
        raise ValueError("Active-learning label budget accounting is inconsistent")
    if metadata["prediction_source"]["single_partition_seed"]!=config["partition_seed"]:
        raise ValueError("Candidate predictions do not match the fixed partition design")
    for field,expected_count in [("n_top30",config["candidate_policy"]["shortlist_size"]),
                                 ("n_priority8",config["candidate_policy"]["priority_size"])]:
        if metadata[field]!=expected_count:
            raise ValueError(f"Candidate handoff incomplete: {field}={metadata[field]}")


def set_style():
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next(name for name in ['Arial', 'Helvetica', 'DejaVu Sans'] if name in available)
    plt.rcParams.update({
        "font.family": [family],
        "font.size": 13, "axes.titlesize": 15, "axes.labelsize": 14,
        "axes.linewidth": 1.8, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def export(fig, directory, stem):
    fig.savefig(directory / f"{stem}.png", dpi=300, bbox_inches="tight", pad_inches=0.15)
    fig.savefig(directory / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.15)
    fig.savefig(directory / f"{stem}.svg", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def taxonomy_figure(coverage, out):
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.3))
    keys = ["is_strict_metabolic_enzyme", "is_direct_histone_enzyme", "is_chromatin_regulator", "is_extended_metabolic_signaling"]
    labels = ["Strict metabolic enzymes", "Direct histone enzymes", "All chromatin-associated", "Metabolic signaling"]
    y = np.arange(len(keys))
    for ax, group, title in [(axes[0], "model_gene_counts", "A  Existing model gene domain"),
                             (axes[1], "ht29_single_gene_counts", "B  Genes with HT29 single effects")]:
        vals = [coverage[group][k] for k in keys]
        ax.barh(y, vals, color=[BLUE, TEAL, LIGHT, RED],
                edgecolor="#272727", lw=0.7, height=0.65)
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.set_xlim(0, max(vals)*1.18)
        ax.set_xlabel("Genes (categories can overlap)")
        ax.set_title(title, pad=12)
        for yi, value in zip(y, vals):
            ax.text(value+max(vals)*0.025, yi, str(value), va="center")
    fig.tight_layout(pad=1.5, w_pad=2)
    export(fig, out, "01_annotation_scope")


def benchmark_figures(metrics, selected_model, out):
    for target, title, stem in [("additive_gi_proxy", "GI proxy", "02_benchmark_gi"),
                                ("joint", "Combined growth phenotype", "03_benchmark_joint")]:
        data = metrics[(metrics.protocol=="pair_holdout") & (metrics.subset=="all") &
                       (metrics.target==target) & (metrics.dataset=="ht29_colo1")]
        models = ordered_models(data.model.unique())
        fig, axes = plt.subplots(1, 2, figsize=(13.6, max(5.5, len(models)*0.48+1.4)), sharey=True)
        for ax, metric, xlabel in zip(axes, ["spearman", "rmse"], ["Spearman correlation", "RMSE (HT29 source scale)"]):
            for i, model in enumerate(models):
                values = data.loc[data.model==model, metric].dropna().to_numpy()
                color = BLUE if model==selected_model else RED if model.startswith("twohead_") else GRAY
                jitter = np.linspace(-0.12, 0.12, len(values))
                ax.scatter(values, i+jitter, s=23, color=color, alpha=0.55, zorder=3)
                ax.errorbar(np.mean(values), i, xerr=np.std(values, ddof=1) if len(values)>1 else 0,
                            fmt="D", ms=6, color=color, lw=2, capsize=4, zorder=4)
            ax.set_xlabel(xlabel)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax.set_yticks(np.arange(len(models)), [MODEL_NAMES[m] for m in models])
            ax.grid(axis="x", color="#E5E5E5", lw=0.7, zorder=0)
        axes[0].invert_yaxis()
        axes[0].set_title("A  Rank prediction")
        axes[1].set_title("B  Error on held-out pairs")
        fig.suptitle(f"{title} · fixed pair split · points: training seeds; mean ± training-seed SD", fontsize=14, y=1.01)
        fig.tight_layout(pad=1.5)
        export(fig, out, stem)

    data = metrics[(metrics.target=="additive_gi_proxy") & (metrics.subset=="all") &
                   (metrics.dataset=="ht29_colo1")]
    matrix = data.pivot_table(index="protocol", columns="model", values="spearman", aggfunc="mean")
    matrix = matrix.reindex(index=list(PROTOCOL_NAMES), columns=ordered_models(matrix.columns))
    cmap = LinearSegmentedColormap.from_list("rho", [RED, "white", BLUE])
    cmap.set_bad("#ECECEC")
    fig, ax = plt.subplots(figsize=(max(12, matrix.shape[1]*1.18), 4.9))
    im = ax.imshow(np.ma.masked_invalid(matrix.to_numpy()), cmap=cmap, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(matrix.shape[1]), [MODEL_NAMES[m] for m in matrix.columns], rotation=35, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]), [PROTOCOL_NAMES[p] for p in matrix.index])
    ax.tick_params(length=0, pad=8)
    ax.spines[["bottom", "left"]].set_visible(False)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i,j]
            ax.text(j, i, "—" if pd.isna(value) else f"{value:.3f}", ha="center", va="center",
                    fontsize=11, color="white" if abs(value)>0.6 else "#272727")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.03, label="Mean GI-proxy Spearman")
    ax.set_title("Generalization protocols · gray: unrun or undefined", pad=16)
    fig.tight_layout(pad=1.5)
    export(fig, out, "04_generalization")


def uncertainty_figure(metrics, selected_model, calibration, out):
    data = metrics[(metrics.model==selected_model) & (metrics.target=="additive_gi_proxy") &
                   (metrics.subset=="all") & (metrics.dataset=="ht29_colo1")].copy()
    data["protocol_order"] = data.protocol.map({p:i for i,p in enumerate(PROTOCOL_NAMES)})
    data = data.sort_values(["protocol_order", "seed"])
    fig, axes = plt.subplots(1, 2, figsize=(12.4, max(5.0, len(data)*0.30+1.5)))
    labels = [f"{PROTOCOL_NAMES[r.protocol]} / {r.seed}" for r in data.itertuples()]
    for ax, metric, xlabel in zip(axes, ["spearman", "rmse"], ["Spearman with 95% dyadic bootstrap CI", "RMSE with 95% dyadic bootstrap CI"]):
        for i, row in enumerate(data.itertuples()):
            estimate = getattr(row, metric)
            lo, hi = getattr(row, metric+"_ci_low"), getattr(row, metric+"_ci_high")
            if np.isfinite([estimate, lo, hi]).all():
                ax.plot([lo, hi], [i,i], color=BLUE, lw=1.8)
                ax.scatter(estimate, i, color=BLUE, s=24, zorder=3)
            else:
                ax.text(.5,i,"not estimable",transform=ax.get_yaxis_transform(),
                        ha="center",va="center",color=GRAY,fontsize=9)
        ax.set_yticks(np.arange(len(data)), labels if ax is axes[0] else [])
        ax.set_ylim(len(data)-.3,-.7)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.grid(axis="x", color="#E5E5E5", lw=0.7)
    fig.suptitle(f"Validation-selected model: {MODEL_NAMES[selected_model]} · each row is one fit", y=1.01, fontsize=14)
    fig.tight_layout(pad=1.5)
    export(fig, out, "05_metric_uncertainty")

    selected = calibration[calibration.model==selected_model]
    protocols = [p for p in PROTOCOL_NAMES if p in set(selected.protocol)]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.7))
    for i, protocol in enumerate(protocols):
        group = selected[selected.protocol==protocol]
        x = np.full(len(group), i)+np.linspace(-0.08, 0.08, len(group))
        axes[0].scatter(x, group.coverage, s=36, color=BLUE)
        axes[1].scatter(x, group.mean_width, s=36, color=TEAL)
    for ax in axes:
        ax.set_xticks(np.arange(len(protocols)), [PROTOCOL_NAMES[p] for p in protocols], rotation=25, ha="right")
        ax.grid(axis="y", color="#E5E5E5", lw=0.7)
    axes[0].axhline(selected.nominal_coverage.iloc[0], color=RED, ls="--", lw=1.7)
    axes[0].set_ylim(0,1.02)
    axes[0].set_ylabel("Held-out empirical interval coverage")
    axes[0].set_title("A  Calibration-set intervals, evaluated on test")
    axes[1].set_ylabel("Mean interval width (HT29 GI units)")
    axes[1].set_title("B  Interval width")
    fig.tight_layout(pad=1.5)
    export(fig, out, "06_prediction_intervals")


def learning_figure(curves, out):
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.9))
    for strategy in STRATEGY_NAMES:
        group = curves[curves.strategy==strategy]
        stats = group.groupby("round").agg(n=("n_observed", "mean"),
                recall=("recall_top100", "mean"), recall_sd=("recall_top100", "std"),
                rho=("holdout_spearman", "mean"), rho_sd=("holdout_spearman", "std")).sort_index()
        for ax, metric, spread in zip(axes, ["recall", "rho"], ["recall_sd", "rho_sd"]):
            ax.plot(stats.n, stats[metric], marker="o", ms=4, lw=2.2,
                    color=STRATEGY_COLORS[strategy], label=STRATEGY_NAMES[strategy])
            ax.fill_between(stats.n, stats[metric]-stats[spread], stats[metric]+stats[spread],
                            color=STRATEGY_COLORS[strategy], alpha=0.13, lw=0)
            ax.set_xlabel("Observed pool labels\n(fixed validation-label cost excluded)")
            ax.grid(axis="y", color="#E5E5E5", lw=0.7)
    axes[0].set_ylim(0,1)
    axes[0].set_ylabel("Recall of lowest 100 GI proxies\n(measured pool)")
    axes[0].set_title("A  Discovery objective")
    axes[1].set_ylabel("Fixed test-set GI Spearman")
    axes[1].set_title("B  Prediction generalization")
    axes[1].legend(loc="best", fontsize=10)
    fig.tight_layout(pad=1.6)
    export(fig, out, "07_active_learning")


def learning_subset_figure(curves,out):
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.1))
    for ax,prefix,title in [(axes[0],"strict","A  Enzyme × broad chromatin"),
                             (axes[1],"core","B  Enzyme × direct histone enzyme")]:
        if int(curves[prefix+"_hit_k"].iloc[0]) == 0:
            ax.text(.5,.5,"Not estimable: measured core pool < 20 pairs",
                    transform=ax.transAxes,ha="center",va="center",fontsize=11,color=GRAY)
        for strategy in STRATEGY_NAMES:
            group=curves[curves.strategy.eq(strategy)]
            stats=group.groupby("round").agg(n=("n_observed","mean"),
                recall=(prefix+"_recall_topk","mean"),sd=(prefix+"_recall_topk","std")).sort_index()
            ax.plot(stats.n,stats.recall,color=STRATEGY_COLORS[strategy],lw=2,marker="o",ms=3,
                    label=STRATEGY_NAMES[strategy])
            ax.fill_between(stats.n,stats.recall-stats.sd,stats.recall+stats.sd,
                            color=STRATEGY_COLORS[strategy],alpha=.12,lw=0)
        target_k = int(curves[prefix+'_hit_k'].iloc[0])
        pool_n = int(curves[prefix+'_pool_n'].iloc[0])
        detail = f'{target_k} lowest-GI targets in {pool_n} measured pairs' if target_k else f'{pool_n} measured pairs; target set not defined'
        ax.set_title(title+'\n'+detail,pad=12)
        ax.set_ylim(0,1)
        budget_low, budget_high = curves.n_observed.min(), curves.n_observed.max()
        budget_margin = (budget_high - budget_low) * 0.04
        ax.set_xlim(budget_low - budget_margin, budget_high + budget_margin)
        ax.set_xlabel("Observed pool labels\n(fixed validation-label cost excluded)")
        ax.set_ylabel("Recall of lowest-GI targets\n(within subset)")
        ax.grid(axis="y",color="#E5E5E5",lw=.7)
    axes[1].legend(loc="upper left",fontsize=9)
    fig.tight_layout(pad=1.6)
    export(fig,out,"13_active_learning_subsets")


def external_coverage_figure(coverage, out):
    data = coverage.sort_values(["cell_line", "dataset"]).copy()
    data["name"] = np.where(data.dataset.str.startswith("dede"), "Dede / ", "In4mer / ") + data.cell_line
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), sharey=True)
    y = np.arange(len(data))
    axes[0].barh(y, data.n_full_measured_pairs, color=LIGHT, label="Full mapped measured cohort", height=0.7)
    axes[0].barh(y, data.n_model_domain_pairs, color=BLUE, label="Both genes in model domain", height=0.7)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Measured gene pairs (log scale)")
    axes[0].set_yticks(y, data.name)
    axes[0].invert_yaxis()
    axes[0].set_title("A  Model coverage")
    axes[0].legend(loc="lower right", fontsize=9)
    axes[1].barh(y-0.17, data.n_model_domain_pairs, color=BLUE, height=0.3, label="All covered")
    axes[1].barh(y+0.17, data.n_unseen_model_domain_pairs, color=TEAL, height=0.3, label="Pair-unseen")
    for i, row in enumerate(data.itertuples()):
        axes[1].text(row.n_model_domain_pairs+2, i-0.17, f"{row.n_model_domain_pairs} / {row.n_model_domain_hits} hits", va="center", fontsize=10)
        noun="hit" if row.n_unseen_model_domain_hits==1 else "hits"
        axes[1].text(row.n_unseen_model_domain_pairs+2, i+0.17, f"{row.n_unseen_model_domain_pairs} / {row.n_unseen_model_domain_hits} {noun}", va="center", fontsize=10)
    axes[1].set_xlim(0,data.n_model_domain_pairs.max()*1.6)
    axes[1].set_xlabel("Evaluable pairs and source-defined hits")
    axes[1].set_title("B  Pair identity overlap separated")
    axes[1].legend(loc="lower right", fontsize=9)
    fig.tight_layout(pad=1.5, w_pad=2)
    export(fig, out, "08_external_coverage")


def external_performance_figure(metrics, selected_model, out):
    data = metrics[(metrics.model==selected_model) & (metrics.subset=="all_covered")]
    datasets = list(data.dataset.unique())
    fig, axes = plt.subplots(1,2,figsize=(12.3,max(4.5,len(datasets)*0.43+1.6)),sharey=True)
    for i,dataset in enumerate(datasets):
        block = data[data.dataset==dataset]
        for ax,metric in zip(axes,["spearman","auprc"]):
            values=block[metric].dropna()
            ax.scatter(values, np.full(len(values),i)+np.linspace(-0.10,0.10,len(values)),color=BLUE,s=27)
            ax.scatter(values.mean(),i,color=BLUE,marker="D",s=44)
        axes[1].scatter(block.baseline_prevalence.mean(),i,color=GRAY,marker="|",s=170,lw=2)
    for ax in axes:
        ax.set_yticks(np.arange(len(datasets)),datasets)
        ax.grid(axis="x",color="#E5E5E5",lw=0.7)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("GI rank Spearman")
    axes[0].set_title("A  Independent-study ranking")
    axes[1].set_xlabel("AUPRC (gray tick: measured-hit prevalence)")
    axes[1].set_xlim(0,1)
    axes[1].set_title("B  Source-defined hit enrichment")
    fig.tight_layout(pad=1.5)
    export(fig,out,"09_external_performance")


def resistance_figure(effects, camera, out):
    primary=effects[effects.contrast=="resistant_vs_parental_untreated"]
    fig,axes=plt.subplots(1,2,figsize=(14.2,6.6))
    neglog=-np.log10(primary.padj_bh.clip(lower=np.nextafter(0,1)))
    significant=(primary.padj_bh<0.05)&(primary.log2fc.abs()>=1)
    colors=np.where(significant,np.where(primary.log2fc>0,RED,BLUE),LIGHT)
    axes[0].scatter(primary.log2fc,neglog,c=colors,s=9,alpha=0.65,lw=0,rasterized=True)
    axes[0].axhline(-np.log10(.05),c=GRAY,ls="--",lw=1)
    for x in [-1,1]:axes[0].axvline(x,c=GRAY,ls=":",lw=1)
    marker_offsets={"CXCL1":(24,-24),"CXCL8":(23,25),"DUSP1":(-78,25)}
    for row in primary[primary.gene.isin(marker_offsets)].itertuples():
        axes[0].annotate(row.gene,(row.log2fc,-np.log10(max(row.padj_bh,np.nextafter(0,1)))),
                         xytext=marker_offsets[row.gene],textcoords="offset points",fontsize=10,
                         arrowprops={"arrowstyle":"-","color":GRAY,"lw":.8},
                         bbox={"facecolor":"white","edgecolor":"none","pad":1.5})
    axes[0].set_xlabel("Log2 fold change: untreated resistant − parental")
    axes[0].set_ylabel("−log10 BH FDR")
    axes[0].set_title("A  GSE76092, independent 3 versus 3")
    paths=camera[camera.contrast=="resistant_vs_parental_untreated"].nsmallest(10,"FDR").iloc[::-1]
    y=np.arange(len(paths));val=-np.log10(paths.FDR.clip(lower=np.nextafter(0,1)))
    axes[1].barh(y,val,color=np.where(paths.Direction=="Up",RED,BLUE),edgecolor="#272727",lw=.5)
    axes[1].set_yticks(y,[textwrap.fill(name,width=38) for name in paths.pathway_name],fontsize=10)
    axes[1].axvline(-np.log10(.05),c=GRAY,ls="--",lw=1)
    axes[1].set_xlabel("−log10 CAMERA pathway FDR")
    axes[1].set_title("B  Top 10 pathways by FDR (red: up, blue: down)")
    fig.tight_layout(pad=1.6,w_pad=2)
    export(fig,out,"10_resistance_expression")


def candidate_figures(top, priority, out, ht29_proxy_median):
    data = top.sort_values("top30_rank").copy()
    y = np.arange(len(data))
    labels = [f"{int(r.top30_rank):02d}  {r.metabolic_gene} + {r.chromatin_gene}"
              for r in data.itertuples()]
    fig, axes = plt.subplots(1, 2, figsize=(14, max(8, len(data)*.34+2)), sharey=True,
                             gridspec_kw={"width_ratios":[1.2, 1]})
    for i, row in enumerate(data.itertuples()):
        color = BLUE if pd.notna(row.priority8_rank) else GRAY
        axes[0].plot([row.pi90_low,row.pi90_high],[i,i], color=color, lw=1.3, alpha=.8)
        axes[0].scatter(row.pred_gi,i,color=color,s=31,
                        marker="D" if pd.notna(row.priority8_rank) else "o",zorder=3)
        m, c = ((row.gene_a_colo1_single_control_scaled,row.gene_b_colo1_single_control_scaled)
                if row.metabolic_gene==row.gene_a else
                (row.gene_b_colo1_single_control_scaled,row.gene_a_colo1_single_control_scaled))
        axes[1].plot([m,c],[i,i],color=LIGHT,lw=1.5)
        axes[1].scatter(m,i,color=TEAL,marker="o",s=29,zorder=3)
        axes[1].scatter(c,i,color=RED,marker="s",s=29,zorder=3)
        absent=[name for name,value in [("M",m),("C",c)] if not np.isfinite(value)]
        if absent:axes[1].text(.98,i,"/".join(absent)+" not measured",transform=axes[1].get_yaxis_transform(),
                               ha="right",va="center",fontsize=8,color=GRAY)
    axes[0].set_yticks(y,labels,fontsize=10)
    axes[0].invert_yaxis()
    axes[0].axvline(0,color=RED,ls="--",lw=1)
    axes[0].axvline(ht29_proxy_median,color=GRAY,ls=":",lw=1.3)
    axes[0].set_xlabel("Predicted GI proxy with marginal 90% interval\n(dotted: measured HT29 median; dashed: zero)")
    axes[0].set_title("A  Diversified Top 30 (blue diamonds: Priority 8)",pad=14)
    axes[1].axvline(data.strong_single_cutoff.iloc[0],color=GRAY,ls="--",lw=1.3)
    axes[1].set_xlabel("COLO1 control-scaled singleton effect\n(negative control 0; positive control −1)")
    axes[1].set_title(f"B  Control-scaled single effects\nCircle: metabolic; square: histone enzyme; dashed: {data.strong_single_cutoff.iloc[0]:g}",pad=14)
    for ax in axes:ax.grid(axis="x",color="#E5E5E5",lw=.7)
    fig.tight_layout(pad=1.5,w_pad=2)
    export(fig,out,"11_candidate_predictions")

    data = priority.sort_values("priority8_rank")
    labels = [f"{int(r.priority8_rank)}  {r.metabolic_gene} + {r.chromatin_gene}"
              for r in data.itertuples()]
    resistance, drug, fdr = [], [], []
    for row in data.to_dict("records"):
        arms = ["a","b"] if row["metabolic_gene"]==row["gene_a"] else ["b","a"]
        resistance.append([row[f"gene_{arm}_GSE76092_log2fc"] for arm in arms])
        fdr.append([row[f"gene_{arm}_GSE76092_padj_bh"] for arm in arms])
        drug.append([row[f"gene_{arm}_GSE254136_screen_log2fc"] for arm in arms])
    cmap=LinearSegmentedColormap.from_list("candidate_evidence",[BLUE,"white",RED])
    cmap.set_bad("#ECECEC")
    fig,axes=plt.subplots(1,2,figsize=(12.7,max(5.6,len(data)*.48+2.4)),sharey=True)
    for panel,ax,values,title in zip(range(2),axes,[resistance,drug],
            ["A  Resistance-state expression log2FC", "B  Single-gene drug-screen effect"]):
        values=np.asarray(values,dtype=float)
        finite=values[np.isfinite(values)]
        limit=max(.1,float(np.max(np.abs(finite)))) if len(finite) else 1
        im=ax.imshow(np.ma.masked_invalid(values),cmap=cmap,vmin=-limit,vmax=limit,aspect="auto")
        ax.set_xticks([0,1],["Metabolic enzyme","Histone enzyme"])
        ax.set_yticks(np.arange(len(data)),labels,fontsize=10)
        ax.tick_params(length=0,pad=8)
        ax.spines[["bottom","left"]].set_visible(False)
        for i in range(len(data)):
            for j in range(2):
                value=values[i,j]
                mark="*" if panel==0 and np.isfinite(fdr[i][j]) and fdr[i][j]<.05 else ""
                label="not measured" if not np.isfinite(value) else f"{value:+.2f}{mark}"
                ax.text(j,i,label,ha="center",va="center",fontsize=10,
                        color="white" if np.isfinite(value) and abs(value)>.65*limit else "#272727")
        fig.colorbar(im,ax=ax,fraction=.045,pad=.04,label="Log2 effect")
        ax.set_title(title,pad=15)
    axes[0].set_xlabel("HTOXAR3 − HT29, untreated; * BH FDR < 0.05",fontsize=11,labelpad=12)
    axes[1].set_xlabel("GSE254136; depletion under OXA; no inferential FDR",fontsize=11,labelpad=12)
    fig.tight_layout(pad=1.7,w_pad=2)
    export(fig,out,"12_priority_evidence")

def benchmark_table(metrics, target='additive_gi_proxy', subset='all', protocol=None):
    data = metrics[(metrics.target == target) & (metrics.subset == subset)]
    if protocol:
        data = data[data.protocol.eq(protocol)]
    rows = []
    for (p, model), group in data.groupby(['protocol', 'model'], sort=False):
        rows.append([p, MODEL_NAMES.get(model, model), count_range(group.n),
                     mean_sd(group.spearman), mean_sd(group.pearson), mean_sd(group.rmse),
                     mean_sd(group.top_k_overlap)])
    return table(['协议', '模型', '测试n', 'Spearman', 'Pearson', 'RMSE', 'Top-k重叠'], rows)


def input_coverage_figure(root, manifest, out):
    data = pd.read_csv(root / 'data/evidence/integrated_gene_evidence.csv')
    mapping = {
        'DepMap': 'DepMap_CRC_median', 'TCGA-COAD': 'TCGA_tumor_median_log2_RSEM',
        'CRC single-cell': 'sc_tumor_epithelial_mean_log1p_CP10K',
        'HCT116 model flux': 'flux_HCT116_Control_mean_log1p_abs_bound',
        'Histone promoter map': 'ChIP_H3K4me3_HCT116_log1p_peak_count',
    }
    names, counts = list(mapping), [int(data[c].notna().sum()) for c in mapping.values()]
    names += ['STRING connected']; counts += [manifest['string']['covered_genes']]
    fig, ax = plt.subplots(figsize=(10, 5))
    y = np.arange(len(names))
    ax.barh(y, counts, color=[BLUE, TEAL, BLUE, TEAL, BLUE, TEAL], edgecolor='#272727', lw=.6)
    ax.set_yticks(y, names); ax.invert_yaxis()
    ax.set_xlim(0, len(data) * 1.12)
    ax.axvline(len(data), color=GRAY, ls=':', lw=1.5, zorder=0)
    for i, n in enumerate(counts):
        ax.text(n + 12, i, str(n), va='center', zorder=3,
                bbox={'facecolor':'white','edgecolor':'none','pad':0.4})
    ax.set_xlabel('Model genes with mapped observations (dotted: full domain)')
    ax.set_title('Independent input coverage; absence remains explicitly missing', pad=14)
    fig.tight_layout(pad=1.7)
    export(fig, out, '14_input_coverage')


def requirement_rows(root, inp, metadata):
    return [
        ['核心对象', '计算已完成', f"38个代谢酶×15个直接组蛋白修饰酶；{metadata['n_strict_core_unmeasured']}对当前HT29未测", 'data/evidence/direct_histone_enzyme_audit.csv', '其他表观调控因子单独导出'],
        ['DepMap单敲', '已接入模型', f"{inp['depmap']['total_screened_models']}个模型，CRC {inp['depmap']['CRC_screened_models']}个", 'data/evidence/depmap_gene_evidence.csv', '亲本细胞单敲；本快照SW480无可用值，不是耐药组合标签'],
        ['TCGA', '已接入模型', f"COAD肿瘤{inp['tcga']['counts']['tumor']['patients']}例、正常{inp['tcga']['counts']['normal']['patients']}例", 'data/evidence/tcga_gene_evidence.csv', '独立组织表达背景'],
        ['单细胞', '已接入模型', f"GSE132465：{inp['single_cell']['n_patients']}患者、{inp['single_cell']['n_cells']}细胞", 'data/evidence/single_cell_gene_evidence.csv', '按患者等权汇总；沿用原作者细胞标签'],
        ['STRING网络', '已接入图模型', f"版本{inp['string']['version']}，{inp['string']['undirected_edges']}边", 'data/evidence/string_edges.csv', '物理网络评分≥0.7；不代表耐药特异直接结合'],
        ['代谢流', '已接入模型估计', f"HCT116/OxR作者模型{inp['metabolic_flux']['reaction_count']}反应范围，映射{inp['metabolic_flux']['covered_genes']}个建模基因", 'data/evidence/flux_gene_evidence.csv', '按生长率归一化；不是实测细胞内通量或KO效应'],
        ['表观调控图谱', '已接入描述性图谱', f"GSE138654三种修饰×两状态；{inp['epigenetic_map']['promoter_mapped_genes']}基因可映射", 'data/evidence/histone_gene_evidence.csv', '无重复差异结合检验；不能指定某修饰酶的直接靶点'],
        ['GNN预测', '已重跑', '验证集选择GNN架构并单列最佳基线', 'results/benchmarks/selection.json', 'GI代理量，不是校准的致死概率'],
        ['Top500', '计算已完成', f"{metadata['n_top500']}对严格核心组合", 'results/candidates/top500.csv', 'HT29原筛选未测；部分无原筛选single；跨背景外推'],
        ['HCT116-R/SW480-R背景', '公开关联已整理；合作方模型待确认', '分别保留HCT116和SW480耐药表达；HCT116/OxR通量及ChIP', 'results/drug_evidence/disease_evidence.csv', '不同研究的耐药株不自动等同'],
        ['聚焦文库衔接', '靶点与容量草案已完成；guide设计未完成', '500对基因靶点、1000条容量预留', 'results/handoff/library_capacity_plan.csv', '不含已验证guide序列，不是可提交合成的订单'],
        ['2～3对代表候选', '计算提名已完成', f"{metadata['n_representative3']}对且基因不重复", 'results/handoff/first_round_representatives.csv', '等待医学院实际验证'],
        ['主动学习', '回顾性计算已完成；前瞻闭环待反馈', '4策略×5种子×12查询轮；真实标签回传模板', 'results/active_learning/curves.csv', '核心实测池不足时不报告命中收益'],
        ['协同致死概率与耐药逆转', '尚无直接证据', '保留可证伪的候选假说', 'results/candidates/candidate_evidence_cards.md', '模型负GI、区间或表达关联均不能替代新实验'],
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    args = parser.parse_args()
    root, results = args.root, args.root / 'results'
    figures = results / 'figures'; figures.mkdir(exist_ok=True)
    config = read_json(root / 'configs/experiment.json')
    annotation = read_json(root / 'data/evidence/annotation_coverage.json')
    inp = read_json(root / 'data/evidence/integrated_input_manifest.json')
    metrics = pd.read_csv(results / 'benchmarks/metrics.csv')
    subset_counts = pd.read_csv(results / 'benchmarks/subset_counts.csv')
    selection = read_json(results / 'benchmarks/selection.json')
    selected = selection['selected_model']; best = selection['best_overall_model']
    calibration = pd.read_csv(results / 'benchmarks/calibration.csv').rename(columns={
        'interval_n': 'n', 'interval_coverage': 'coverage', 'interval_mean_width': 'mean_width'})
    curves = pd.read_csv(results / 'active_learning/curves.csv')
    external = pd.read_csv(results / 'external/metrics.csv')
    ext_coverage = pd.read_csv(root / 'data/evidence/external/coverage_summary.csv')
    effects = pd.read_csv(results / 'disease/GSE76092_gene_effects.csv')
    contrasts = pd.read_csv(results / 'disease/contrast_summary.csv')
    camera = pd.read_csv(results / 'disease/reactome_camera.csv')
    core = pd.read_csv(results / 'candidates/core_candidates.csv', low_memory=False)
    top500 = pd.read_csv(results / 'candidates/top500.csv')
    top = pd.read_csv(results / 'candidates/top30.csv', low_memory=False)
    priority = pd.read_csv(results / 'candidates/priority8.csv', low_memory=False)
    reps = pd.read_csv(results / 'candidates/representative3.csv', low_memory=False)
    metadata = read_json(results / 'candidates/ranking_metadata.json')
    diagnostics = read_json(results / 'diagnostics.json')
    runtime = read_json(results / 'runtime.json')
    audit = read_json(results / 'integrity_audit.json')
    validate_report_inputs(config, metrics, curves, metadata)
    if not audit['passed'] or 'completed_utc' not in runtime:
        raise ValueError('Refuse to report an incomplete or failed run as finished')
    evidence_review = review_evidence(root)
    requirements = requirement_rows(root, inp, metadata)
    pd.DataFrame(requirements, columns=['requirement', 'status', 'actual_delivery', 'evidence_file', 'boundary']).to_csv(results / 'data_requirements_audit.csv', index=False)
    set_style()
    taxonomy_figure(annotation, figures)
    benchmark_figures(metrics, selected, figures)
    uncertainty_figure(metrics, selected, calibration, figures)
    learning_figure(curves, figures)
    learning_subset_figure(curves, figures)
    external_coverage_figure(ext_coverage, figures)
    external_performance_figure(external, selected, figures)
    resistance_figure(effects, camera, figures)
    candidate_figures(top, priority, figures, metadata['source_proxy_distribution']['median'])
    input_coverage_figure(root, inp, figures)
    def img(stem, caption):
        return f'![{caption}](figures/{stem}.png)\n\n{caption}'
    def choose(model, protocol='pair_holdout', subset='all'):
        return metrics[metrics.model.eq(model) & metrics.protocol.eq(protocol) & metrics.subset.eq(subset) & metrics.target.eq('additive_gi_proxy')]
    selected_primary = choose(selected)
    final = curves[curves['round'].eq(config['active_learning']['query_rounds'])]
    primary = effects[effects.contrast.eq('resistant_vs_parental_untreated')]
    full_pairs = pd.read_csv(root / 'data/processed/model_pairs.csv')
    proxy = metadata['source_proxy_distribution']
    graph_summary = []
    for name, edges in dict(np.load(root / 'data/evidence/graphs.npz')).items():
        graph_summary.append([name, len({tuple(sorted(pair)) for pair in edges.T.tolist()}), len(np.unique(edges))])
    partition_rows = []
    for protocol in config['protocols']:
        example = next((results / 'benchmarks/runs').glob(f'{protocol}_*/partition_indices.npz'))
        part = dict(np.load(example))
        partition_rows.append([protocol, len(part['train']), len(part['val']), len(part['calibration']), len(part['test']), len(part['pooled_train']), len(part['discarded_ht29'])])
    performance_rows = []
    for protocol in config['protocols']:
        group = choose(selected, protocol)
        performance_rows.append([protocol, count_range(group.n) if len(group) else '此架构未注册该协议',
                                 mean_sd(group.spearman), mean_sd(group.rmse)])
    al_rows = []
    for strategy in STRATEGY_NAMES:
        group = final[final.strategy.eq(strategy)]
        al_rows.append([STRATEGY_NAMES[strategy], count_range(group.n_observed), count_range(group.labels_observed_plus_validation),
                       mean_sd(group.recall_top100 * 100, 1), mean_sd(group.holdout_spearman),
                       mean_sd(group.core_recall_topk) if int(group.core_hit_k.iloc[0]) else '核心池不足20对，未计算'])
    ext_rows = []
    for (dataset, subset), group in external[external.model.eq(selected)].groupby(['dataset', 'subset']):
        ext_rows.append([dataset, subset, count_range(group.n), count_range(group.n_hits),
                         mean_sd(group.spearman), mean_sd(group.auprc), number(group.baseline_prevalence.iloc[0])])
    candidate_rows = []
    for r in reps.itertuples():
        candidate_rows.append([r.representative_rank, f'{r.metabolic_gene} + {r.chromatin_gene}', number(r.pred_gi),
            f'[{number(r.pi90_low)}, {number(r.pi90_high)}]', number(r.pred_joint), number(r.source_relative_GI_percentile, 1),
            int(r.strong_single_count) if pd.notna(r.strong_single_count) else '未知',
            f'{int(r.target_resistance_association_count)}/2背景达到阈值', '探索性提名；耐药组合未验证'])
    top_rows = [[int(r.top30_rank), f'{r.metabolic_gene} + {r.chromatin_gene}', number(r.pred_gi), number(r.pred_joint),
                int(r.strong_single_count) if pd.notna(r.strong_single_count) else '未知',
                int(r.priority8_rank) if pd.notna(r.priority8_rank) else '—'] for r in top.itertuples()]
    n_runs = metrics[['protocol', 'model', 'seed']].drop_duplicates().shape[0]
    main_core_test = subset_counts[subset_counts.protocol.eq('pair_holdout') & subset_counts.subset.eq('core_metabolism_histone')].n_test.iloc[0]
    text = [
        '# MetaGene实验汇报：AI主动学习驱动的代谢—表观遗传跨维度合成致死候选筛选与奥沙利铂耐药逆转研究',
        f"报告更新时间：{datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='minutes')}；数值实验完成时间：{datetime.fromisoformat(runtime['completed_utc']).astimezone(ZoneInfo('Asia/Shanghai')).isoformat(timespec='minutes')}。本报告引用已完成的**{n_runs}个基准拟合、20条12轮主动学习序列、5个最终候选模型**。运行设备：{', '.join(runtime['gpus'])}。计算产物完整性检查通过{audit['passed_count']}项。",
        '## 1. 研究目标与当前结论',
        '本研究的核心问题是：代谢重编程与表观遗传重塑是否形成代偿关系，使双基因干预揭示单基因筛选不能发现的脆弱性，并在结直肠癌奥沙利铂耐药背景中恢复敏感性。计算侧负责提出、排序和解释这一假说；新实验与机制验证由医学院承担。',
        f"研究对象已严格限定为**代谢酶 × 直接组蛋白修饰酶**。现有{annotation['n_model_genes']}个建模基因包含38个严格代谢酶、15个直接组蛋白修饰酶，570对中19对已有HT29组合测量，剩余**{len(core)}对核心候选**。由同一GNN评分形成**Top500→Top30→Priority8→3对首轮代表候选**。BRD3等reader、MEN1等辅因子和其他表观调控因子仅在扩展表中汇报。",
        f"主模型为 **{MODEL_NAMES[selected]}**，由候选GNN架构在验证集上选择；全部方法中验证表现最好的方法为 **{MODEL_NAMES[best]}**。主GNN的未见配对测试Spearman为 **{mean_sd(selected_primary.spearman)}**，RMSE为 **{mean_sd(selected_primary.rmse)}**。该GNN未稳定优于全部基线，模型选择与泛化评价分别报告。",
        f"更严格的泛化结果需要优先注意：主GNN未见基因测试Spearman为 **{mean_sd(choose(selected,'gene_holdout').spearman)}**，跨细胞迁移为 **{mean_sd(choose(selected,'transfer_ht29').spearman)}**，接近零的结果不支持可靠外推。这限制了Top500作为预测命中名单的证据强度；当前仅作为待检验假说目录与后续采集真实反馈的起点。",
        f"严格核心的已测空间仅19对，本次主测试包含{main_core_test}对，主动学习池包含{int(curves.core_pool_n.iloc[0])}对。核心域的样本量不足限制了对新组合发现能力的评价；全量域性能不能替代该核心域的独立验证。**当前成果适合干实验阶段汇报与探索性候选交接，不证明已发现合成致死机制或逆转耐药。**",
        '**交付判断：可用于计算阶段汇报与探索性候选交接；现有计算证据尚不足以支撑高可信耐药逆转候选的最终定稿。** 这里区分计算工作是否完成与候选是否获得足够支持，不将合作方尚未开展的新实验记作计算流程故障。',
        '## 2. 数据、分析与产出概览',
        table(['研究环节', '状态', '实际内容', '结果文件', '解释范围'], requirements),
        '当前交付包含500个基因组合和1000条文库容量草案，guide序列设计仍待衔接。模型输出GI排序分数与不确定性；合成致死概率的校准需要适用背景中的真实标签，目前尚未完成。',
        '## 3. 对象定义、来源与多源表示',
        img('01_annotation_scope', '图1．严格研究对象与HT29单基因测量覆盖。类别可能重叠；广义表观调控集合不是核心组蛋白修饰酶集合。'),
        '15个直接组蛋白修饰酶：'+ '、'.join(annotation['direct_histone_genes'])+'。每个基因的催化功能、UniProt来源与多功能蛋白注意事项保存在[直接组蛋白酶审查表](../data/evidence/direct_histone_enzyme_audit.csv)。KDM8保留机制争议，不因名称直接纳入经典组蛋白去甲基化酶。未纳入严格代谢集合并不意味着不存在任何代谢作用；本项目不是全基因组催化酶目录。',
        img('14_input_coverage', '图14．实际模型输入覆盖。缺少通量映射不表示零通量；缺少表达值保留为缺失；组蛋白可映射但没有作者调用峰的情况另记为0。'),
        f"整合后每个基因有{inp['features']}维特征，包含Reactome通路成员、严格功能标记、独立证据数值及显式缺失标记。连续特征的填补和标准化统计量在每个训练基因集合内估计；缺失标准化为0，并保留缺失指示。TCGA按患者汇总；单细胞先按全部基因UMI做CP10K与log1p，再在患者内、患者间汇总，避免细胞数多的患者主导。肿瘤组织中的上皮细胞沿用作者标签，不等同本项目重新鉴定的恶性细胞。",
        '**通量输入。** [Herrmann等的原研究](https://doi.org/10.3390/cancers13164130)与[固定代码来源](https://github.com/HAHerrmann/Hct116_DrugRes)提供HCT116亲本/OxR热力学模型的上下界。对每个状态按其HCC_biomass下界归一化，取log1p最大绝对界，再按反应—基因映射等权汇总。共有4723个反应范围，192个反应映射到134个建模基因。这里重分析的是clean快照，未声称复现论文全部图；GPR成员关系未被当成敲除模拟。通量上下界也不是统计置信区间。',
        '**表观图谱。** [GSE138654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138654)含H3K4me3、H3K27ac、H3K27me3在HCT116/HCT116OxR中的六个已处理峰表。使用hg19 RefGene各转录本TSS±2kb窗口并集，提取不重复峰数与最大富集。仅保留描述性状态差，不做无重复的差异结合检验；峰数及富集差异可能受样本与峰调用影响，没有跨样本spike-in校准。图谱不能指定某候选修饰酶的直接靶基因。该研究的耐药株、通量研究的耐药株与合作方的HCT116-R分别保留身份。',
        '图网络实际覆盖：\n\n'+table(['网络','无向边数','连通基因数'],graph_summary),
        'STRING物理网络按0.7阈值二值化；Reactome为共通路关联；BioPlex为HCT116 AP-MS关联。它们不提供候选双扰动结果。全部外部注释图在训练与推理中固定；未见基因仍可拥有可公开获取的特征和网络边，因此这里是标签留出的泛化评估，不是无任何注释的新基因零样本实验。',
        '## 4. 建模与评估协议',
        f"训练来源及测量量：\n\n"+table(['来源','测量记录'], [[name,len(group)] for name,group in full_pairs.groupby('dataset')]),
        '监督目标是来源标尺上的`GI = joint − single_a − single_b`。不同来源先用训练记录的均值和标准差做标准化；未知来源仅比较排序，不强行换算原始误差。该加性残差不是作者统计GI，也不是已验证SL标签。',
        f"HT29已测GI代理量中，{100*proxy['fraction_below_zero']:.2f}%低于0，中位数为{proxy['median']:.3f}。因此负值本身几乎没有筛选特异性。候选同时报告相对于HT29来源分布的经验百分位（不是P值）、已测single、重构joint与独立证据。两个独立预测头的代数不一致RMS为{diagnostics['twohead_inconsistency_rms']:.3f}；主GNN只学习GI，joint在具备两臂single时代数重构。",
        '主GNN将独立基因特征映射到共享表示，对外部图做一跳传播并使用可学习门控；以两个基因表示的和、绝对差、逐元素积组成交换对称pair表示，再接细胞背景与GI回归头。它不消费HT29原筛选single值作为GI输入，因而能以同一预测器覆盖551对核心候选。NAIAD为适配本任务的重实现，线性、NAIAD、双头和XGBoost基线可使用已测single；因此跨类模型差异包含输入及结构差异，不作为纯图因果消融。',
        f"隐藏维度{config['hidden_dim']}，dropout {config['dropout']}，学习率{config['learning_rate']}，batch {config['batch_size']}，最多{config['max_epochs']}轮、早停耐心{config['patience']}。pooled训练平衡K562/Jurkat/HT29来源，权重为0.25/0.25/0.5。分割种子{config['partition_seed']}；主训练种子{config['primary_seeds']}，次级协议种子{config['secondary_seeds']}。随机种子SD主要反映优化随机性。",
        table(['协议','目标训练','验证','校准','测试','pooled训练','被隔离的HT29跨组记录'],partition_rows),
        '配对留出同时排除其他训练来源中相同的保留pair；基因留出要求pair两臂属于同一分组，跨组pair被隔离；锚点留出保持保留锚点不进入训练标签。跨细胞迁移允许源与目标pair身份重叠，另报告pair-unseen子集。模型选择只用主验证集均值RMSE；主校准和主测试标签不用于该选择。公开数据此前已被查看，这次重跑仍属于重复内部验证。',
        '## 5. 实际模型表现',
        img('02_benchmark_gi','图2．未见配对GI预测。小点为训练种子；菱形及横线为均值±训练种子SD。'),
        benchmark_table(metrics, protocol='pair_holdout'),
        img('03_benchmark_joint','图3．联合表型表现。主GNN的joint由实测single和预测GI重构，不能把高joint相关性全部归因于组合交互学习。'),
        benchmark_table(metrics, target='joint', protocol='pair_holdout'),
        img('04_generalization','图4．不同留出协议的排序相关性。协议测试样本不同，灰格表示该架构未注册或指标未定义，不代表零性能。'),
        '所选主GNN各协议：\n\n'+table(['协议','n','Spearman','RMSE'],performance_rows),
        '严格核心子集的实际测试量：\n\n'+table(['协议','严格核心n','是否达到5对指标计算门槛'],[[p,int(g.n_test.iloc[0]),bool(g.metrics_evaluated.iloc[0])] for p,g in subset_counts[subset_counts.subset.eq('core_metabolism_histone')].groupby('protocol')]),
        '核心子集可计算的全部结果：\n\n'+benchmark_table(metrics,subset='core_metabolism_histone'),
        '广义表观调控子集及pair-unseen结果保留于[完整指标表](benchmarks/metrics.csv)。小样本中相关系数与Top-k容易波动；没有可计算的行时不补零。图模型的收益应查看同一数据范围下STRING或多图与对应去图模型的比较，不能将pooled与target的差异只归因于图。新增多源特征在本轮统一接入，未以本轮结果单独证明每一类特征的因果贡献。',
        img('05_metric_uncertainty','图5．主GNN每个拟合的95% dyadic基因簇bootstrap指标区间，与训练种子SD分别报告。'),
        img('06_prediction_intervals','图6．名义90%预测区间的实际测试覆盖与宽度；这与相关系数的置信区间不同。'),
        '最终候选模型将初始训练与已报告测试集合并重新拟合，保留原验证和校准集；该合并集不再被称为独立测试。五个最终拟合预测均值构成候选评分，固定校准集的绝对残差给出共同区间半径。双非锚点设计、新基因组合和耐药状态改变可能破坏覆盖率。',
        '## 6. 主动学习回顾性小闭环',
        f"模型为{MODEL_NAMES[selection['selected_al_model']]}，由验证集在target STRING与target去图模型中选择。训练池初始揭示5%，每轮新增2%，共12查询轮；4策略×5种子，每条序列包含初始拟合及12次更新。不同策略共享相同种子下的初始样本、验证集、校准集与测试集。随机/预测效应/不确定性/预期改进四策略同预算比较；测试指标仅在查询后记录，不决定下一轮或模型。",
        img('07_active_learning','图7．全量已测池低GI目标的发现与固定留出表现。横轴不含固定验证集成本，表格额外列出该成本。阴影为训练/初始抽样种子SD。'),
        table(['策略','已揭示池标签','加验证标签总数','找到最低100个GI的个数','测试Spearman','严格核心召回'],al_rows),
        img('13_active_learning_subsets','图13．宽表观子集与直接组蛋白核心子集分开评价。核心池少于20对时，按预设门槛不计算命中率；空图不表示失败率为100%。'),
        '预期改进使用MC-dropout均值与SD、已观测GI的低10%分位点，计算数值稳定的Gaussian log-EI排序；24次dropout的差异是模型启发式不确定性。真实未查询标签在训练输入中屏蔽；标签揭示文件、查询分数和每轮索引均可复核。当前池是基础生长测量，回顾性全量域收益不能证明核心域或耐药背景的前瞻性收益。',
        '## 7. 独立研究外测',
        img('08_external_coverage','图8．独立研究的完整实测覆盖、可评价模型基因域与未见pair子集。未覆盖不是阴性；实测未过hit阈值也不是已证明无效。'),
        img('09_external_performance','图9．主GNN在独立研究已覆盖pair上的排序与AUPRC。灰刻度为各来源hit比例。'),
        table(['来源','子集','n','来源hit数','Spearman','AUPRC','hit比例'],ext_rows),
        '外测使用[In4mer](https://doi.org/10.1038/s41467-024-47795-3)和[Dede](https://doi.org/10.1186/s13059-020-02173-2)完整实测表，保留各自效应量hit定义。冻结预测器不以这些外测分数调参。不同来源的标准化尺度不同，因此不报告跨源原始RMSE。COLO1已进入训练，不能当作独立外测。当前外测都是基础生长条件，不提供候选在奥沙利铂耐药背景中的双扰动结果。',
        '## 8. 耐药状态与机制证据',
        img('10_resistance_expression','图10．GSE76092未处理HTOXAR3−HT29主对比的基因表达与CAMERA通路关联。通路方向是集合层面，不等于每个成员的方向。'),
        table(['对比','基因数','BH<0.05','BY<0.05','BH<0.05且|log2FC|≥1'],[[r.contrast,r.n_genes,r.n_bh_005,r.n_by_005,r.n_high_effect_bh] for r in contrasts.itertuples()]),
        'GSE76092主对比是独立3对3未处理样本，采用limma稳健经验贝叶斯，BH与BY校正均保留；CAMERA和分方向ORA以实际被检验基因为背景。急性用药对比单列，亲本与耐药株剂量不同，交互项不解释为同剂量敏化。HCT116/SW480的GSE299428用各背景3对3计数模型；GSE289766与GSE42387仅报告其可支持的效应量。GSE254136只有一个PBS对照，单基因药筛结果无推断性FDR。',
        '表达上调和下调都可提供状态关联；没有根据方向直接给出“抑制即可逆转耐药”的结论。候选证据卡同时保留DepMap亲本单敲（本快照SW480无可用值）、TCGA/单细胞表达、HCT116/OxR通量与组蛋白峰、HT29及HCT116/SW480耐药表达、单基因药筛与外部同pair证据。',
        '## 9. Top500、首轮代表组合与交接',
        f"Top500全部来自{len(core)}对严格核心候选，按预测GI从低到高排序，不用扩展类别补足500。Top500中两臂HT29 single均有测量的组合为{int(top500.observed_HT29_single_count.eq(2).sum())}对；其他组合的raw joint保持缺失。Top30进一步限定两臂single已测，并按预先写入的证据层级排序、限制每基因最多3次。Priority8沿同一规则选取，首轮3对要求基因不重复。",
        f"551对核心候选中，{int(core.anchor_gene_count.eq(0).sum())}对为双非锚点配对，超出了原HT29筛选主要配对设计。原筛选{100*diagnostics['ht29_anchor_pair_fraction']:.2f}%的已测组合包含锚点，因此这一设计外推是单独的主要限制。",
        '排序依次审查：control-scaled single是否缺失、强single数量、raw joint是否为负、DepMap亲本单敲负效应、GI区间层、种子SD层、HCT116/SW480耐药关联、HT29耐药表达关联、GI数值及pair key。强single以作者负控0/正控−1标尺的≤−0.5为默认，并报告−0.3/−0.7敏感性。该阈值是动态范围筛查，不是临床毒性或正式essentiality概率。DepMap≤−0.5为辅助分层，来源不同，不与作者scaled值混为同一测量。',
        '首轮3对计算提名：\n\n'+table(['顺序','组合','预测GI','名义90%区间','raw joint','来源GI百分位','强single数','目标耐药表达关联','证据层'],candidate_rows),
        f"首轮3对中，达到预设HCT116/SW480耐药表达关联条件的有{int(reps.target_resistance_association_count.gt(0).sum())}对。该条件为至少一臂在相应GSE299428背景达到FDR<0.05且|log2FC|≥1。未达到不等于无作用，但这些提名的依据主要是严格机制分类、single动态范围、raw joint与预测稳定性，不能称为已有耐药逆转支持的组合。",
        f"Top30中有独立组合测量覆盖的候选为{evidence_review['top30_external_measured_pairs']}/{len(top)}。" + ('由于没有实测覆盖，外部命中率当前不可评价，不能把它解读成30对全部验证失败。' if evidence_review['top30_external_measured_pairs'] == 0 else '来源hit仅按已覆盖候选评价。') + '所有提名仍为待验证组合假说。',
        '### 证据计数的含义',
        table(['项目','实际结果','如何解释'], [
            ['独立组合测量覆盖', f"{evidence_review['top30_external_measured_pairs']}/{len(top)}对", '没有实测覆盖时，外部命中率不可评价；不能计作阴性验证'],
            ['首轮候选的目标耐药表达关联', f"{evidence_review['representatives_with_target_expression_association']}/{len(reps)}对达到阈值", '当前关联证据不足；表达未达阈值也不能证明功能无效'],
            ['Top30强单基因效应筛查', f"{evidence_review['top30_without_strong_single_effect']}/{len(top)}对的两臂均未触发阈值", '这是动态范围筛查项，不是无效组合数，也不证明肿瘤选择性'],
            ['核心主动学习池', f"{evidence_review['active_learning_core_measured_pool']}对实测", '不足预设20对评价门槛，核心召回率未计算，不是0%'],
        ]),
        f"进一步复核Top500，{evidence_review['top500_with_target_expression_association']}对达到目标耐药表达阈值的组合均涉及{'、'.join(evidence_review['associated_metabolic_genes'])}；其中{evidence_review['associated_pairs_excluded_by_missing_HT29_singles']}对因缺少HT29单基因测量，被首轮入选规则排除。当前Top30中每对的两臂在HCT116/SW480的HD/MD四个对比均有表达与FDR值，未达阈值并非这些表达字段缺失所致。",
        '这体现了现有筛选规则的取舍：优先保证HT29单基因效应可审查，同时失去了现有目标耐药表达线索。相关组合已列入[单基因证据待补充的组合](handoff/resistance_expression_followup.csv)，供补充信息后复核；此次计数审查保持原阈值、预测和候选顺序。15对共享HCCS在同一研究中的基因表达线索，不构成15项独立的组合效应证据；表达方向也不能直接推出干预方向或耐药逆转。完整计数见[证据复核记录](handoff/evidence_review.json)。',
        img('11_candidate_predictions','图11．Top30的GI区间与作者控制标定single。负GI不能单独确认SL；同时审查single动态范围和raw joint。'),
        table(['Top30','组合','预测GI','raw joint','强single数','Priority8'],top_rows),
        img('12_priority_evidence','图12．Priority8的HT29耐药表达和单基因药物筛选线索，均不等于双基因相互作用。'),
        '候选文件：[Top500](candidates/top500.csv)、[全部551对核心候选及完整证据](candidates/core_candidates.csv)、[Top30](candidates/top30.csv)、[Priority8](candidates/priority8.csv)、[逐对证据卡](candidates/candidate_evidence_cards.md)、[阈值敏感性](candidates/single_threshold_sensitivity.csv)、[排序规则](candidates/ranking_policy.json)。Top500是便于浏览的紧凑表，完整原始字段保留在核心表。',
        '文库衔接交付：[500对基因靶点](handoff/mini_pool_gene_pairs.csv)、[1000条容量草案](handoff/library_capacity_plan.csv)、[首轮代表组合](handoff/first_round_representatives.csv)、[真实测量回传模板](handoff/measured_feedback_template.csv)。容量草案包含500主候选、450第二设计预留与50对照预留，**不含guide序列，也不宣称这些预留足以完成全部必要对照**。若每基因多guide、重复与单基因对照需要更多容量，应由合作方确定最终组合覆盖与文库大小；当前不是合成订单。',
        '## 10. 可复现交付与尚待完成事项',
        '完整数据、模型权重与逐次训练产物保存在计算服务器。本地保留阅读所需的代码、报告、图表、核心候选、证据和汇总表；全流程复现使用服务器完整目录。报告中的阅读材料使用相对链接。依赖见[requirements](../requirements.txt)，实际版本与GPU见[runtime](runtime.json)，输入及训练代码校验值见[input_fingerprints](input_fingerprints.json)，检查明细见[integrity_audit](integrity_audit.json)。数值分析没有调用LLM评分。',
        '复现顺序见[README](../README.md)。源文件URL、字节数和SHA256记录位于`data/evidence/requested_sources/*.provenance.json`及各来源manifest。模型推理所需数组、独立特征/图和最终权重保存在服务器；从原始文件重建还需要Python、R/limma及所列数据快照。本地精简目录不包含完整训练数据。',
        '**已完成的计算工作：** 研究对象定义、公开输入接入、模型拟合、留出与独立外测、回顾性主动学习、候选排序和证据整理。**计算证据仍需补强：** 新基因与跨背景预测尚不可靠，严格核心评价样本少，首轮筛选依赖HT29单基因测量并排除了现有目标耐药表达线索；当前也未校准真实SL概率。因此可交付阶段性成果与探索性假说，尚不足以将首轮名单作为高可信耐药逆转候选定稿。**合作方后续工作：** 明确耐药株身份、衔接guide设计并提供匹配条件的测量，以支持前瞻性更新与机制验证。当前判断不以湿实验是否已经成功替代对计算结果本身的评价。',
    ]
    destination = results / 'MetaGene实验汇报.md'
    destination.write_text('\n\n'.join(text) + '\n')
    print(destination)


if __name__ == '__main__':
    import sys
    if '--legacy-gi' in sys.argv:
        sys.argv.remove('--legacy-gi')
        main()
    else:
        from make_sl_report import main as current_main
        current_main()
