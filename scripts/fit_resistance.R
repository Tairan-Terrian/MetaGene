#!/usr/bin/env Rscript
# Primary inference uses only the six untreated samples.
args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
if (length(args) >= 2 && nzchar(args[2])) .libPaths(c(args[2], .libPaths()))
suppressPackageStartupMessages(library(limma))
expr <- as.matrix(read.csv(file.path(outdir, "GSE76092_log2_intensity.csv"), row.names=1, check.names=FALSE))
meta <- read.csv(file.path(outdir, "GSE76092_sample_mapping.csv"), row.names=1, check.names=FALSE)
intensity <- read.csv(file.path(outdir, "GSE76092_gene_intensity.csv"), row.names=1, check.names=FALSE)
pathways <- read.csv(file.path(outdir, "reactome_tested_membership.csv"), stringsAsFactors=FALSE)
path_index <- lapply(split(pathways$gene, pathways$pathway_id), function(g) match(g, rownames(expr)))
path_name <- setNames(pathways$pathway_name, pathways$pathway_id)
all_results <- list()
all_camera <- list()
all_ora <- list()
summaries <- list()

analyse <- function(samples, contrasts, analysis_role) {
    y <- expr[, samples, drop=FALSE]
    groups <- factor(meta[samples, "group"])
    design <- model.matrix(~0+groups)
    colnames(design) <- levels(groups)
    ct <- makeContrasts(contrasts=unname(contrasts), levels=design)
    colnames(ct) <- names(contrasts)
    fit <- eBayes(contrasts.fit(lmFit(y, design), ct), trend=TRUE, robust=TRUE)
    for (j in seq_len(ncol(ct))) {
        key <- colnames(ct)[j]
        tab <- topTable(fit, coef=j, number=Inf, sort.by="none", confint=TRUE, adjust.method="BH")
        result <- data.frame(
            gene=rownames(tab), study="GSE76092", analysis_role=analysis_role,
            contrast=key, log2fc=tab$logFC, ci95_low=tab$CI.L, ci95_high=tab$CI.R,
            average_log2_intensity=tab$AveExpr, moderated_t=tab$t,
            pvalue=tab$P.Value, padj_bh=tab$adj.P.Val, padj_by=p.adjust(tab$P.Value, "BY"),
            log_odds_de=tab$B, n_samples=ncol(y), residual_df=fit$df.residual,
            stringsAsFactors=FALSE)
        result$significant_bh_005 <- result$padj_bh < 0.05
        result$significant_by_005 <- result$padj_by < 0.05
        result$high_effect_bh <- result$significant_bh_005 & abs(result$log2fc) >= 1
        result <- cbind(result, intensity[result$gene, , drop=FALSE])
        write.csv(result, file.path(outdir, paste0(key, "_gene_effects.csv")), row.names=FALSE)
        all_results[[key]] <<- result
        summaries[[key]] <<- data.frame(
            contrast=key, role=analysis_role, n_samples=ncol(y), n_genes=nrow(result),
            n_bh_005=sum(result$significant_bh_005), n_by_005=sum(result$significant_by_005),
            n_high_effect_bh=sum(result$high_effect_bh),
            n_up_high_effect_bh=sum(result$high_effect_bh & result$log2fc>0),
            n_down_high_effect_bh=sum(result$high_effect_bh & result$log2fc<0))

        camera_result <- camera(y, path_index, design=design, contrast=ct[,j],
                                inter.gene.cor=0.01, sort=FALSE)
        camera_result <- data.frame(pathway_id=rownames(camera_result), camera_result,
                                    row.names=NULL, check.names=FALSE)
        camera_result$pathway_name <- unname(path_name[camera_result$pathway_id])
        camera_result$contrast <- key
        camera_result$method <- "CAMERA competitive, fixed inter-gene correlation 0.01"
        all_camera[[key]] <<- camera_result

        for (direction in c("up", "down")) {
            selected <- result$gene[result$high_effect_bh &
                                   if (direction=="up") result$log2fc>0 else result$log2fc<0]
            ora <- do.call(rbind, lapply(names(path_index), function(pid) {
                genes <- rownames(expr)[path_index[[pid]]]
                hit <- intersect(genes, selected)
                data.frame(pathway_id=pid, pathway_name=unname(path_name[pid]),
                           contrast=key, direction=direction, n_pathway=length(genes),
                           n_selected=length(selected), n_universe=nrow(expr), n_overlap=length(hit),
                           pvalue=phyper(length(hit)-1, length(genes), nrow(expr)-length(genes),
                                        length(selected), lower.tail=FALSE),
                           overlap_genes=paste(hit, collapse=";"), stringsAsFactors=FALSE)
            }))
            ora$padj_bh <- p.adjust(ora$pvalue, "BH")
            ora$selection_rule <- "gene BH<0.05 and absolute log2FC>=1"
            all_ora[[paste(key, direction)]] <<- ora
        }
    }
}

analyse(rownames(meta)[meta$treatment=="untreated"],
        c(resistant_vs_parental_untreated="HTOXAR3_untreated-HT29_untreated"),
        "primary_acquired_resistance_association")
analyse(rownames(meta)[meta$treatment %in% c("untreated", "OXA")],
        c(parental_oxa_response="HT29_OXA-HT29_untreated",
          resistant_oxa_response="HTOXAR3_OXA-HTOXAR3_untreated",
          state_by_oxa_interaction="(HTOXAR3_OXA-HTOXAR3_untreated)-(HT29_OXA-HT29_untreated)"),
        "secondary_unequal_oxa_doses")
write.csv(do.call(rbind, all_results), file.path(outdir, "GSE76092_gene_effects.csv"), row.names=FALSE)
write.csv(do.call(rbind, all_camera), file.path(outdir, "reactome_camera.csv"), row.names=FALSE)
write.csv(do.call(rbind, all_ora), file.path(outdir, "reactome_ora.csv"), row.names=FALSE)
write.csv(do.call(rbind, summaries), file.path(outdir, "contrast_summary.csv"), row.names=FALSE)
capture.output(sessionInfo(), file=file.path(outdir, "R_session_info.txt"))
print(do.call(rbind, summaries))
