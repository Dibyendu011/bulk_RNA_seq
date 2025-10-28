setwd("/path/to/gene_counts.txt")
getwd()

# calculating RPKM
library(tidyverse)
out_FeatCount_mut_rep_count = read.delim("gene_counts.txt", skip=1, check.names = FALSE)

# edgeR
library(edgeR)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# Prepare count matrix
counts_data = out_FeatCount_mut_rep_count[, 7:10]
colnames(counts_data) = c("conA_rep1", "conA_rep2", "conB_rep1", "conB_rep2")
rownames(counts_data) = out_FeatCount_mut_rep_count$Geneid
counts_data = counts_data[!is.na(rownames(counts_data)), ]
counts_data = counts_data[!duplicated(rownames(counts_data)), ]

# Create group information
group = factor(c("conA", "conA", "conB", "conB"))
# Build DGEList
y = DGEList(counts = counts_data, group = group)
# Filter low-expression genes
keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes = FALSE]
# Normalize
y = calcNormFactors(y)
# Estimate dispersions
y = estimateDisp(y)
# Exact test
et = exactTest(y)

# Results
results_edgeR = topTags(et, n = Inf)$table
write.csv(results_edgeR, file = "edger_results.csv")

# Variance-stabilized logCPM for visualization
logCPM = cpm(y, log = TRUE, prior.count = 2)

# Heatmap of significant genes (FDR < 0.05)
sig_genes = rownames(results_edgeR[results_edgeR$FDR < 0.05, ])
pheatmap(logCPM[sig_genes, ],
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = data.frame(condition = group, row.names = colnames(logCPM)),
           show_rownames = FALSE,
           fontsize_col = 12,
           main = "edgeR - Heatmap of Significant Genes")

# Heatmap of top 20 genes by FDR
topN = 20
top_genes = rownames(results_edgeR)[1:min(topN, nrow(results_edgeR))]
length(top_genes) >= 2
pheatmap(logCPM[top_genes, ],
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = data.frame(condition = group, row.names = colnames(logCPM)),
           fontsize_row = 8,
           fontsize_col = 12,
           main = paste("edgeR - Top", topN, "Genes"))

# volcano plot
res_df = results_edgeR
p_threshold = 0.05
fc_threshold = 1
res_df$neg_log10_p = -log10(res_df$FDR)
res_df$significant = ifelse(res_df$FDR < p_threshold & abs(res_df$logFC) >= fc_threshold,
                             "Significant", "Not significant")
xlim_range = c(-max(abs(res_df$logFC), na.rm = TRUE),
                 max(abs(res_df$logFC), na.rm = TRUE))
plot(res_df$logFC, res_df$neg_log10_p,
     type = "n",
     xlab = "log2 Fold Change",
     ylab = "-log10 FDR",
     main = "edgeR Volcano Plot",
     cex.lab = 1.2,
     cex.main = 1.5,
     xlim = xlim_range)
points(res_df$logFC[res_df$significant == "Not significant"],
       res_df$neg_log10_p[res_df$significant == "Not significant"],
       col = "gray60", pch = 16, cex = 0.6)
points(res_df$logFC[res_df$significant == "Significant"],
       res_df$neg_log10_p[res_df$significant == "Significant"],
       col = "red", pch = 16, cex = 0.8)
abline(h = -log10(p_threshold), col = "blue", lty = 2, lwd = 2)
abline(v = c(-fc_threshold, fc_threshold), col = "blue", lty = 2, lwd = 2)
legend("topright",
       legend = c("Significant", "Not significant"),
       col = c("red", "gray60"),
       pch = 16,
       bty = "n")
sig_count = sum(res_df$significant == "Significant")
total_count = nrow(res_df)
text(par("usr")[1], par("usr")[4],
     paste("Significant:", sig_count, "/", total_count),
     pos = 4, col = "darkred", cex = 0.9)
