# create a new folder "R" in /home/ibab/Desktop/DIBYENDU/s3/RA/rnaseq/trimmed_reads/bam_files/ & copy gene_counts.txt, gene_counts.txt.summary and gene_length.txt there
setwd("path/to/gene_counts.txt")
getwd()

# calculating RPKM just for practice
library(tidyverse)
out_FeatCount_mut_rep_count = read.delim("gene_counts.txt", skip=1, check.names = FALSE)

# Calculate library sizes and scaling factors
LibSize1 = sum(out_FeatCount_mut_rep_count$conA_rep1_sorted)
LibSize2 = sum(out_FeatCount_mut_rep_count$conA_rep2_sorted)
LibSize3 = sum(out_FeatCount_mut_rep_count$conB_rep1_sorted)
LibSize4 = sum(out_FeatCount_mut_rep_count$conB_rep2_sorted)

sf1 = LibSize1 / 1000000
sf2 = LibSize2 / 1000000
sf3 = LibSize3 / 1000000
sf4 = LibSize4 / 1000000

# RPM just for practice
out_FeatCount_mut_rep_count$WT1_rpm = out_FeatCount_mut_rep_count$conA_rep1_sorted / sf1
out_FeatCount_mut_rep_count$WT2_rpm = out_FeatCount_mut_rep_count$conA_rep2_sorted / sf2
out_FeatCount_mut_rep_count$WT3_rpm = out_FeatCount_mut_rep_count$conB_rep1_sorted / sf3
out_FeatCount_mut_rep_count$WT4_rpm = out_FeatCount_mut_rep_count$conB_rep2_sorted / sf4

# Gene length calculation
gene_lengths = read.delim("/path/to/gene_length.txt", header = TRUE, sep = "\t")

# in case of data point mismatch, consider gene lengths with ids matching with count matrix
out_FeatCount_mut_rep_count$Matched_Gene = NA  # Initialize with NA
for (i in 1:nrow(out_FeatCount_mut_rep_count)) {
  current_gene = out_FeatCount_mut_rep_count$Geneid[i]  
  if (current_gene %in% gene_lengths$Gene) {
    out_FeatCount_mut_rep_count$Matched_Gene[i] = gene_lengths[,2][i]/1000
  }  
}

out_FeatCount_mut_rep_count$RPKM_WT1 = out_FeatCount_mut_rep_count$WT1_rpm / out_FeatCount_mut_rep_count$Matched_Gene
out_FeatCount_mut_rep_count$RPKM_WT2 = out_FeatCount_mut_rep_count$WT2_rpm / out_FeatCount_mut_rep_count$Matched_Gene
out_FeatCount_mut_rep_count$RPKM_WT3 = out_FeatCount_mut_rep_count$WT3_rpm / out_FeatCount_mut_rep_count$Matched_Gene
out_FeatCount_mut_rep_count$RPKM_WT4 = out_FeatCount_mut_rep_count$WT4_rpm / out_FeatCount_mut_rep_count$Matched_Gene

# Save table with RPKM values
write.csv(out_FeatCount_mut_rep_count, "RPKM_table.csv", row.names = FALSE)


# DESeq2
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# loading only raw count matrix & renaming columns
counts_data = out_FeatCount_mut_rep_count[, 7:10]
colnames(counts_data) = c("conA_rep1", "conA_rep2", "conB_rep1", "conB_rep2")

# Use Geneid as rownames
rownames(counts_data) = out_FeatCount_mut_rep_count$Geneid

# Remove rows with missing or duplicated IDs
counts_data = counts_data[!is.na(rownames(counts_data)), ]
counts_data = counts_data[!duplicated(rownames(counts_data)), ]

# creating a dataset object
col_data = data.frame(
  sample = colnames(counts_data),
  condition = factor(c("conA", "conA", "conB", "conB")),
  row.names = colnames(counts_data)
)

dds = DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = ~ condition
)

# Prefilter low-count genes
keep = rowSums(counts(dds)) >= 10
dds = dds[keep, ]

# Run DESeq2
dds = DESeq(dds)
results_table = results(dds)
results_ordered = results_table[order(results_table$padj), ]
write.csv(as.data.frame(results_ordered), file = "deseq2_results.csv")

# variance stabilization
vsd = vst(dds, blind = FALSE)

# Heatmap of significant genes
# Colours in the heatmap represent log₂(TMM-normalized counts per million)
sig_genes = rownames(subset(results_ordered, padj < 0.05)) 
pheatmap(assay(vsd)[sig_genes, ],
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = col_data,
           show_rownames = FALSE,
           fontsize_col = 12,
           main = "Heatmap of Significant Genes")

# Heatmap of top 20 genes
topN = 20
top_genes = rownames(results_ordered)[1:min(topN, nrow(results_ordered))]

length(top_genes) >= 2
pheatmap(assay(vsd)[top_genes, ],
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = col_data,
           fontsize_row = 8,
           fontsize_col = 12,
           main = paste("Top", topN, "Genes"))

# volcano plot
res_df = as.data.frame(results_table)
# Replace NA padj with 1 (so they won't appear as sig)
res_df$padj[is.na(res_df$padj)] = 1
p_threshold = 0.05
fc_threshold = 1
volcano_data = data.frame(
  logFC = res_df$log2FoldChange,
  neg_log10_p = -log10(res_df$padj),
  significant = ifelse(res_df$padj < p_threshold & abs(res_df$log2FoldChange) >= fc_threshold,
                       "Significant", "Not significant")
)

xlim_range = c(-max(abs(volcano_data$logFC), na.rm = TRUE),
                 max(abs(volcano_data$logFC), na.rm = TRUE))

plot(volcano_data$logFC, volcano_data$neg_log10_p,
     type = "n",
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted p-value",
     main = "Volcano Plot - DESeq2 Results",
     cex.lab = 1.2,
     cex.main = 1.5,
     xlim = xlim_range)
     
# Non-significant points
points(volcano_data$logFC[volcano_data$significant == "Not significant"],
       volcano_data$neg_log10_p[volcano_data$significant == "Not significant"],
       col = "gray60", pch = 16, cex = 0.6)

# Significant points
points(volcano_data$logFC[volcano_data$significant == "Significant"],
       volcano_data$neg_log10_p[volcano_data$significant == "Significant"],
       col = "red", pch = 16, cex = 0.8)

# Threshold lines
abline(h = -log10(p_threshold), col = "blue", lty = 2, lwd = 2)
abline(v = c(-fc_threshold, fc_threshold), col = "blue", lty = 2, lwd = 2)

# Legend
legend("topright",
       legend = c("Significant", "Not significant"),
       col = c("red", "gray60"),
       pch = 16,
       bty = "n")

# Add counts
sig_count = sum(volcano_data$significant == "Significant")
total_count = nrow(volcano_data)
text(par("usr")[1], par("usr")[4],
     paste("Significant:", sig_count, "/", total_count),
     pos = 4, col = "darkred", cex = 0.9)
     
