setwd("/path/to/bam-files")
getwd()

# Load packages
library(SGSeq)
library(pheatmap)

# 1. Define sample info for BAM files
file_bam = c("Con_A1_sorted.bam","Con_A2_sorted.bam",
              "Con_B1_sorted.bam","Con_B2_sorted.bam")
sample_name = c("conA_rep1","conA_rep2","conB_rep1","conB_rep2")
v = data.frame(sample_name, file_bam)

# 2. Extract BAM metadata
Baminfo = getBamInfo(v)

# 3. Import annotation & convert to transcript features
x = importTranscripts(file = "/path/to/GCF_000146045.2_R64_genomic_copy.gtf")
TxFeat = convertToTxFeatures(x)

# 4. Generate splice graph features (reference-based model)
sgf_ucsc = convertToSGFeatures(TxFeat)
head(sgf_ucsc)

# 5. Quantify splicing features using BAM files
sgfc_ucsc = analyzeFeatures(Baminfo, features = TxFeat)
sgfc_ucsc

# 6. Alternative: discover novel features (prediction mode)
sgfc_pred = analyzeFeatures(Baminfo)

# 7. Annotate predicted features against reference
sgfc_annot = annotate(sgfc_pred, TxFeat)

# 8. Inspect results
colData(sgfc_annot)        # sample metadata
rowRanges(sgfc_annot)      # feature genomic coordinates
counts(sgfc_annot)[1:5, ]  # raw counts (exons/junctions)
head(FPKM(sgfc_annot))     # normalized values

# 9. Visualization: splice graph for a gene
df = plotFeatures(sgfc_annot, geneID = 10, color_novel = "red")

# check how many features from each chromosome
table(seqnames(sgfc_annot))

# try
View(sgfc_annot)
sgfc_annot@rowRanges@seqinfo
sgfc_annot@rowRanges@ranges

# 10. Optional heatmap of splicing feature expression
expr_mat = assay(sgfc_annot)   # feature-level counts
top_feats = head(order(rowMeans(expr_mat), decreasing = TRUE), 50)
pheatmap(expr_mat[top_feats, ],
         cluster_rows = TRUE, cluster_cols = TRUE,
         fontsize_row = 6, fontsize_col = 10,
         main = "Top Splicing Features")
