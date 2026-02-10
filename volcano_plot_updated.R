#############################################
# Volcano Plot: Het-NC vs WT-NC (NC only)
# Matches collaborator: raw log2FC, fixed axis ranges
#############################################

if (!requireNamespace("DESeq2", quietly=TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly=TRUE)) install.packages("ggrepel")

library(DESeq2)
library(ggplot2)
library(ggrepel)

# ---- Parameters ----
counts_file <- "gene_count_WT-NC_Het-NC_GeneID.csv"
outdir <- "RNAseq_results_volcano_raw_fixedaxis"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

padj_cutoff <- 0.05
log2fc_cutoff <- 1
top_labels <- 30

# ---- Load counts & set gene_name as rownames ----
counts_raw <- read.csv(counts_file)
rownames(counts_raw) <- make.unique(counts_raw$gene_name)

# Keep only needed columns
col_order <- c("WT21","WT31","Het11","Het21","Het31")
counts_df <- counts_raw[, col_order]

# ---- Sample metadata ----
sample_info <- data.frame(
  row.names = col_order,
  condition = c("WT_NC","WT_NC","Het_NC","Het_NC","Het_NC")
)
sample_info$condition <- factor(sample_info$condition, levels=c("WT_NC","Het_NC"))

# ---- Run DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = sample_info,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# ---- Get raw DESeq2 results (no shrinkage) ----
res_raw <- results(dds, contrast=c("condition","Het_NC","WT_NC"))

# ---- Prepare dataframe ----
res_df <- as.data.frame(res_raw)
res_df$gene <- rownames(res_df)
res_df$padj[is.na(res_df$padj)] <- 1
res_df$neglog10padj <- -log10(res_df$padj)

# Define significance
res_df$significant <- ifelse(res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) > log2fc_cutoff,
                             "Significant", "NotSig")

# Pick top N genes for labeling (most significant by padj)
top_genes <- res_df$gene[order(res_df$padj)][1:min(top_labels, nrow(res_df))]

# ---- Volcano plot ----
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = significant), alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("NotSig" = "black", "Significant" = "red")) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff),
             linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "gray40") +
  geom_text_repel(data = subset(res_df, gene %in% top_genes),
                  aes(label = gene),
                  size = 3, max.overlaps = Inf) +
  labs(
    title = "Het vs WT (NC only) - Raw DESeq2 log2FC",
    x = "Log2 Fold Change (raw)",
    y = "-Log10 Adjusted P-Value",
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face="bold")
  ) +
  xlim(-6, 6) +        # FIXED X-AXIS RANGE
  ylim(0, 9)           # FIXED Y-AXIS RANGE

# ---- Save plot ----
png(file.path(outdir, "volcano_Het_vs_WT_NC_raw_fixedaxis.png"),
    width = 800, height = 1200, res = 150)
print(volcano_plot)
dev.off()

message("✅ Volcano plot with RAW log2FC + fixed axis range saved in ", outdir)