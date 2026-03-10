# Author: Yiming Shi

# ---- Load Libraries ----

if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly=TRUE)) install.packages("ggrepel")

library(ggplot2)
library(ggrepel)

# ---- Parameters ----
counts_file <- "DEG_list_Het_basal vs WT_basal.csv"

padj_cutoff <- 0.05
log2fc_cutoff <- 0
top_labels <- 30

# ---- Load counts ----
counts_res <- read.csv(counts_file)
counts_res$padj[is.na(counts_res$padj)] <- 1

# ---- Define significance ----
counts_res$significant <- ifelse(counts_res$padj < padj_cutoff, "Significant", "NotSig")
count_significant <- sum(counts_res$significant == "Significant", na.rm = TRUE)
print(paste("Significant Gene Number: ", count_upregulated))

# ---- Define Upregulated & Downregulated ----
counts_res$regulated <- ifelse(
  counts_res$significant == "Significant",
  ifelse(
    counts_res$padj < padj_cutoff & counts_res$log2FoldChange > log2fc_cutoff,
    "Upregulated",
    "Downregulated"
  ),
  "Others"
)
count_upregulated <- sum(counts_res$regulated == "Upregulated", na.rm = TRUE)
print(paste("Upregulated Gene Number: ", count_upregulated))

# ---- Calculate -log10padj ----
counts_res$neglog10padj <- -log10(counts_res$padj)

# Pick top N genes for labeling (most significant by padj)
top_symbols <- counts_res$symbol[order(counts_res$padj)][1:min(top_labels, nrow(counts_res))]

# ---- Volcano plot ----
volcano_plot <- ggplot(counts_res, aes(x = log2FoldChange, y = neglog10padj)) +
  geom_point(aes(color = regulated), alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("Upregulated" = "blue", "Downregulated" = "pink", "Others" = "grey")) +
  geom_vline(xintercept = log2fc_cutoff,
             linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "gray40") +
  geom_text_repel(data = subset(counts_res, symbol %in% top_symbols),
                  aes(label = symbol),
                  size = 3, max.overlaps = Inf) +
  labs(
    title = "Het vs WT (NC only) - Upregulated Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Upregulated"
  ) +
  xlim(-6, 6) +        # FIXED X-AXIS RANGE
  ylim(0, 10)           # FIXED Y-AXIS RANGE

# ---- Save File ----

png(file.path(getwd(), "volcano_Het_vs_WT_NC_Upregulated.png"),
    width = 900, height = 1600, res = 150)
print(volcano_plot)
dev.off()