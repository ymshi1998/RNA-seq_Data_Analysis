### ReadMe
#### Overview
This is an RNA-seq data analysis R script, which takes the DEG (Differentially Expressed Genes) List csv as input to draw a volcano plot that displays significant & upregulated or downregulated genes.

The DEG list csv file named "DEG_list_Het_basal vs WT_basal.csv" compares Het (Heterozygous) and WT (Wild Type) mouse samples, whereas Het samples carry a knockout of genes involved in pulmonary functions. The data will not be provided, but the code can be used as sharing and learning.

#### Functional Details
The code uses padj < 0.05 as cutoff to define significant genes. Based on that, the code uses log2FoldChange > 0 as cutoff to define upregulated genes, and log2FoldChange < 0 as cutoff to define downregulated genes. The top 30 high padj valued genes' names will be displayed.

#### Running Result
This is a sample output we got from running the R script.
<img width="800" height="1200" alt="volcano_Het_vs_WT_NC_Upregulated" src="https://github.com/user-attachments/assets/2dda4431-4cb3-4f47-a1b1-b6ba666180ed" />
