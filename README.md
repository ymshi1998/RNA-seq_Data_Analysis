### ReadMe
#### Overview
This is an RNA-seq data analysis R script, which takes the DEG (Differentially Expressed Genes) List csv as input to draw a volcano plot that displays significant & upregulated or downregulated genes.

The DEG list csv file named "DEG_list_Het_basal vs WT_basal.csv" compares Het (Heterozygous) and WT (Wild Type) mouse samples, whereas Het samples carry a knockout of genes involved in pulmonary functions. The data will not be provided, but the code can be used as sharing and learning.

#### Functional Details
The code uses padj < 0.05 as cutoff to define significant genes. Based on that, the code uses log2FoldChange > 0 as cutoff to define upregulated genes, and log2FoldChange < 0 as cutoff to define downregulated genes. The top 30 high padj valued genes' names will be displayed.

#### Running Result
This is a sample output we got from running the R script.
<img width="900" height="1600" alt="volcano_Het_vs_WT_NC_Upregulated 2" src="https://github.com/user-attachments/assets/ad80283d-2852-4bc6-a9f7-f4ea460e8729" />
