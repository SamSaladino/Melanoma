## 🔬 Differential Expression Analysis Summary

This analysis identifies **Differentially Expressed Genes (DEGs)** using RNA-seq count data and the **DESeq2** R package. The goal is to compare gene expression between conditions (e.g., *Primary Tumor* vs. *Metastatic Tumor*) to detect genes with statistically significant expression changes.

### 📦 Tools Used
- R
- DESeq2
- ggplot2

### 🧪 Input Data
- **Gene list** (`FPKM.txt`): Extracted and cleaned to get gene names.
- **Count matrix**: Raw expression counts across samples.
- **Sample metadata**: Annotated with sample types (e.g., `Primary Tumor`, `Metastatic Tumor`).

### 🔧 Analysis Steps
1. Prepare count matrix with gene names as row names.
2. Create `DESeqDataSet` object for modeling.
3. Filter genes with zero counts.
4. Set reference condition (`Primary Tumor`) for comparisons.
5. Run DESeq2 to normalize data and test for differential expression.
6. Summarize results, including significant genes (`alpha = 0.01`).
7. Visualize results using a volcano plot.
8. Export DEGs to `DEGs.csv`.

### 📊 Output
- CSV file with DEGs and statistics (log2 fold change, adjusted p-values).
- Volcano plot highlighting significant genes.
