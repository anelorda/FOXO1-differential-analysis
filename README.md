---
title: "RNA-seq Analysis of FOXO1 Stimulation Dataset"
author: "Anel Ordabayeva"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

# Introduction

This study investigates the impact of FOXO1 stimulation on gene expression in human cells. FOXO1 is a transcription factor known to play crucial roles in cellular processes such as apoptosis, stress resistance, and metabolism. Analysis aims to uncover the differential gene expression patterns between FOXO1-stimulated and control conditions in various states: non-stimulated, stimulated, and stimulated rest.

# Data Description

The dataset consists of RNA-seq count data from human anti-Lewis Y CAR T cells overexpressing FOXO1 or mCherry (control), generated from healthy donor peripheral blood mononuclear cells (PBMC). Control or FOXO1-expressing CD8+ CAR T cells were assessed by RNA-Sequencing under the following conditions:

-   **Control (non-stimulated, stimulated, and stimulated rest)**

    -   **Non-stimulated:** Control CD8+ CAR T cells prior to stimulation.

    -   **Stimulated:** Control CD8+ CAR T cells after 24 hours activation with 0.1 μg/ml plate-bound anti-Lewis Y.

    -   **Stimulated rest:** Control CD8+ CAR T cells after 7 days recovery post-stimulation with anti-Lewis Y.

-   **FOXO1-treated (non-stimulated, stimulated, and stimulated rest)**

    -   **Non-stimulated:** FOXO1-expressing CD8+ CAR T cells prior to stimulation.

    -   **Stimulated:** FOXO1-expressing CD8+ CAR T cells after 24 hours activation with 0.1 μg/ml plate-bound anti-Lewis Y.

    -   **Stimulated rest:** FOXO1-expressing CD8+ CAR T cells after 7 days recovery post-stimulation with anti-Lewis Y.

Each condition has three replicates, resulting in a total of 18 samples. The count data represents the expression levels of genes across these conditions.

# Results

## Data Preparation and Quality Control

We begin by loading the necessary libraries and preparing our data for analysis.

```{r}
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)  # Human database
library(enrichplot)
library(dplyr)
```

## Prepare Data

``` bash
# Download the data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263257/suppl/GSE263257%5FProceed%5Fdata%5Fhuman%5FFOXO1OE%5Fstimulation%5Frest%5FRNAseq.csv.gz

# Unzip the downloaded file
gunzip GSE263257_Proceed_data_human_FOXO1OE_stimulation_rest_RNAseq.csv.gz
```

```{r prepare_data}
#load data
data <- read.table("data_analysis/GSE263257_Proceed_data_human_FOXO1OE_stimulation_rest_RNAseq.csv", header = TRUE, sep = ",")
# Set column names using the first row for columns 1-6
colnames(data)[1:6] <- as.character(data[1, 1:6])

# Keep original column names for the rest
# Remove the first row
data <- data[-1, ]

# Reset row names
rownames(data) <- NULL

# Extract count data
countData <- data[, 7:24]  # Assuming counts start from column 7 to 24

# Use mgi_symbol as row names
rownames(countData) <- make.unique(as.character(data$mgi_symbol))

# Create sample information
sampleInfo <- data.frame(
  row.names = colnames(countData),
  condition = rep(c("control_non_stimulated", "FOXO1_non_stimulated",
                    "control_stimulated", "FOXO1_stimulated",
                    "control_stimulated_rest", "FOXO1_stimulated_rest"), each = 3),
  batch = rep(c(1, 2), each = 9)
)

sampleInfo$condition <- factor(sampleInfo$condition, levels = c("control_non_stimulated", "FOXO1_non_stimulated",
                                                                "control_stimulated", "FOXO1_stimulated",
                                                                "control_stimulated_rest", "FOXO1_stimulated_rest"))
sampleInfo$batch <- factor(sampleInfo$batch)

# Convert count data to numeric without changing row names
countData <- as.data.frame(lapply(countData, function(x) as.numeric(gsub(",", "", x))))
rownames(countData) <- make.unique(as.character(data$mgi_symbol))

# Check the structure of the data
print(str(countData))
print(head(rownames(countData)))
```

In this dataset, gene IDs appear to be already as gene symbols.

## Create DESeq DataSet and Run DESeq2

We use DESeq2 to perform differential expression analysis.

```{r deseq2}
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(countData),
  colData = sampleInfo,
  design = ~ condition
)

dds <- DESeq(dds)
```

## Differential Expression Analysis with Preserved Gene IDs

```{r diff_expression_preserved_ids}
# Function to get results for a specific comparison
get_condition_results <- function(condition1, condition2) {
  results_df <- results(dds, contrast = c("condition", condition1, condition2))
  results_df <- as.data.frame(results_df)
  results_df$gene <- rownames(results_df)  # Ensure gene IDs are preserved
  results_df <- results_df[!is.na(results_df$padj), ]
  return(results_df)
}

# Get results for each comparison
results_non_stimulated <- get_condition_results("FOXO1_non_stimulated", "control_non_stimulated")
results_stimulated <- get_condition_results("FOXO1_stimulated", "control_stimulated")
results_rest <- get_condition_results("FOXO1_stimulated_rest", "control_stimulated_rest")

# Check gene IDs
print(head(results_non_stimulated$gene))
print(head(results_stimulated$gene))
print(head(results_rest$gene))
```

## Volcano Plots

Volcano plots help visualize the distribution of differentially expressed genes.

```{r volcano_plots, fig.width=10, fig.height=8}
create_volcano_plot <- function(results_df, title) {
  ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), size = 1) +
    scale_color_manual(values = c("black", "red")) +
    labs(x = "Log2 Fold Change", y = "-log10(Adjusted p-value)", title = title) +
    theme_bw() +
    theme(legend.position = "none")
}

# Create volcano plots
plot_non_stimulated <- create_volcano_plot(results_non_stimulated, "Volcano Plot: FOXO1 vs Control (Non-stimulated)")
plot_stimulated <- create_volcano_plot(results_stimulated, "Volcano Plot: FOXO1 vs Control (Stimulated)")
plot_rest <- create_volcano_plot(results_rest, "Volcano Plot: FOXO1 vs Control (Stimulated Rest)")

# Display plots
print(plot_non_stimulated)
print(plot_stimulated)
print(plot_rest)
```
![volcano_plot_non_stimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/de93a100-12f4-4a07-a014-b173bccbe9a4)
![volcano_plot_stimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/9615e21a-4c7a-4238-8201-e8ba6ef49622)
![volcano_plot_rest](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/29290cdf-7377-4a6d-b88e-9c1cffc2503e)


According to the volcano plot and the performed DESeq differential analysis, the most differential genes are found in the non-stimulated and stimulated-rest conditions, while in the stimulated condition, the genes appear to be very close to the cutoff values.

## Principal Component Analysis (PCA)

PCA helps us visualize the overall variance in the dataset and identify potential batch effects or outliers.

```{r pca_plot, fig.width=10, fig.height=8}
# Perform variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

# Save PCA plot
ggsave("PCA_plot.png")
```
![PCA_plot](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/5f6a56c8-b33a-47df-bdb1-a663aba340c0)

According to the PCA plot, the first principal component (PC1), which explains 71% of the variance, distinctly separates the stimulated samples (both FOXO1 and control) from the other conditions. This indicates that stimulation has a major effect on gene expression profiles. The second principal component (PC2), explaining 16% of the variance, separates the non-stimulated and stimulated-rest conditions, suggesting that these states have additional, but lesser, effects on gene expression. In addition to this, they are separated by samples in each condition, meaning that there can be other factors such as experimental conditions. This interpretation highlights the dominant role of stimulation in driving gene expression changes while also acknowledging the contributions of the other conditions.

## Heatmap of Top Differentially Expressed Genes

Heatmaps provide a visual representation of expression patterns across samples for the top differentially expressed genes.

```{r heatmap, fig.width=10, fig.height=8}
create_heatmap <- function(results_df, title, dds) {
  # Select top genes
  top_genes <- head(order(results_df$padj), 50)
  top_gene_names <- results_df$gene[top_genes]
  
  # Get normalized counts
  counts <- counts(dds, normalized = TRUE)
  
  # Subset counts for top genes
  counts_subset <- counts[rownames(counts) %in% top_gene_names, ]
  
  # Log2 transform counts, adding a small value to avoid log(0)
  log_counts <- log2(counts_subset + 1)
  
  # Remove any rows with infinite values
  log_counts <- log_counts[is.finite(rowSums(log_counts)), ]
  
  # Check if we have any genes left
  if (nrow(log_counts) == 0) {
    print(paste("No valid genes for heatmap in", title))
    return(NULL)
  }
  
  # Scale rows
  scaled_counts <- t(scale(t(log_counts)))
  
  # Replace any remaining infinite values with NA
  scaled_counts[is.infinite(scaled_counts)] <- NA
  
  # Define the desired order of conditions
  desired_order <- c("control_non.stimulated_1", "control_non.stimulated_2", "control_non.stimulated_3",
                     "FOXO1_non.stimulated_1", "FOXO1_non.stimulated_2", "FOXO1_non.stimulated_3",
                     "control_stimulated_1", "control_stimulated_2", "control_stimulated_3",
                     "FOXO1_stimulated_1", "FOXO1_stimulated_2", "FOXO1_stimulated_3",
                     "control_stimulated.rest_1", "control_stimulated.rest_2", "control_stimulated.rest_3",
                     "FOXO1_stimulated.rest_1", "FOXO1_stimulated.rest_2", "FOXO1_stimulated.rest_3")
    
  # Reorder the columns of the scaled_counts matrix
  col_order <- match(desired_order, colnames(scaled_counts))
  scaled_counts <- scaled_counts[, col_order, drop = FALSE]
  
  # Create annotation dataframe
  annotation_col <- data.frame(Condition = sampleInfo$condition)
  rownames(annotation_col) <- colnames(scaled_counts)
  
  # Ensure the annotation dataframe is in the same order
  annotation_col <- annotation_col[desired_order, , drop = FALSE]
  
  # Create the heatmap
  pheatmap(scaled_counts, 
           scale = "none",  # We've already scaled the data
           show_rownames = FALSE,
           annotation_col = annotation_col,
           main = paste("Top 50 Differentially Expressed Genes -", title),
           na_col = "grey")  # Color NA values as grey
}

# Create heatmaps
create_heatmap(results_non_stimulated, "Non-stimulated", dds)
create_heatmap(results_stimulated, "Stimulated", dds)
create_heatmap(results_rest, "Stimulated Rest", dds)

```
![heatmap_non-stimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/fdc4e986-dc31-4731-8f1f-50009acb8c1d)
![heatmap_stimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/db352cac-3c3f-454a-9066-d557b1c5cbeb)
![heatma_rest](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/329ae4aa-2ecb-48aa-9ad2-7ddb1e320ae3)


According to the heatmap of the deferentially expressed genes in non-stimulated condition are mostly under-expressed in control, while in FOXO1 are moderately expressed. These genes are also seemed to be also expressed in control-stimulated condition, which has some similarities with FOXO1 non-stimulated. The deferentially expressed genes in stimulated are differently between other conditions between control and FOXO1. While, stimulated rest, shows similar expression in stimulated.

## GO Term Analysis

GO term analysis helps us understand the biological processes, molecular functions, and cellular components associated with the differentially expressed genes.

```{r go_analysis, fig.width=10, fig.height=8}
run_go_analysis <- function(results_df, title, ont) {
  diff_genes <- results_df$gene[results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1]
  
  ego <- enrichGO(gene = diff_genes,
                  OrgDb = org.Hs.eg.db,  # Human database
                  keyType = "SYMBOL",
                  ont = ont,  # Ontology type (BP, MF, CC)
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    print(dotplot(ego, showCategory = 20, title = paste("GO Enrichment Analysis -", title, "-", ont)))
  } else {
    print(paste("No significant GO terms found for", title, "-", ont))
  }
}

# Run GO analysis for Biological Process (BP)
run_go_analysis(results_non_stimulated, "Non-stimulated", "BP")
run_go_analysis(results_stimulated, "Stimulated", "BP")
run_go_analysis(results_rest, "Stimulated Rest", "BP")
```
![GO_BP_nonstimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/a2d50ffa-1dea-4152-8ed7-d1b2fbeb6661)
![GO_BP_stimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/fbba8e06-bc53-499c-b83a-188a16308950)
![GO_BP_rest](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/22800a5f-4518-483e-aabb-d568e8fe57c8)


The Biological Processes in the top differential genes are associated mostly with regulation of cell-cell adhesion; positive regulation of hydrolase activity in non-stimulated condition. In stimulated, it affects mostly on genes related to positive regulation of cytokine activity, while in rest to mononuclear differentiation.

```{r go_analysis, fig.width=10, fig.height=8}
# Run GO analysis for Molecular Function (MF)
run_go_analysis(results_non_stimulated, "Non-stimulated", "MF")
run_go_analysis(results_stimulated, "Stimulated", "MF")
run_go_analysis(results_rest, "Stimulated Rest", "MF")
```
![GO_MF_nonstimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/b5ad2c0e-f084-40ae-b347-95752a88024b)
![GO_MF_stimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/c647a399-2deb-4f8e-8114-6322d5670b30)
![GO_MF_rest](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/7b38bf2e-c354-4082-86a4-f0d1b4c3d8db)

The Molecular Function of the top deferentially expressed genes are associated with cytokine receptor binding and cytokine activity in both non-stimulated and stimulated rest, while in stimulated, genes have function of immune receptor activity and MHC complex binding.

```{r go_analysis, fig.width=10, fig.height=8}
# Run GO analysis for Cellular Component (CC)
run_go_analysis(results_non_stimulated, "Non-stimulated", "CC")
run_go_analysis(results_stimulated, "Stimulated", "CC")
run_go_analysis(results_rest, "Stimulated Rest", "CC")
```
![GO_CC_nonstimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/def5c3b8-4781-4b69-88fd-0e9e1bf8321a)
![GO_CC_stimulated](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/d04c918a-e415-477e-98b8-bfa7e3eaf9aa)
![GO_CC_rest](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/9949d944-62ab-467a-8ffb-7044fca52898)


The majority of deferentially expressed genes are associated to the cell leading edge, lamellipodium in non-stimulated and rest, while in stimulated they are located in different organelle membranes such as lysosomal, lytic vacuolar, late endosome and others.

## KEGG Pathway Analysis

KEGG pathway analysis provides insights into the molecular interaction and reaction networks associated with the differentially expressed genes.

```{r kegg_analysis, fig.width=10, fig.height=8}
run_kegg_analysis <- function(results_df, title) {
  diff_genes <- results_df$gene[results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1]
  
  gene_df <- bitr(diff_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)  # Changed to human database
  
  if (nrow(gene_df) > 0) {
    kegg_result <- enrichKEGG(gene = gene_df$ENTREZID,
                              organism = 'hsa',  # Changed to human
                              pvalueCutoff = 0.05)
    
    if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
      print(dotplot(kegg_result, showCategory = 20, title = paste("KEGG Pathway Enrichment -", title)))
    } else {
      print(paste("No significant KEGG pathways found for", title))
    }
  } else {
    print(paste("No genes could be mapped to Entrez IDs for KEGG analysis -", title))
  }
}

run_kegg_analysis(results_non_stimulated, "Non-stimulated")
run_kegg_analysis(results_stimulated, "Stimulated")
run_kegg_analysis(results_rest, "Stimulated Rest")
```
![KEGG_nonstim](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/7daa429f-11e1-4bb2-babb-d266664f826b)
![KEGG_STIM](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/b5086b35-5c5e-4f0c-8344-8b6b73b25baf)
![KEGG_rest](https://github.com/anelorda/FOXO1-differential-analysis/assets/98325381/576fb79c-b195-4ddf-a547-b12008774fc6)


KEGG Pathways show the molecular interactions of the most expressed genes in three conditions between control and FOXO1. As it can be seen from the graphs, the most differentiated genes are mostly associated with cytokine-cytokine receptor interaction and chemokine signaling pathway in both non-stimulated and stimulated rest conditions, while in stimulated the mostly differentiated genes are associated with Th17 cell differentiation and interactions related to the specific pathogens or diseases.

## Discussion

Based on the RNA-seq analysis of FOXO1 stimulation in human anti-Lewis Y CAR T cells, findings highlight distinct gene expression patterns across various experimental conditions: non-stimulated, stimulated, and stimulated rest. Principal Component Analysis (PCA) revealed significant separation primarily along PC1, driven by stimulation, indicating FOXO1's substantial influence on gene expression dynamics. Heatmap analysis further illustrated condition-specific expression patterns, where non-stimulated FOXO1 samples exhibited lower gene expression compared to controls, whereas stimulated conditions showed comparable expression levels with minor variations. Functional enrichment analyses, including GO and KEGG pathway analyses, identified enriched terms and pathways associated with immune modulation and response. Notably, in non-stimulated conditions, terms such as regulation of cell-cell adhesion and cytokine signaling pathways underscored FOXO1's role in basal immune surveillance. During stimulation, pathways like Th17 cell differentiation and cytokine production regulation highlighted FOXO1's involvement in immune activation. Post-stimulation, pathways related to immune regulation persisted, suggesting FOXO1's role in maintaining immune homeostasis. These findings collectively underscore FOXO1's pivotal role in immune function regulation, with implications for understanding immune-related diseases and therapeutic strategies aimed at modulating FOXO1 activity for clinical benefit. '

``` markdown
In this `Rmd` file:

1. **Introduction** and **Data Description** sections introduce the study and dataset.
2. **Data Preparation and Quality Control** section includes shell commands to download and unzip the data, as well as R code to load and prepare the data.
3. **Differential Expression Analysis with Preserved Gene IDs** uses DESeq2 to analyze differential expression.
4. **Volcano Plots** and **Principal Component Analysis (PCA)** sections create and display relevant plots.
5. **Heatmap of Top Differentially Expressed Genes** creates heatmaps for the top differentially expressed genes.
6. **GO Term Analysis** and **KEGG Pathway Analysis** sections perform functional enrichment analyses.
7. **Discussion** section is intended for summarizing findings and discussing the experiment.

This file provides a comprehensive workflow for RNA-seq analysis, from data preparation to visualization and enrichment analysis.
```

