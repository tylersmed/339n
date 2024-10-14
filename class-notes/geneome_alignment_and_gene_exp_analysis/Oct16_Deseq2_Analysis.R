#10/16/24
#We are analyzing the gene expression profiles of individuals 
#carrying the Diamond Blackfan Anemia mutation.
#Two individuals are more severely affected than the other carriers.
#Our goal is to determine if there is a specific gene expression pattern 
#associated with carrying the mutation and being more affected.
#We will cluster and visualize the data, and use the DESeq2 package 
#to analyze gene expression patterns. 


# Load required packages
install.packages("BiocManager")
install.packages("tidytable")
install.packages("tibble")
BiocManager::install("DESeq2")
library(tidytable)
library(tibble)
library(DESeq2)
library(RColorBrewer)
install.packages("pheatmap")
library(pheatmap)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggplot2")
library(ggplot2)
library(dplyr)


# Reading the raw data
DBA_data <- read.table("~/Downloads/unprocessed_counts.tab")
DBA_raw_counts <- DBA_data[, c(1, 22, 11, 16, 23, 5, 36, 17, 14, 38, 15, 42)]
colnames(DBA_raw_counts) <- c(
  "carrier1_s1", "carrier1_s2", "carrier2_s1", "carrier2_s2",
  "noncarrier1_s1", "noncarrier1_s2", "noncarrier2_s1", "noncarrier2_s2",
  "sick1_s1", "sick1_s2", "sick2_s1", "sick2_s2"
)

# Plotting the raw count distribution as histograms
# You can customize these plots as needed
ggplot(DBA_raw_counts) + geom_histogram(aes(x = carrier1_s1), bins = 200) +
  xlab("Raw expression counts") + ylab("Number of genes")
ggplot(DBA_raw_counts) + geom_histogram(aes(x = noncarrier1_s1), bins = 200) +
  xlab("Raw expression counts") + ylab("Number of genes") +
  xlim(0, 2000) + ylim(0, 2000)

ggplot(DBA_raw_counts, aes(x = carrier1_s1)) +
  geom_dotplot(binwidth = 10, method = "histodot") +  # Adjust binwidth as needed
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  ggtitle("Dot Plot of Raw Expression Counts for carrier1_s1")+
  xlim(0, 2000) + ylim(0, 2000)


vignette ("DESeq2")

# Preparing the metadata table
genotype <- c(
  "carrier", "carrier", "carrier", "carrier",
  "noncarrier", "noncarrier", "noncarrier", "noncarrier",
  "carrier", "carrier", "carrier", "carrier"
)
condition <- c(
  "normal", "normal", "normal", "normal",
  "normal", "normal", "normal", "normal",
  "DBA", "DBA", "DBA", "DBA"
)

DBA_metadata <- data.frame(genotype, condition)
rownames(DBA_metadata) <- c(
  "carrier1_s1", "carrier1_s2", "carrier2_s1", "carrier2_s2",
  "noncarrier1_s1", "noncarrier1_s2", "noncarrier2_s1", "noncarrier2_s2",
  "sick1_s1", "sick1_s2", "sick2_s1", "sick2_s2"
)

# Check if row names in DBA_metadata and column names in DBA raw data are in the same order
all(colnames(DBA_raw_counts) == rownames(DBA_metadata))

# Creating the DESeq2 object
dds_DBA <- DESeqDataSetFromMatrix(countData = DBA_raw_counts, colData = DBA_metadata, design = ~condition + genotype)

# Normalization of counts, why do we need this? 

dds_DBA <- estimateSizeFactors(dds_DBA)

# scaling raw counts ,by the size factors
# Normalized counts: extraction
normalized_DBA_counts <- counts(dds_DBA, normalized = TRUE)


# Unsupervised clustering analysis, log transformation
vsd_DBA <- vst(dds_DBA, blind = TRUE)
vsd_mat_DBA <- assay(vsd_DBA)

# Compute pairwise correlation values
vsd_cor_DBA <- cor(vsd_mat_DBA)
pheatmap(vsd_mat_DBA, 
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Heatmap of vsd_mat_DBA")
# Principal component Analysis (PCA): 
#PCA is a statistical technique used to transform high-dimensional data into 
#a lower-dimensional representation while preserving as much of the original 
#data's variability as possible
plotPCA(vsd_DBA, intgroup = "genotype")
plotPCA(vsd_DBA, intgroup = "condition")

# DESeq2 analysis
dds_DBA <- DESeq(dds_DBA)

# Calculating mean for each gene (each row)
mean_counts <- apply(DBA_raw_counts[, 1:12], 1, mean)

# Calculating variance for each gene (each row)
variance_counts <- apply(DBA_raw_counts[, 1:12], 1, var)

# Plotting relationship between mean and variance
df <- data.frame(mean_counts, variance_counts)
ggplot(df) +
  geom_point(aes(x = mean_counts, y = variance_counts)) +
  scale_y_log10() + scale_x_log10() +
  xlab("Mean Counts per gene") + ylab("Variance per gene")

#What does this mean? 
#the variance of measurements tends to increase 
#as the mean of the measurements increases. 

#Deseq2 handles mean-variance dependence in a couple of ways: 
#1-Normalization: scales raw counts using size factors to account 
#for sequencing depth differences.
#2-Modeling Mean-Variance: uses the negative binomial distribution 
#to model the mean-variance relationship.
#Variance Stabilization: applies variance stabilization (which rlog transformation)
#to make variance approximately constant across means.
#Differential Expression: identifies differentially expressed genes 
#while handling increased variance at higher means.

#Variance stabilization instapoll

#What is dispersion?
#variability in gene expression levels across samples or conditions.
# Plot dispersion estimates
plotDispEsts(dds_DBA)
#DESeq2estimates the dispersion for each gene as part of 
#its statistical modeling. This estimation is crucial 
#because it helps account for the inherent variability in RNA-Seq data. 

#Mean-dispersion plots are valuable for quality control and 
#understanding the characteristics of the data before performing 
#differential expression analysis. 

#Mean-Variance Relationship, Dispersion Estimation instapoll

# DESeq2 results
DBA_res <- results(dds_DBA, contrast = c("condition", "normal", "DBA"), alpha = 0.05)

# DESeq2 results table
mcols(DBA_res)
head(DBA_res, n = 10)

# Differentially expressed genes
summary(DBA_res)

# Log2 fold change threshold
DBA_res <- results(dds_DBA, contrast = c("condition", "normal", "DBA"), alpha = 0.05, lfcThreshold = 0.35)
summary(DBA_res)

# Creating a data frame with DESeq2 results
DBA_res_all <- data.frame(DBA_res)
View(DBA_res_all)


