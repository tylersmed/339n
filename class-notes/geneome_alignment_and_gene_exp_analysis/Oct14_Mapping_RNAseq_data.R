## 10/14/24
# In this class, we will practice the following: 
#1- Plotting quality score for FASTQ files
#2- Alignment of raw reads using RSubread package
#3- From Alignment to Gene Counts using featureCounts
#4- Downloading published counts data from NCBI SRA using recounts package for differential analysis

##Bioconductor
#https://www.bioconductor.org/
#Advantages: open source software,  genomic data analysis tools, integration with R
#Biological data structures, example: "Granges" object handle genomic ranges
#Annotation Resources: genes, gene sets, genomic coordinates
#Datacamp access to Introduction to Bioconductor in R Course

setwd("~/Documents/BCH339N_teaching/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

# Install the 'ShortRead' package from Bioconductor.
# This package provides functionality to read, process, and output next-generation sequencing data.
#https://bioconductor.org/packages/release/bioc/html/ShortRead.html
BiocManager::install("ShortRead")
BiocManager::install("ggplot2")

library (reshape)
library (ggplot2)
library (ShortRead)




# Read the file
sread <- readFastq("subsampled.fastq")

# Plot quality
# Convert the quality scores of reads in 'sread' to a matrix format.
# Each row of this matrix represents a read, and each column represents a base position.
# The matrix entries are the quality scores.
q <- as(quality(sread), "matrix")
# Construct a data frame 'df' for plotting.
# - 'base' indicates the base position.
# - 'read' indicates the read number.
# - 'quality' indicates the quality score of the base at that position in that read.
df <- data.frame(base = factor(rep(1:ncol(q), each=nrow(q))),
                 read = factor(rep(1:nrow(q), ncol(q))),
                 quality = as.vector(q))
# Create a boxplot visualizing the distribution of quality scores for each base position across all reads.
# The x-axis represents base positions, the y-axis represents quality scores, 
# and the boxplots provide a summary of the quality scores' distributions at each position.
quality_plot= ggplot(df, aes(x = base, y = quality)) +
  geom_boxplot() +
  labs(title = "Quality Boxplots by Base Position", y = "Quality Score", x = "Base Position")

# --------------------------------------------------------------
# RETRIEVE GENOMIC DATA & METADATA
# --------------------------------------------------------------

# Link to the NCBI page for the chromosome 22 reference genome in the GenBank format.
# https://www.ncbi.nlm.nih.gov/nuccore/NC_000022.11?report=genbank

# Use NCBI eutils to download the sequence for chromosome 22 in FASTA format.
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000022.11&rettype=fasta&retmode=text

# Reference to a sample on NCBI's GEO (Gene Expression Omnibus) database.
# GEO Sample: GSM1609442

# Use sratools to fetch a sample from the SRA database.
# You can find sratools at https://github.com/ncbi/sra-tools
# The following command downloads 500,000 reads from the sample SRR1803211 in a compressed format.
# ~/sratoolkit.2.9.6-1-mac64/bin/fastq-dump --split-files --gzip -X 500000 SRR1803211

# --------------------------------------------------------------
# BUILD INDEX AND PERFORM ALIGNMENT
# --------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")
library("Rsubread")
# Building an index for the reference genome.
buildindex(basename="chr22",reference="chr22.fasta")

# Perform alignment on the sample data against the reference genome index.
align(index="chr22", readfile1="SRR1803211_1.fastq.gz", readfile2="SRR1803211_2.fastq.gz", type="rna",
      output_file="chr22.bam", minFragLength=50, maxFragLength=600)

align(index="chr22", readfile1="SRR1803211_1.fastq.gz", readfile2="SRR1803211_2.fastq.gz", type="rna",
      output_file="chr22.sam", minFragLength=50, maxFragLength=600, output_format = "SAM", unique = T)

# For an explanation on SAM/BAM flags, see:
# https://broadinstitute.github.io/picard/explain-flags.html

# To understand the function 'subjunc', view its documentation.
?subjunc

# --------------------------------------------------------------
# ANNOTATION AND FEATURE COUNTING
# --------------------------------------------------------------

# Count features in the aligned BAM file.
?featureCounts
my_counts = featureCounts(files="chr22.bam",isPairedEnd=TRUE,requireBothEndsMapped=TRUE, annot.inbuilt="hg38")

# Link to NCBI Gene database.
# https://www.ncbi.nlm.nih.gov/gene

# Fetch the Gencode annotation for the human genome (hg38).
annotation_human = getInBuiltAnnotation("hg38")

# If there are any installation issues, the following can be tried.
BiocManager::install("nloptr", update = TRUE, ask = FALSE)

# Install and load the 'recount' package from Bioconductor for accessing processed datasets.
BiocManager::install('recount')
library('recount')
browseVignettes('recount')

# Fetch and load a study from the 'recount' database.
study_id = "SRP055009"
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65912
browse_study(study_id)
url <- download_study(study_id, download = FALSE)
url
url2 <- download_study(study_id)
load(file.path(study_id, 'rse_gene.Rdata'))

# Extract and view specific data from the study.
gm19240_polyA_rep1_idx = rse_gene$sample == "SRS845361"
rse_gene$mapped_read_count[gm19240_polyA_rep1_idx]
rowData(rse_gene)
colData(rse_gene)
assays(rse_gene)$counts