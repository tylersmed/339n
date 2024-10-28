## We learned about NGS and RNA-Seq in class.
## In this assignment, we will analyze RNA-Seq data using DESeq2. 
## If you want to learn more about this package, check out the vignette
## https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
## We are going to carry out exploratory and differential expression analyses. 

#BiocManager::install("DESeq2")
library('DESeq2')

## In this homework, we will use the fread function from the data.table package. 
## fread allows us to extract .csv.gz files directly from the web. 
## If you are interested in learning more, data.table is a very popular R package.
#install.packages("R.utils")
library('R.utils')
#install.packages('data.table')
library('data.table')

## We are going to analyze the following dataset
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65912
## Note that there is a link to .csv.gz file with RNA Expression on this page.
## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65912/suppl/GSE65912_RNASeq.csv.gz

rnaseq_counts = fread ( 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65912/suppl/GSE65912_RNASeq.csv.gz', 
                        data.table = F)
row.names(rnaseq_counts) = rnaseq_counts$V1
rnaseq_counts = as.matrix (rnaseq_counts[,-1] ) 
as.matrix(rnaseq_counts)

## Q1 - 3 pt 
## Explore the data. 
## Note that we assigned gene names to the row.names attribute above. 
## The column headers correspond to samples and look like: "GM12878_Rep1_Counts"
## GM12878 corresponds to the cell line 
## Rep1 corresponds to the replicate number

## Use quantile function to determine how many total reads are assigned to genes 
quantile(rowSums(rnaseq_counts), probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99))

## A typical exploratory analysis is to inspect a scatterplot of read counts
## Plot two scatterplots using log10 of the gene read counts. 
## For the first plot compare GM12878_Rep1_Counts to GM12878_Rep2_Counts
## For the second plot compare GM12878_Rep1_Counts to GM12891_Rep1_Counts
## Briefly describe your observations regarding these two plots in a few sentences. 
plot(log10(rnaseq_counts[, "GM12878_Rep1_Counts"]), log10(rnaseq_counts[, "GM12878_Rep2_Counts"]))
plot(log10(rnaseq_counts[, "GM12878_Rep1_Counts"]), log10(rnaseq_counts[, "GM12891_Rep1_Counts"]))
# There is a positive correlation between both sets of variables. However, the plot shows a stronger
# correlation between GM12878_Rep1_Counts and GM12878_Rep2_Counts then GM12878_Rep1_Counts and GM12891_Rep1_Counts


## Q2 - 3 pt
## Use cor function to calculate all pairwise Spearman correlations between all the samples. 
## Draw a heatmap of the calculated pairwise Spearman correlations.
## To do this, you can use "pheatmap" package.
## Write a few lines describing the pattern(s) you observe in the heatmap

#install.packages('pheatmap')
library('pheatmap')
?pheatmap
?cor

cor_matrix = cor(rnaseq_counts, method = 'spearman')
pheatmap(cor_matrix)
# The strongest correlations exist between the repeats of the same gene. Here these are all in groups of 3.
# The correlation matrix also shows a weakening correlation between genes as they move further apart.

## We will  extract the cell line names from the column names. 
## We will then define a factor with two levels of length equal to the number of columns of rnaseq_counts
##  the reference level of this factor will match the cell lines:
## c("GM19238", "GM19240", "GM12892", "GM12878" )
## The alternative level of the factor will match: c("GM19239", "GM12891" ) 

cell_line_ids = sapply (strsplit(colnames(rnaseq_counts), "_"), "[[" , 1 ) 
Factor_level_1 = c("GM19238", "GM19240", "GM12892", "GM12878" )
Factor_level_2 =  c("GM19239", "GM12891" )
factor_of_interest = cell_line_ids %in% Factor_level_1
factor_of_interest = as.factor(factor_of_interest)
factor_of_interest = relevel(factor_of_interest, ref = "TRUE")

## Use DESeqDataSetFromMatrix to prepare the data for differential expression analysis.
dds <- DESeqDataSetFromMatrix (countData = rnaseq_counts,
                               colData = DataFrame(Variable = factor_of_interest ),
                               design = ~ Variable)



##  We will then run DESeq to identify differentially expressed genes. Here, you will also need to install "apeglm" package. Please fill in coef 
BiocManager::install("apeglm")

## Q3 - 3  pt
## We will  run DESeq to identify differentially expressed genes
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef= "Variable_FALSE_vs_TRUE")
resLFC <- resLFC[order(resLFC$pvalue),]
summary(resLFC, alpha = 0.05)

## In the above code, what does lfcShrink achieve? Explain in a few sentences.
## Generate an MA-plot using resLFC

?DESeq2::plotMA

## Q4 -3  pt
## Next we will create an interactive html to explore our results
#BiocManager::install("ReportingTools")
library("ReportingTools")
des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
                         title = 'RNA-seq analysis of differential expression using DESeq2',
                         reportDirectory = "./report")
## This might take a while. 
publish(dds,des2Report, pvalueCutoff=0.05,
        annotation.db="org.Hs.eg", factor = factor_of_interest,
        reportDir="./reports")
finish(des2Report)

## If everything went well, you should be able to inspect an html file containing your results.
## If successful, include one boxplot generated from this report. 
## Even if the above code fails, you can see the most differentially expressed genes with the following code
resLFC
## The top five differentially expressed genes by adjusted p-value should be: 
## KDM5C, PSMA8, KDM6A, ZFX, SMC1A
## Search what each of these genes are. 
## Do you notice any shared feature(s)? If so, can you speculate what our factor_of_interest was? 

