# Intro to matricies and dataframes

# matrix
numeric_v1 = c(6, 4)
numeric_v2 = c(2, 3, 11)
numeric_v3 = c(61, 42, 100)

# cbind/rbind allows us to put together vectors as columns or rows of a matrix
mtx1 = cbind(numeric_v1, numeric_v2, numeric_v3)
mtx2 = rbind(numeric_v1, numeric_v2, numeric_v3)

mtx1
mtx2

mtx1 + 5
mtx1 + mtx2
# indexing a matrix
mtx1[1:2, 3]

# Row of Column Sums or Means
mtxA = matrix(c(1, 2, 3, 4, 4, 6), nrow = 2, byrow = T)
mtxA
rowSums(mtxA)
rowMeans(mtxA)

# Dataframes can hold different types of vectors
mtx3 = matrix(ncol = 2, nrow = 4) # fills in matrix with NAs

num_v4 = c(2.4, 6, 12, -5)
char_v1 = c('a', 'b', 'c', 'd')
mtx3[,1] = num_v4
mtx3[,2] = char_v1
mtx3

# modulo operator
5%%3 # gives the remainder

test_mtrx = matrix(1:50, nrow = 5)
test_mtrx

mtrx_val_change = function(m1) {
  for (r in 1:nrow(m1)) {
    for (c in 1:ncol(m1)) {
      if ((r + c) %% 3 == 0) {
        m1[r, c] = m1[r, c] + 5
      }
    }
  }
  return(m1)
}

test_mtrx_changed = mtrx_val_change(test_mtrx)
test_mtrx_changed

## dataframe
df1 = data.frame(CharCol = char_v1, NumCol = num_v4)

gene_exp = data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
                      Liver = c(12.5, 9.8, 15.3, 8.7, 10.1),
                      Kidney = c(11.5, 10.8, 12.3, 9.7, 1),
                      Heart = c(12.7, 8.8, 15, 7, 10.8))
gene_exp
gene_exp$Gene
gene_exp$Kidney

# Calculate the avg exp level for each gene across all tissues
numeric_rows = sapply(gene_exp, class) == 'numeric'
gene_exp$Avg_Expression = rowMeans(gene_exp[numeric_rows])
gene_exp

# identify the gene with the highest expression in the liver
max_liver_gene = gene_exp$Gene[which.max(gene_exp$Liver)]
max_liver_gene

# Add a column to classify genes as "highly expressed"
# if any tissue has an expression > 15
gene_exp$Highly_expressed = apply(gene_exp[numeric_rows], 1, function(x) {any(x>15)})
gene_exp

dna_seqs = data.frame(
  Species = c("Species1", "Species2", "Species3"),
  Sequence = c("ATGAGATGCTG", "GACTAGGGATCT", "GACTAGATCTAGATCGATC")
)
# Create a function to calculate the GC content of each sequence
gc_content = function(seq) {
  
}
