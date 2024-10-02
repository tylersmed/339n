# Big public databases:
# 1. PubMed
# 2. Biomart
# 3. UniProt

# Both of these can be accessed through R packages.

# biomaRt package
renv::install('BiocManager')
renv::install("Bioc::biomaRt")
BiocManager::install('biomaRt')
library('biomaRt')
ensembl = useEnsembl(biomart = 'ensembl',
                     dataset = 'hsapiens_gene_ensembl')

genes = c("BRCA1", "TP53", "EGFR")
gene_info = getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = genes,
  mart = ensembl
)
print(gene_info)
