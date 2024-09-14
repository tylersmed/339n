## Please add your code in between the comments
## For the first two questions on this assignment, please only use base R and do not rely on any packages.

## Q1 - 4pt
## GC-content of a piece of DNA is defined as the percentage of Gs + Cs. 
## For example "ACCTGCA" has a 57.1% GC-content
## Write a function that will take a sequence 
## The function should return a GC ratio

get_gc_ratio = function(s) {
  split_s = unlist(strsplit(s, ""))
  gc_count = sum(split_s == "G" | split_s == "C")
  return(gc_count/length(split_s))
}

print(get_gc_ratio("ACCTGCA")) # 0.571

## Q2 - 4pt
## Every amino acid in a protein is encoded by three nucleotides. 
## Execute the following two lines to get a list of all codons and corresponding amino acids
codons = c('UUU','UUC','UUA','UUG','UCU','UCC','UCA','UCG','UAU','UAC','UAA','UAG','UGU','UGC','UGA','UGG','CUU','CUC','CUA','CUG','CCU','CCC','CCA','CCG','CAU','CAC','CAA','CAG','CGU','CGC','CGA','CGG','AUU','AUC','AUA','AUG','ACU','ACC','ACA','ACG','AAU','AAC','AAA','AAG','AGU','AGC','AGA','AGG','GUU','GUC','GUA','GUG','GCU','GCC','GCA','GCG','GAU','GAC','GAA','GAG','GGU','GGC','GGA','GGG')
amino_acids = c('F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G' )

## Write a function that will take a coding region sequence as input. You can assume the sequence is starting with AUG and is a multiple of three nucleotides. 
## The output should be the corresponding protein sequence.
translate_codons = function(s) {
  protein = ""
  for (i in seq(1, nchar(s), 3)) {
    codon = substr(s, i, i+2)
    amino_acid = amino_acids[codons == codon]
    if (amino_acid == '*') {break}
    protein = paste(protein, amino_acid, sep="")
  }
  return(protein)
}

print(translate_codons("AUGUUUUCUUAGUCU")) # MFS (stop codon UAG in middle)

## Q3 - 4pt
## R has powerful plotting options available. 
## You may have already used ggplot package before 
## even the base R language has many simple plotting functions
## Let's first install the dslabs package
install.packages("dslabs")
library(dslabs)
## This package has a number of interesting datasets
## We can load an example for simple plotting
data(gapminder)
## We can explore some of the features in this dataset with basic plots
?gapminder
## We can look at the first few entries of this dataset with "head"
head ( gapminder)
str (gapminder)
## Given this data frame let's extract the data from 2015 into a new data frame
gapminder_2015 = gapminder[gapminder$year == 2015 ,  ]

### plot the life expectancy (x axis) and infant mortality (y axis) using the data from 2015(2 pts).  
### make a histogram of distribution of Fertility using the data from 2015 (2pts).
plot(gapminder_2015$life_expectancy, 
     gapminder_2015$infant_mortality,
     main = "Infant Mortality vs Life Expectancy in 2015",
     ylab = "Infant Mortality",
     xlab = "Life Expectancy")

hist(gapminder_2015$fertility,
     main = "Fertility Rates in 2015",
     xlab = "Fertility Rate",
     ylab = "Number of Countries",
     breaks = 10)

