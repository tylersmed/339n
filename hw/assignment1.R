## Q1- 4pts
## Write a function that will simulate the number of tails in a series of coin tosses. 
## The input to the function is the number of coin tosses
## The output should be a single number equivalent to the number of tails in the series
number_of_tails = function(num_tosses) {
  coin = c('heads', 'tails')
  flips = sample(coin, size = num_tosses, replace = TRUE)
  num_tails = table(flips)['tails']
  return(num_tails)
}

## Q2 - 4pts
## Using the function you wrote in Q2, generate 5000 experiments each with 40 coin tosses
## Plot the distribution of the outputs as a histogram
hist(replicate(5000, number_of_tails(40)), 20, 
     main = "Distribution of Number of Tails",
     xlab = "Number of Tails",
     ylab = "Numbe of Trials") 

## Q3- 4pts
## Values and objects have different data types in R 
## You are given two DNA sequences of equal length
## find out how many bases are different between the two
dna1 = "TTAGCCTAAGCTA"
dna2 = "TTAGCCGATGCTA"
dna1_split = unlist(strsplit(dna1, ''))
dna2_split = unlist(strsplit(dna2, ''))
num_bases_diff = sum(dna1_split != dna2_split)
print(num_bases_diff)


