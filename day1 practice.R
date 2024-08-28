# Data types

# numeric
my_value = 3.3

# string/text
my_dna = "GTACCCTA"

is.numeric(my_value)
is.numeric(my_dna)

# boolean values --logical
str(my_dna)
str(my_value)

?str

# basic data types can be combined using 'c' into vectors
my_numeric_vertor = c(2, 4, 6, -5)
str(my_numeric_vertor)
my_char_vector = c('a', 'b', 'c', 'd')
str(my_char_vector)

# can use functions with vectors
my_numeric_vertor/2
sum(my_numeric_vertor)
my_logical_vector = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
sum(my_logical_vector)

?sample
# create random sample
sample(my_numeric_vertor, replace = T, size=5)

?replicate
# basically can do what pandas does
?strsplit
my_dna_split = strsplit(my_dna, split='')[[1]]
my_dna_pasted = paste(my_dna_split, collapse = '')

# functions and loops
my_vector_sum = function(my_vector) {
  s = 0
  for (i in 1:length(my_vector)) {
    s = s + i
  }
  return(s)
}

my_sum = my_vector_sum(my_numeric_vertor)
my_sum

# indexing a string
?substr
substr(my_dna, 2, 2)
substr(my_dna, 2, 4)

# indexing a vector. R starts indexing at 1
my_numeric_vertor[3]
my_numeric_vertor[2:4]

# %in% operator
# to check overlap between two vectors

vals_to_search = c(6, 24)
vals_to_search %in% my_numeric_vertor

my_function = function(s, ss) {
  indicies = c()
  for (i in 1:nchar(s)) {
    temp_ss = substr(s, i, i+nchar(ss)-1)
    if (ss == temp_ss) {
      indicies <- c(indicies, i)
    }
  }
  return(indicies)
}

x = my_function("GATATATGCATATACTT", "ATAT")
print(x)
