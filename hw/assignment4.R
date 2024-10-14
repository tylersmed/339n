##The Viterbi and Forward algorithms discussed in class are implemented in the R package HMM.

##You can access this package with:

library(HMM)
?initHMM
?forward
?viterbi

##In this and the next homework assignment, we will develop a Hidden Markov Model (HMM) to predict the secondary structure of proteins.
##Your answers from this homework will be utilized in the next one.
##Kaggle (https://www.kaggle.com) hosts various data science and machine learning challenges, including several in the field of bioinformatics.
##We will focus on the Protein Secondary Structure Prediction challenge, which can be found here: ##https://www.kaggle.com/alfrandom/protein-secondary-structure/data.
##Please review the introduction to become familiar with this dataset.


## Q1 (4 points)

##Understanding Protein Secondary Structure

## What is protein secondary structure? Briefly explain.
##  - Protein secondary structure is the local spatial conformation of backbone of the protein
##    excluding the amino acid side chains

## We will be working with the three-state secondary structure encoding (sst3):
## E: Strand of a β-pleated sheet, H: Right-handed α-helix, C: Coil
## Explain the biochemical features of each of the three states (C, E, and H) that we will model.
##  - C (coils): These consist of the 'loops' in a protein structure. They are sections of the secondary
##               structure without a highly orgainzed, repeating structure.
##  - E (sheets): These are formed with hydrogen bonds.
##                In sheets, they consist of multiple strands of amino acids that alternate in conformation
##                and 'stack' on top of eachother.
##  - H (helix): Like sheets, these are formed with hydrogen bonding between the NH and carbonyl groups in the backbone.
##               However, the amino acids all adopt the same conformation and assemble into a repeating helical pattern.

##Q2 (4 points)

## For this homework, use the dataset with strict quality controls.
##Data Cleaning: Remove all sequences that contain non-standard amino acids. Make sure to carefully check all columns in the data file.
df = read.csv('hw3_dataset.csv')

symbols = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
symbols_list = paste(symbols, collapse = "")

df_cleaned = df[grepl(paste0("^[" , symbols_list, "]*$"), df$seq),]

## Print the number of rows and columns remaining after cleaning.
print(paste("Rows:", nrow(df_cleaned), "Columns:", ncol(df_cleaned)))

##Data Splitting: Split the cleaned dataset into two parts: The training set should contain 10% of the sequences.
##The test set should contain the remaining 90% to assess the performance of your HMM in predicting secondary structure.

## To ensure reproducibility, use the set.seed function with a seed value of 3:

?set.seed
set.seed(3)

sample = sample(c(TRUE, FALSE), nrow(df_cleaned), replace = TRUE, prob = c(0.1, 0.9))
train = df_cleaned[!sample,]
test = df_cleaned[sample,]

## Q3 (4 points)

## Inferring HMM Parameters:  Use the training data to infer the parameters of the HMM.
## Specifically, use the sst3 and seq columns of the dataset to determine: Initial probabilities, Transition probabilities, Emission probabilities
## Write a function that takes the training data as input and returns a list with the following named elements:

##"initial_probs"
##"emission_probs" 
##"transition_probs"

##You will use the following symbols and states:
symbols = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
states = c("C", "E", "H")

## Your function should output a list structured as:
## list ( initial_probs = ...,       # A vector of initial state probabilities
## emission_probs = ...,      # A matrix of emission probabilities
## transition_probs = ...     # A matrix of transition probabilities)

## Note: Make sure to properly clean, preprocess, and validate your data to ensure accurate results in both this and future assignments

# Function to calculate the initial, emission, and transition probabilities
get_hmm_params = function(train, symbols, states) {
  # Initialize the initial, emission, and transition probabilities
  initial_probs = c(C = 0, E = 0, H = 0)
  emission_probs = matrix(0, nrow = length(states), ncol = length(symbols), dimnames = list(states, symbols))
  transition_probs = matrix(0, nrow = length(states), ncol = length(states), dimnames = list(states, states))
  
  # Calculate the initial probabilities
  initial_state_counts = table(substr(train$sst3, 1, 1)) # get first char of sst3
  for (state in states) {
    if (state %in% names(initial_state_counts)) {
      initial_probs[state] = initial_state_counts[state] / nrow(train)
    }
  }

  # Calculate the emission probabilities
  for (i in 1:nrow(train)) { # go over each row and split the seq and sst3 strings
    seq = strsplit(as.character(train$seq[i]), "")[[1]]
    sst3 = strsplit(as.character(train$sst3[i]), "")[[1]]
    for (j in 1:length(seq)) { # go over each character in the seq and sst3 strings and add to the emission_probs matrix
      emission_probs[sst3[j], seq[j]] = emission_probs[sst3[j], seq[j]] + 1
    }
  }
  emission_probs = emission_probs / rowSums(emission_probs) # convert total counts to probabilities
  
  # Calculate the transition probabilities
  for (i in 1:nrow(train)) { # go over each row and split the sst3 string
    sst3 = strsplit(as.character(train$sst3[i]), "")[[1]]
    for (j in 2:length(sst3)) { # go over each character in the sst3 string and add to the transition_probs matrix
      transition_probs[sst3[j - 1], sst3[j]] = transition_probs[sst3[j - 1], sst3[j]] + 1
    }
  }
  transition_probs = transition_probs / rowSums(transition_probs) # convert total counts to probabilities
  
  return(list(initial_probs = initial_probs, emission_probs = emission_probs, transition_probs = transition_probs))
  
}
hmm_params = get_hmm_params(train, symbols = symbols, states = states)
print(hmm_params)



##############################################################
##                         Part 2                           ##
##############################################################



##  In this assignment, we will apply the Viterbi and Forward Algorithms as discussed in class.
##You will be working with the R package HMM to implement these concepts.
## Please note that you will need to use the emission and transition probability tables that were generated in Homework 5.
##If you missed that homework, please attend Yue's TA session on Monday for assistance.
##Important: Ensure that your emission_probs and transition_probs matrices are structured correctly and are compatible with the initHMM function in the HMM package.

## Q1 (4 points) Building and Evaluating the HMM
## Review the Documentation: Refer to the documentation for the HMM package in R.
## You can use ?HMM and ?initHMM for reference.
## Using the emission and transition probabilities inferred in Homework 5, create an HMM using the initHMM function

?initHMM
my_hmm = initHMM(states, symbols, startProbs = hmm_params$initial_probs, transProbs = hmm_params$transition_probs, emissionProbs = hmm_params$emission_probs)

##  For each sequence in the test dataset (use the seq column), use the HMM to predict the most likely path of hidden states.
for (i in 1:nrow(test)) {
  seq = strsplit(test$seq[i], "")
  path = viterbi(my_hmm, seq[[1]])
  test$predicted_sst3[i] = paste(path, collapse = "")
}

## Compare this predicted path to the actual values from the sst3 column (i.e., the true secondary structure states).
print(head(test$predicted_sst3))
print(head(test$sst3))
## For each sequence in the test set, calculate the percentage of amino acids for which the predicted secondary structure matches the actual structure.  This will give you the prediction accuracy.
for (i in 1:nrow(test)) {
  test$accuracy[i] = sum(strsplit(test$predicted_sst3[i], "")[[1]] == strsplit(test$sst3[i], "")[[1]]) / nchar(test$sst3[i])
}

## Q2 (4 points) Analyzing Prediction Accuracy Distribution
##   Calculate the distribution of the prediction accuracies across all examples in the test set.

##Plot the distribution of these accuracy percentages.
hist(test$accuracy, main = "Prediction Accuracy Distribution", xlab = "Accuracy Percentage", ylab = "Frequency")
## Provide the following statistical measures for the distribution:Mean, Median, 90th percentile
print(paste("Mean:", mean(test$accuracy)))
print(paste("Median:", median(test$accuracy)))
print(paste("90th percentile:", quantile(test$accuracy, 0.9)))

## Q3 (4 points)  Low-Accuracy Predictions and Model Limitations
## Identify any test examples where the prediction accuracy is below 1%.
test$low_accuracy = test$accuracy < 0.05
low_accuracy_examples = test[which(test$low_accuracy == TRUE),]

## List the pdb_ids (Protein Data Bank identifiers) of these examples.
## Examine the predicted and actual secondary structure sequences for these examples.

## For each of these low-accuracy examples, calculate the percentage of amino acids labeled as "E" (strand of a β-pleated sheet) in the actual secondary structure.
## Plot the relationship between this "E" content percentage and the overall accuracy percentage of the prediction.
## Calculate the correlation coefficient to assess the relationship between "E" content and prediction accuracy.

## Provide insights into why the HMM might perform poorly for certain proteins.
## Consider the biological characteristics of the "E" state and how it may affect the model's predictions.
## Discuss any limitations or assumptions of the Markov chain model that might contribute to these prediction challenges.