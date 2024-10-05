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
print(paste("Rows:", dim(df_cleaned)[1], "Columns:", dim(df_cleaned)[2]))

##Data Splitting: Split the cleaned dataset into two parts: The training set should contain 10% of the sequences.
##The test set should contain the remaining 90% to assess the performance of your HMM in predicting secondary structure.

## To ensure reproducibility, use the set.seed function with a seed value of 3:

?set.seed
set.seed(3)

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