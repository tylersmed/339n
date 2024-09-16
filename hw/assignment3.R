## In class, we have discussed dynamic programming for globally aligning two sequences. 
## Please implement this algorithm. 
## Specifically, the function should take two strings as well as scores for match, mismatch, gap penalty
## Your algorithm should work with the two given strings and match, mismatch, gap penalty scores.
## Input example: needleman_wunsch_pairwise("ATTAGC","ATTCAGG", match_score=3, mismatch_score=-10, gap_open=-3)
## You can validate your algorithm using simple sequences by manually drawing a matrix and calculating the alignment score
## Note that all mismatches will have the same score for this question (unlike Blosum-50). 
## We will use a linear gap penalty (we will not use affine score: no separate penalty for gap open). 
## The output should be optimal alignment and the associated score

## Expected Output Example: 
## ATT-AGC
## ATTCAGG
## Score: 18

##The score 18 isn't necessarily right with the above given needleman_wunch parameters

global_align = function(seq1, seq2, match_score = 3, mismatch_score = -10, gap_open = -3) {
  score_matrix = matrix(-Inf, nrow = nchar(seq1) + 1, ncol = nchar(seq2) + 1)
  print(score_matrix)
}


global_align("ATTAGC","ATTCAGG")
