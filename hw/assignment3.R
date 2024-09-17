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
  m = nchar(seq1)
  n = nchar(seq2)
  score_matrix = matrix(-Inf, nrow = m+1, ncol = n+1)
  # fill out the matrix
  for (i in 1:(m+1)) {
    for (j in 1:(n+1)) {
      if (i == 1) {
        score_matrix[i,j] = gap_open * (j-1)
      } else if (j == 1) {
        score_matrix[i, j] = gap_open * (i-1)
      } else {
        match = score_matrix[i-1, j-1] + ifelse(substr(seq1, i-1, i-1) == substr(seq2, j-1, j-1), match_score, mismatch_score)
        gap_right = score_matrix[i, j-1] + gap_open
        gap_down = score_matrix[i-1, j] + gap_open
        score_matrix[i, j] = max(match, gap_right, gap_down)
      }
    }
  }
  print(score_matrix)
  # trace the matrix back
  top_seq = traceback_matrix(score_matrix, seq1, m+1, n+1, T)
  bottom_seq = traceback_matrix(score_matrix, seq2, m+1, n+1, F)
  print(top_seq)
  print(bottom_seq)
  
}

traceback_matrix = function(score_matrix, seq, i, j, top_seq) {
  if (i==1 && j==1) {
    return('')
  }
  
  current = score_matrix[i, j]
  if (i == 1) {
    # Only possible to come from the left
    return(paste(traceback_matrix(score_matrix, seq, i, j-1, top_seq), "_"))
  } else if (j == 1) {
    # Only possible to come from above
    return(paste(traceback_matrix(score_matrix, seq, i-1, j, top_seq), "_"))
  }
  
  # Possible directions: diagonal, left, up
  diag = score_matrix[i-1, j-1]
  left = score_matrix[i, j-1]
  up = score_matrix[i-1, j]
  print(paste(current,diag, left, up))
  
  if (current == diag+3 | current == diag-10) {
    return(paste(traceback_matrix(score_matrix, seq, i-1, j-1, top_seq), substr(seq, i-1, i-1)))
  } else if (top_seq) {
    if (current == left+3) {
      return(paste(traceback_matrix(score_matrix, seq, i, j-1, top_seq), substr(seq, i-1, i-1)))
    } else {
      return(paste(traceback_matrix(score_matrix, seq, i, j-1, top_seq), "_"))
    } else if (current == up+3) {
      return(paste(traceback_matrix(score_matrix, seq, i-1, j, top_seq), substr(seq, j-1, j-1)))
    } else {
      return(paste(traceback_matrix(score_matrix, seq, i-1, j, top_seq), "_"))
    }
  }
}

global_align("ATTAGC","ATTCAGG")
