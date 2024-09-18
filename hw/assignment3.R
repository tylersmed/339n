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
  directions = traceback_matrix(score_matrix, m+1, n+1)
  directions = gsub(" ", "", directions)
  print(directions)
  translated = translate_directions(directions, seq1, seq2)
  print(translated)
  
}

# Recursively trace back the matrix and return a string of the directions to take
traceback_matrix = function(score_matrix, i, j) {
  if (i==1 && j==1) {
    return()
  }
  current = score_matrix[i, j]
  if (i == 1) {
    # Only possible to come from the left
    return(paste(traceback_matrix(score_matrix, i, j-1), " L"))
  } else if (j == 1) {
    # Only possible to come from above
    return(paste(traceback_matrix(score_matrix, i-1, j), " D"))
  }
  
  # Possible directions: diagonal, left, up
  diag = score_matrix[i-1, j-1]
  left = score_matrix[i, j-1]
  up = score_matrix[i-1, j]

  if (current == diag-10 | current == diag+3) {
    return(paste(traceback_matrix(score_matrix, i-1, j-1), "X"))
  } else if (current == left-3) {
    return(paste(traceback_matrix(score_matrix, i, j-1), "L"))
  } else {
    return(paste(traceback_matrix(score_matrix, i-1, j), "D"))
  }
}

#translate directions to alignment
translate_directions = function(directions, seq1, seq2) {
  seq1_aligned = ""
  seq2_aligned = ""
  i = 1
  n = 1
  directions = unlist(strsplit(directions, ""))
  for (direction in directions) {
    if (direction == "X") {
      seq1_aligned = paste(seq1_aligned, substr(seq1, i, i), sep="")
      seq2_aligned = paste(seq2_aligned, substr(seq2, n, n), sep="")
      i = i + 1
      n = n + 1
    } else if (direction == "L") {
      seq1_aligned = paste(seq1_aligned, "_", sep="")
      seq2_aligned = paste(seq2_aligned, substr(seq2, n, n), sep="")
      n = n + 1
    } else if (direction == "D") {
      seq1_aligned = paste(seq1_aligned, substr(seq1, i, i), sep="")
      seq2_aligned = paste(seq2_aligned, "_", sep="")
      i = i + 1
    }
  }
  
  return(list(seq1_aligned, seq2_aligned))
}

global_align("ATTAGC","ATTCAGG")

