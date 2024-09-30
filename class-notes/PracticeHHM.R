
# States and observables
states = c("sunny", "rainy")
observables = c("paint", "clean", "shop", "bike")

# Initial  and transition probabilities
initProbs = c(0.6, 0.4)
transProbs = matrix(c(0.8, 0.4, 0.2, 0.6), nrow=2, byrow=TRUE)
# Emission probabilities
emissionProbs = matrix(c(0.4, 0.1, 0.2, 0.3, 0.3, 0.45, 0.2, 0.05), ncol=4, byrow=TRUE)

# Observations
obsSeq = c("paint", "clean", "shop", "bike")
obsMapped = list('paint'=1, 'clean'=2, 'shop'=3, 'bike'=4)

fwdProbTbl = function(obsSeq) {
  fwdProbMatrix = matrix(0, nrow=length(states), ncol=length(obsSeq))
  # Initialization
  for (i in 1:length(states)) {
    fwdProbMatrix[i, 1] = initProbs[i] * emissionProbs[i, obsMapped[[obsSeq[1]]]]
  }
  # Fill out the table
  for (t in 2:length(obsSeq)) {
    for (s in 1:length(states)) {
      fwdProbMatrix[s, t] = sum(fwdProbMatrix[, t-1] * transProbs[, s]) * emissionProbs[s, obsMapped[[obsSeq[t]]]]
    }
  }
  return(list(fwdProbMatrix, sum(fwdProbMatrix[, length(obsSeq)])))
}

fwdProb = fwdProbTbl(obsSeq)
fwdProb

# Viterbi algorithm
viterbiManual = function(obsSeq) {
  R = length(obsSeq)
  N = length(states)
  
  # create matricies
  vProbMatrix = matrix(0, nrow=N, ncol=R)
  vPathMatrix = matrix(0, nrow=N, ncol=R)
  
  # Initialization
  vProbMatrix[, 1] = initProbs * emissionProbs[, obsMapped[[obsSeq[1]]]]
  
  for (t in 2:R) {
    for (s in 1:N) {
      tmp = vProbMatrix[, t-1] * transProbs[, s]
      vProbMatrix[s, t] = max(tmp) * emissionProbs[s, obsMapped[[obsSeq[t]]]]
      vPathMatrix[s, t] = which.max(tmp)
    }
  }
  optimalPath  = numeric(R)
  optimalPath[R] = which.max(vProbMatrix[, R])
  
  for (t in (R-1):1) {
    optimalPath[t] = vPathMatrix[optimalPath[t+1], t+1]    
  }
  
  names(optimalPath) = 1:R
  return(list(states[optimalPath], vProbMatrix, vPathMatrix))
}
  
viterbiManual(obsSeq)

# Using Libraries
#install.packages("HMM")
library(HMM)

weatherHHM = HHM.init(states, observables, initProbs, transProbs, emissionProbs)
forward(weatherHHM, obsSeq)
