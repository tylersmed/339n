# Fibonacci numbers (recursion, memoization)

# Fibonaci numbers: 1, 1, 2, 3, 5, 8, 13, 21, ...
# F(n) = F(n-1) + F(n-2)

# Write a function that calculates the n-th Fibonacci number iteratively
fib_iter = function(n) {
  if(n == 1 | n == 2) {
    return(1)
  }
  else {
    fibo_prev = 1
    fibo_current = 1
    for (i in 3:n) {
      fibo_next = fibo_prev + fibo_current
      fibo_prev = fibo_current
      fibo_current = fibo_next
    }
  }
  return(fibo_current)
}

print(fib_iter(5))
print(fib_iter(8))

# Write a function that calculates the n-th Fibonacci number recursively
fib_rec = function(n) {
  if (n==1 | n==2) {
    return(1)
  }
  else {
    return(fib_rec(n-1) + fib_rec(n-2))
  }
}

print(fib_rec(5))
print(fib_rec(8))

# Write a function that calculates the n-th Fibonacci number recursively with memoization

fibo_nums = numeric(0)
fibo_nums[1] = 1
fibo_nums[2] = 1
fib_memo = function(n) {
  if (is.na(fibo_nums[n])) {
    # store result if the number is not already calculated
    fibo_nums[n] = fib_memo(n-1) + fib_memo(n-2)
  }
  return(fibo_nums[n])
}

print(fib_memo(5))
print(fib_memo(8))

# Sequence Alignments
# Recursive approach to compute the number of alignments between sequences of length m and n
num_alignments = function (n, m) {
  if (n==0) {return(1)}
  if (m==0) {return(1)}
  return(num_alignments(n-1, m) + num_alignments(n, m-1) + num_alignments(n-1, m-1))
}

print(num_alignments(2, 2))
print(num_alignments(3, 3))

# Combinatorial formula
closed_formula = function(m, n) {
  total_alignments = 0
  for (k in 0:min(n,m)) {
    total_alignments = total_alignments + 2^k * choose(m, k) * choose(n, k)
  }
  return(total_alignments)
}
closed_formula(3, 3)
