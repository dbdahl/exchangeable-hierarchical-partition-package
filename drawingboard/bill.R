# p_nk: number of ways to partition n items into exactly k nonempty, unordered subsets
# according to the recurrence p(n,k) = p(n-1,k-1) + p(n-k,k)
p_nk <- function(n, k) {
  # handle edge cases
  if (n < 0 || k < 0)       return(0)
  if (n == 0 && k == 0)     return(1)
  if (n == 0 || k == 0)     return(0)   # (except the 0,0 case)
  if (k > n)                return(0)
  
  # dp[i+1, j+1] will hold p(i, j)
  dp <- matrix(0L, n + 1, k + 1)
  dp[1, 1] <- 1L   # p(0,0) = 1
  
  for (i in 1:n) {
    for (j in 1:min(i, k)) { 
      # recurrence:
      #   p(i, j) = p(i-1, j-1) + p(i-j, j)
      dp[i + 1, j + 1] <- dp[i, j] +
                          dp[i - j + 1, j + 1]
    }
  }
  
  dp[n + 1, k + 1]
}

# Examples
p_nk(5, 2)  # number of partitions of 5 items into exactly 2 subsets: should be 15
p_nk(7, 3)  # partitions of 7 items into 3 subsets

