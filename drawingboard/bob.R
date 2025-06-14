# Function to compute p(S | K=k) for CRP, in base R
p_crp_sizes_given_k <- function(sizes, k = length(sizes)) {
  # Ensure sizes sum to n:
  n <- sum(sizes)
  stopifnot(k == length(sizes), k <= n)
  # Compute |s(n,k)| via dynamic programming (Stirling numbers of first kind)
  # We use a vector where dp[j] = |s(i, j)| for current i.
  dp <- numeric(k+1)
  dp[1] <- 1            # |s(0,0)| = 1  (dp index 1 represents 0 clusters)
  for (i in 1:n) {
    # update for each new element i
    # Use next_dp for |s(i, j)| based on previous |s(i-1, j)|
    next_dp <- numeric(k+1)
    # j = 0 case implicitly 0 (no way to have 0 cycles if i>0)
    for (j in 1:min(i, k)) {
      next_dp[j+1] <- dp[j] + (i-1) * dp[j+1]
    }
    dp <- next_dp
  }
  c_nk <- dp[k+1]       # this is |s(n,k)|
  
  # Compute N(S) = n! / (prod_{j} s_j! * prod_{u} m_u!)
  # and multiply by prod (s_j - 1)! for the numerator count
  sizes <- as.integer(sizes)
  freq  <- as.integer(table(sizes))      # frequencies of each size
  # Compute numerator in a stable way:
  # We'll accumulate in numerator and denominator to avoid overflow.
  num <- 1        # will hold n! * prod(s_j - 1)!
  den <- 1        # will hold prod(s_j!) * prod(m_u!)
  # Multiply num by n! and den by prod(s_j!) in parts
  for (s in 1:n) {
    num <- num * s            # multiply by all factors up to n (forming n!)
    # If s equals any cluster size, incorporate that cluster's factorial in denom:
    # (We do this incrementally to cancel factors early if possible)
    count_index <- which(names(freq) == as.character(s))
    if (length(count_index) > 0) {
      # multiply denominator by (s!)^freq if cluster size s appears freq times
      # But doing s! directly would be large; instead multiply by s stepwise.
      # Since we are in the loop for n!, when s == cluster size value, we should 
      # multiply denominator by s (for each occurrence of cluster of that size).
      # Actually, better: accumulate cluster factorial fully after loop to avoid confusion.
    }
  }
  # (For simplicity, we'll compute numerator and denominator fully then divide, 
  # relying on R's big integers if n is not too large. For larger n, a prime factor or 
  # incremental cancellation approach would be needed.)
  num <- factorial(n)                        # n!
  den <- 1
  for (sz in sizes) den <- den * factorial(sz)   # prod(s_j!)
  for (f in freq)  den <- den * factorial(f)     # prod(m_u!)
  # Now numerator count of permutations for S:
  numerator_count <- num / den * prod(factorial(sizes - 1))
  
  # Finally, p(S | K=k) = numerator_count / |s(n,k)|
  p_val <- numerator_count / c_nk
  return(p_val)
}

all_possible <- list(c(7,1,1), c(6,2,1), c(5,2,2), c(5,3,1), c(4,3,2), c(4,4,1), c(3,3,3))
sum(sapply(all_possible, \(x) p_crp_sizes_given_k(x)))

















# ---- CRP Ordered PMF ----
crp_ordered_pmf <- function(sizes, alpha) {
  n <- sum(sizes)
  k <- length(sizes)
  m <- table(sizes)

  # Rising factorial: alpha^{(n)} = alpha * (alpha+1) * ... * (alpha+n-1)
  alpha_rising_n <- prod(alpha + 0:(n - 1))

  num <- factorial(n) * alpha^k
  denom <- alpha_rising_n * prod(factorial(as.integer(m))) * prod(sizes)

  return(num / denom)
}

# ---- Log-space Stirling number of the second kind (exact inclusion-exclusion) ----
log_stirling2 <- function(n, k) {
  s <- 0
  for (j in 0:k) {
    sign <- (-1)^(k - j)
    term <- choose(k, j) * j^n
    s <- s + sign * term
  }
  return(log(s / factorial(k)))
}

# ---- Log-space marginal probability of k clusters ----
library(gmp)
crp_p_k <- function(n, k, alpha = 1.0) {
  log_rising_factorial <- function(alpha, n) {
    lgamma(alpha + n) - lgamma(alpha)
  }

  log_stirling2 <- function(n, k) {
    log(as.numeric(Stirling2(n, k)))
  }

  log_p <- k * log(alpha) + log_stirling2(n, k) - log_rising_factorial(alpha, n)
  return(exp(log_p))
}

# ---- Final: Conditional PMF given k clusters ----
crp_s_given_k <- function(sizes, alpha = 1.0) {
  n <- sum(sizes)
  k <- length(sizes)

  pS <- crp_ordered_pmf(sizes, alpha)
  pk <- crp_p_k(n, k, alpha)

  return(pS / pk)
}

all_possible <- list(c(7,1,1), c(6,2,1), c(5,2,2), c(5,3,1), c(4,3,2), c(4,4,1), c(3,3,3))
sum(sapply(all_possible, \(x) crp_s_given_k(x, alpha = 1.0)))

