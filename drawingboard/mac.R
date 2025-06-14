# 1) Numerically‐stable log‐sum‐exp of two numbers
logSumExp <- function(a, b) {
  if (is.infinite(a) && a < 0) return(b)
  if (is.infinite(b) && b < 0) return(a)
  m <- max(a, b)
  m + log(exp(a - m) + exp(b - m))
}

# 2) Safe log‐subtraction: log(x - y) given log(x)=lx, log(y)=ly, returns -Inf if y>=x
logSubtract <- function(lx, ly) {
  if (ly >= lx) return(-Inf)
  lx + log1p(-exp(ly - lx))
}

# 3) Compute log d_r(n,k) via the two‐term recurrence in log‐space
log_d_r_nk <- function(n, k, r) {
  # Precompute logs of the falling‐factorials (i-1)!/(i-r)! for i>=r
  ff <- rep(-Inf, n + 1)
  if (r == 1) {
    ff[] <- 0
  } else {
    for (i in r:n) {
      ff[i] <- sum(log((i-1):(i-r+1)))
    }
  }
  
  # Circular buffer to hold last 'r' rows of dp
  buf <- matrix(-Inf, nrow = r, ncol = k + 1)
  
  # dp_cur[j+1] = log d_r(i, j)
  dp_cur <- rep(-Inf, k + 1)
  dp_cur[1] <- 0  # log d_r(0,0) = 0
  buf[1, ] <- dp_cur
  
  for (i in 1:n) {
    dp_prev <- dp_cur
    dp_cur  <- rep(-Inf, k + 1)
    
    for (j in 1:k) {
      # term1: (i-1) * d_r(i-1, j)
      t1 <- log(i - 1) + dp_prev[j + 1]
      # term2: ((i-1)!/(i-r)!) * d_r(i-r, j-1)
      if (i >= r) {
        row_idx <- ((i - r) %% r) + 1
        t2 <- ff[i] + buf[row_idx, j]
      } else {
        t2 <- -Inf
      }
      dp_cur[j + 1] <- logSumExp(t1, t2)
    }
    
    buf[((i) %% r) + 1, ] <- dp_cur
  }
  
  dp_cur[k + 1]  # log d_r(n,k)
}

# 4) Main function: returns vector of log‐weights log w[r], r=1..floor(n/k)
min_size_log_weights <- function(n, k) {
  Rmax <- floor(n / k)
  # compute log d_r(n,k) for r = 1..Rmax+1
  logd <- sapply(1:(Rmax + 1), function(r) log_d_r_nk(n, k, r))
  # compute log w[r] = log(d_r - d_{r+1})
  sapply(1:Rmax, function(r) logSubtract(logd[r], logd[r + 1]))
}

# -- Example usage and sanity check --
n <- 25; k <- 4
lw <- min_size_log_weights(n, k)
# normalize to get probabilities
logZ <- Reduce(logSumExp, lw)
p <- exp(lw - logZ)
print(round(p, 4))  # should be c(0.5882, 0.2941, 0.1176)



library(gourd)
crp <- CRPPartition(n, 1.0)
x <- samplePartition(crp, 1000000)
w <- apply(x, 1, \(x) length(unique(x)) == k)
table(apply(x[w, ], 1, \(x) min(table(x)))) / sum(w)

