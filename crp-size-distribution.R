library(gourd)
library(gupd)


# Function to compute the unsigned Stirling numbers of the first kind recursively
stirling <- function(n, k) {
  # Base cases
  if(n == 0 && k == 0) return(1)
  if(n == 0 || k == 0) return(0)
  if(n == k) return(1)
  
  # Recursive step: S(n, k) = S(n-1, k-1) + (n-1) * S(n-1, k)
  return(stirling(n - 1, k - 1) + (n - 1) * stirling(n - 1, k))
}

# Recursive sampling function for cluster sizes using the corrected weights
sample_cluster_sizes <- function(n, k) {
  # Base case: one cluster must contain all items.
  if (k == 1) {
    return(n)
  }
  
  # The first cluster size s1 can range from 1 to n - k + 1,
  # because each of the remaining clusters must have at least one item.
  possible_s1 <- 1:(n - k + 1)
  
  # Compute weights:
  # Weight(s1) = choose(n-1, s1-1) * (s1-1)! * S(n-s1, k-1)
  weights <- sapply(possible_s1, function(s1) {
    choose(n - 1, s1 - 1) * factorial(s1 - 1) * stirling(n - s1, k - 1)
  })
  
  # Sample the first cluster size s1 using these weights
  s1 <- sample(possible_s1, size = 1, prob = weights)
  
  # Recursively sample the sizes for the remaining k-1 clusters from the remaining items
  remaining <- sample_cluster_sizes(n - s1, k - 1)
  
  # Return the vector of cluster sizes (ordered as they are generated)
  return(c(s1, remaining))
}

# Example usage:
set.seed(123)  # For reproducibility
n <- 10        # Total number of items
k <- 3         # Desired number of clusters

z <- table(sapply(seq_len(10000), \(i) paste0(rev(sort(sample_cluster_sizes(n, k))),collapse="")))
z / sum(z)


n_items <- 10
concentration <- 2
crp <- CRPPartition(nItems = n_items, concentration = concentration)
x <- samplePartition(crp, 100000)
w <- apply(x, 1, \(y) length(unique(y))) == 3
mean(w)
z <- table(apply(x[w,], 1, \(y) paste0(rev(sort(table(y))), collapse="")))
z / sum(z)





