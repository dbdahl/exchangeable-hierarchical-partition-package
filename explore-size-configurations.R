# Function to generate unique size configurations recursively with enforced non-increasing order
generate_size_configurations <- function(n, k, max_val=n) {
  # Base case: If only one cluster left, return the remaining size as a single partition
  if (k == 1) {
    if (n <= max_val) {
      return(list(c(n)))
    } else {
      return(list())  # No valid partition
    }
  }
  
  # List to store valid size configurations
  configurations <- list()
  
  # Iterate from max_val downward to ensure non-increasing order
  for (s in seq(min(n - k + 1, max_val), 1, by=-1)) {
    # Recursively generate partitions for remaining (n-s) into (k-1) groups
    sub_configs <- generate_size_configurations(n - s, k - 1, s)
    
    # Append each result as a new partition
    for (sub in sub_configs) {
      configurations <- append(configurations, list(c(s, sub)))
    }
  }
  
  return(configurations)
}

# Example usage
n <- 60
k <- 10
size_configurations <- generate_size_configurations(n, k)

# Print results in a readable format
print(size_configurations)
length(size_configurations)

size_configurations <- generate_size_configurations(1000, 3)
print(size_configurations)
length(size_configurations)

entropies <- sapply(size_configurations, \(x) { p <- x/sum(x); -sum(p*log(p)) })
plot(density(entropies))


system.time(size_configurations <- generate_size_configurations(1000, 4))
print(size_configurations)
length(size_configurations)

entropies <- sapply(size_configurations, \(x) { p <- x/sum(x); -sum(p*log(p)) })
plot(density(entropies))
