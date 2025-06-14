# Function to enumerate all possible size configurations for partitioning n items into k clusters
# Returns configurations in non-increasing order to avoid duplicates
enumerate_partition_configs <- function(n, k) {
  if (n <= 0 || k <= 0) {
    stop("Both n and k must be positive integers")
  }
  
  if (k > n) {
    return(matrix(nrow = 0, ncol = k))  # Cannot partition n items into more than n non-empty clusters
  }
  
  # Recursive helper function to generate partitions in non-increasing order
  generate_partitions <- function(remaining, num_parts, max_allowed, current_partition) {
    # Base case: if we need exactly 1 more part, it must contain all remaining items
    if (num_parts == 1) {
      if (remaining >= 1 && remaining <= max_allowed) {
        return(list(c(current_partition, remaining)))
      } else {
        return(list())
      }
    }
    
    # If remaining items can't fill the remaining parts with at least 1 each
    if (remaining < num_parts) {
      return(list())
    }
    
    results <- list()
    
    # Calculate bounds for current cluster size
    # Must be at least 1
    min_current <- 1
    # Must be at most max_allowed (to maintain non-increasing order)
    # Must leave at least 1 item for each remaining cluster
    max_current <- min(max_allowed, remaining - (num_parts - 1))
    
    # Try sizes from largest to smallest to build non-increasing sequences
    for (size in max_current:min_current) {
      # Recursively generate partitions for the remaining items and parts
      # Next max_allowed is current size (to ensure non-increasing order)
      sub_results <- generate_partitions(
        remaining - size, 
        num_parts - 1, 
        size,  # Next cluster must be <= current size
        c(current_partition, size)
      )
      results <- c(results, sub_results)
    }
    
    return(results)
  }
  
  # Start the recursive generation
  # max_allowed = n (first cluster can have up to all items)
  partitions <- generate_partitions(n, k, n, c())
  
  # Convert to matrix for cleaner output (each row is a configuration)
  if (length(partitions) > 0) {
    do.call(rbind, partitions)
  } else {
    matrix(nrow = 0, ncol = k)
  }
}

# Example usage and testing
# configs <- enumerate_partition_configs(9, 3)
# print(configs)

