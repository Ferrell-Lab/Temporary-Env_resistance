#' Calculate Earth Mover's Distance (EMD)
#'
#' @param f Vector of densities for the first histogram
#' @param m Vector of densities for the second histogram
#' @return The Earth Mover's Distance between the two histograms
#' @export
emdR <- function(f, m) {
  # Normalize the vectors
  f_norm <- f / sum(f)
  m_norm <- m / sum(m)

  # Initialize variables
  n <- length(f_norm)  # Length of the vectors
  d <- numeric(n + 1)  # Create vector d, starting with d0 = 0
  d[1] <- 0  # d0 = 0

  # Iterate through each element to calculate d_i+1
  for (i in 1:n) {
    d[i + 1] <- f_norm[i] + d[i] - m_norm[i]
  }

  # Calculate the EMD as the sum of absolute differences
  emd_value <- sum(abs(d[-1]))  # Exclude d0 when summing |di|

  return(emd_value)
}
