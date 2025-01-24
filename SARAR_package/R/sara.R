#' Calculate SARA Score using Earth Mover's Distance
#'
#' @param unstim_data Numeric vector of unstimulated data
#' @param stim_data Numeric vector of stimulated data
#' @param n_permutations Number of permutations for significance testing
#' @param num_bins Number of bins for histograms
#' @return The SARA score
#' @export
sara <- function(unstim_data, stim_data, n_permutations = 1000, num_bins = 100) {
  common_breaks <- seq(min(c(unstim_data, stim_data)), max(c(unstim_data, stim_data)), length.out = num_bins + 1)

  hist_unstim <- hist(unstim_data, breaks = common_breaks, plot = FALSE)
  hist_stim <- hist(stim_data, breaks = common_breaks, plot = FALSE)

  unstim_density <- hist_unstim$density
  stim_density <- hist_stim$density

  if (length(unstim_density) != length(stim_density)) stop("Histograms must have the same number of bins.")

  emd_value <- emdR(unstim_density, stim_density)

  emd_permutations <- numeric(n_permutations)
  for (i in 1:n_permutations) {
    permuted_unstim <- sample(unstim_data, length(unstim_data), replace = TRUE)
    permuted_stim <- sample(stim_data, length(stim_data), replace = TRUE)

    hist_permuted_unstim <- hist(permuted_unstim, breaks = common_breaks, plot = FALSE)
    hist_permuted_stim <- hist(permuted_stim, breaks = common_breaks, plot = FALSE)

    permuted_unstim_density <- hist_permuted_unstim$density
    permuted_stim_density <- hist_permuted_stim$density

    emd_permutations[i] <- emdR(permuted_unstim_density, permuted_stim_density)
  }

  num_extreme <- sum(emd_permutations >= emd_value)
  p_value <- num_extreme / n_permutations

  mean_diff <- mean(stim_data) - mean(unstim_data)
  SARA_score <- emd_value * sign(mean_diff) * (1 - p_value)

  return(SARA_score)
}
