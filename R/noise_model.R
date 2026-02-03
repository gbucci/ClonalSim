#' Apply biological noise using Beta distribution
#'
#' @description Simulates intra-tumor heterogeneity using a Beta distribution.
#'   The Beta distribution is more realistic than Gaussian for frequency data
#'   as it is naturally bounded between 0 and 1.
#'
#' @param true_freq numeric, the true allele frequency (0-1)
#' @param n_mutations integer, number of mutations to generate
#' @param concentration numeric, concentration parameter for Beta distribution
#'   (alpha + beta). Higher values = less variability. Default: 50
#'
#' @return numeric vector of VAF values with biological noise
#' @export
#'
#' @details
#' For a clone with true frequency \code{f}, VAF values are sampled from
#' Beta(alpha, beta) where:
#' \itemize{
#'   \item alpha = f * concentration
#'   \item beta = (1-f) * concentration
#' }
#'
#' Higher concentration values result in VAF closer to the expected frequency,
#' while lower values increase spread, simulating spatial/temporal heterogeneity.
#'
#' @examples
#' # Generate 100 mutations with true frequency 0.3 and moderate heterogeneity
#' vaf <- applyBiologicalNoise(0.3, 100, concentration = 50)
#' hist(vaf, main = "VAF with biological noise", xlab = "VAF")
#'
#' # High heterogeneity (low concentration)
#' vaf_high_het <- applyBiologicalNoise(0.3, 100, concentration = 20)
#'
#' # Low heterogeneity (high concentration)
#' vaf_low_het <- applyBiologicalNoise(0.3, 100, concentration = 100)
#'
applyBiologicalNoise <- function(true_freq, n_mutations, concentration = 50) {
  if (n_mutations == 0) return(numeric(0))
  if (true_freq < 0 || true_freq > 1) {
    stop("true_freq must be between 0 and 1")
  }
  if (concentration <= 0) {
    stop("concentration must be positive")
  }

  # Beta distribution parameters
  alpha <- true_freq * concentration
  beta <- (1 - true_freq) * concentration

  # Avoid degenerate cases
  alpha <- pmax(alpha, 0.1)
  beta <- pmax(beta, 0.1)

  # Sample from Beta distribution
  noisy_freq <- rbeta(n_mutations, shape1 = alpha, shape2 = beta)

  # Limit to reasonable range
  noisy_freq <- pmax(0.01, pmin(0.99, noisy_freq))

  return(noisy_freq)
}

#' Simulate realistic sequencing depth
#'
#' @description Simulates sequencing coverage with realistic overdispersion
#'   using negative binomial distribution, which better captures real NGS
#'   variability compared to Poisson.
#'
#' @param n_mutations integer, number of mutations
#' @param mean_depth numeric, mean sequencing coverage. Default: 100
#' @param distribution character, distribution type: "negative_binomial"
#'   (realistic, default), "poisson" (simpler), or "uniform" (for testing)
#' @param dispersion numeric, size parameter for negative binomial. Lower values
#'   = more variable coverage. Default: 20
#'
#' @return integer vector of depth values
#' @export
#'
#' @details
#' The negative binomial distribution allows overdispersion (variance > mean),
#' which occurs in real sequencing due to:
#' \itemize{
#'   \item GC content bias
#'   \item Mappability differences
#'   \item PCR amplification artifacts
#' }
#'
#' The dispersion parameter controls variability: lower values produce more
#' variable coverage across positions.
#'
#' @examples
#' # Realistic depth with overdispersion
#' depth_realistic <- simulateDepth(1000, mean_depth = 100, dispersion = 20)
#' hist(depth_realistic, main = "Realistic depth", xlab = "Depth")
#'
#' # More variable coverage (low dispersion)
#' depth_variable <- simulateDepth(1000, mean_depth = 100, dispersion = 10)
#'
#' # Uniform coverage (for testing)
#' depth_uniform <- simulateDepth(1000, mean_depth = 100, distribution = "uniform")
#'
simulateDepth <- function(n_mutations, mean_depth = 100,
                          distribution = "negative_binomial",
                          dispersion = 20) {
  if (mean_depth <= 0) {
    stop("mean_depth must be positive")
  }

  if (distribution == "negative_binomial") {
    if (dispersion <= 0) {
      stop("dispersion must be positive")
    }
    depths <- rnbinom(n_mutations, size = dispersion, mu = mean_depth)

  } else if (distribution == "poisson") {
    depths <- rpois(n_mutations, lambda = mean_depth)

  } else if (distribution == "uniform") {
    depths <- rep(as.integer(mean_depth), n_mutations)

  } else {
    stop("distribution must be 'negative_binomial', 'poisson', or 'uniform'")
  }

  # Ensure minimum depth of 10
  depths <- pmax(depths, 10L)

  return(as.integer(depths))
}

#' Simulate sequencing reads with binomial sampling
#'
#' @description Simulates stochastic read counts using binomial distribution,
#'   which models the independent sampling of each sequencing read.
#'
#' @param true_vaf numeric vector, true variant allele frequencies
#' @param depth integer vector, sequencing depth at each position
#' @param error_rate numeric, base miscall rate (0-1). Default: 0.001 (0.1%,
#'   typical for Illumina)
#'
#' @return list with three components:
#'   \itemize{
#'     \item vaf: observed VAF values
#'     \item alt_reads: alternative allele read counts
#'     \item ref_reads: reference allele read counts
#'   }
#' @export
#'
#' @details
#' Each sequencing read is independently sampled from the DNA pool. Given a
#' true VAF and depth D, the number of alternative reads follows a binomial
#' distribution: alt_reads ~ Binomial(D, VAF).
#'
#' This introduces natural sampling variation, with higher uncertainty at:
#' \itemize{
#'   \item Low sequencing depth
#'   \item Low variant frequency
#' }
#'
#' The error_rate parameter simulates base miscalls from sequencer optical/
#' chemistry errors.
#'
#' @examples
#' # Simulate reads for mutations at VAF 0.3 with depth 100
#' true_vaf <- rep(0.3, 100)
#' depth <- rep(100, 100)
#' reads <- simulateSequencingReads(true_vaf, depth)
#'
#' # Observed VAF varies due to sampling
#' hist(reads$vaf, main = "Observed VAF", xlab = "VAF")
#' abline(v = 0.3, col = "red", lty = 2)
#'
#' # Check alt_reads distribution
#' hist(reads$alt_reads, main = "Alt read counts", xlab = "Alt reads")
#'
simulateSequencingReads <- function(true_vaf, depth, error_rate = 0.001) {
  if (length(true_vaf) != length(depth)) {
    stop("true_vaf and depth must have the same length")
  }
  if (any(true_vaf < 0 | true_vaf > 1, na.rm = TRUE)) {
    stop("true_vaf values must be between 0 and 1")
  }
  if (any(depth <= 0, na.rm = TRUE)) {
    stop("depth values must be positive")
  }
  if (error_rate < 0 || error_rate > 1) {
    stop("error_rate must be between 0 and 1")
  }

  # Observed VAF includes sequencing errors
  observed_vaf <- true_vaf + error_rate
  observed_vaf <- pmin(observed_vaf, 0.99)

  # Sample alt reads from binomial distribution
  alt_reads <- rbinom(length(true_vaf), size = depth, prob = observed_vaf)

  # Calculate observed VAF from read counts
  final_vaf <- alt_reads / depth

  # Ensure alt_reads doesn't exceed depth (should never happen with rbinom, but safety check)
  alt_reads <- pmin(alt_reads, depth)

  return(list(
    vaf = final_vaf,
    alt_reads = as.integer(alt_reads),
    ref_reads = as.integer(depth - alt_reads)
  ))
}
