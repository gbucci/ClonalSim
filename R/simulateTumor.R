#' Simulate tumor clonal evolution with realistic noise
#'
#' @description Main function to simulate a heterogeneous tumor sample with
#'   hierarchical clonal structure, including realistic biological and technical
#'   sequencing noise.
#'
#' @param subclone_freqs numeric vector, frequencies of each subclone (must sum to <= 1).
#'   Default: c(0.15, 0.25, 0.30, 0.30)
#' @param n_mut_per_clone integer vector, number of private mutations per clone.
#'   Default: c(20, 25, 30, 15)
#' @param n_mut_founder integer, number of founder mutations (present in all clones).
#'   Default: 10
#' @param n_mut_shared named list, shared mutations between clone groups.
#'   Names indicate which clones (e.g., "2 3 4"), values are mutation counts.
#'   Default: list("2 3 4" = 15, "3 4" = 12, "1 2" = 8)
#' @param biological_noise list with biological noise parameters:
#'   \itemize{
#'     \item enabled: logical, enable biological noise (default: TRUE)
#'     \item concentration: numeric, Beta distribution concentration (default: 50)
#'   }
#' @param sequencing_noise list with sequencing noise parameters:
#'   \itemize{
#'     \item enabled: logical, enable sequencing noise (default: TRUE)
#'     \item mean_depth: numeric, average coverage (default: 100)
#'     \item depth_variation: character, distribution type (default: "negative_binomial")
#'     \item depth_dispersion: numeric, dispersion parameter (default: 20)
#'     \item error_rate: numeric, base miscall rate (default: 0.001)
#'     \item binomial_sampling: logical, use binomial sampling (default: TRUE)
#'   }
#' @param germline_variants list with germline variant parameters:
#'   \itemize{
#'     \item enabled: logical, include germline variants (default: FALSE)
#'     \item n_variants: integer, number of germline variants to simulate (default: 50)
#'     \item vaf_expected: numeric, expected VAF for heterozygous germline (default: 0.5)
#'   }
#'   Germline variants will have VAF adjusted by tumor purity:
#'   observed_VAF = germline_vaf * tumor_purity + 0.5 * (1 - tumor_purity)
#' @return ClonalSimData object containing mutations, parameters, and metadata
#' @note For reproducibility, set the random seed before calling this function
#'   using \code{set.seed()}.
#' @export
#'
#' @examples
#' # Basic simulation with default parameters
#' sim <- simulateTumor()
#' show(sim)
#' plot(sim)
#'
#' # Custom clonal structure
#' sim <- simulateTumor(
#'   subclone_freqs = c(0.2, 0.3, 0.5),
#'   n_mut_per_clone = c(30, 40, 30)
#' )
#'
#' # Low purity tumor (40% purity)
#' sim <- simulateTumor(
#'   subclone_freqs = c(0.1, 0.15, 0.15),
#'   n_mut_per_clone = c(20, 25, 30)
#' )
#'
#' # High coverage sequencing
#' sim <- simulateTumor(
#'   sequencing_noise = list(
#'     enabled = TRUE,
#'     mean_depth = 500,
#'     depth_dispersion = 100
#'   )
#' )
#'
#' # No noise (ideal data for testing)
#' sim <- simulateTumor(
#'   biological_noise = list(enabled = FALSE),
#'   sequencing_noise = list(enabled = FALSE)
#' )
#'
#' # Include germline variants
#' sim <- simulateTumor(
#'   subclone_freqs = c(0.3, 0.4),  # 70% purity
#'   n_mut_per_clone = c(40, 50),
#'   germline_variants = list(enabled = TRUE, n_variants = 100)
#' )
#'
simulateTumor <- function(
    subclone_freqs = c(0.15, 0.25, 0.30, 0.30),
    n_mut_per_clone = c(20, 25, 30, 15),
    n_mut_founder = 10,
    n_mut_shared = list("2 3 4" = 15, "3 4" = 12, "1 2" = 8),
    biological_noise = list(enabled = TRUE, concentration = 50),
    sequencing_noise = list(
      enabled = TRUE,
      mean_depth = 100,
      depth_variation = "negative_binomial",
      depth_dispersion = 20,
      error_rate = 0.001,
      binomial_sampling = TRUE
    ),
    germline_variants = list(
      enabled = FALSE,
      n_variants = 50,
      vaf_expected = 0.5
    )
) {

  # Input validation
  if (any(subclone_freqs < 0) || any(subclone_freqs > 1)) {
    stop("subclone_freqs must be between 0 and 1")
  }
  if (sum(subclone_freqs) > 1) {
    stop("Sum of subclone_freqs cannot exceed 1")
  }
  if (length(n_mut_per_clone) != length(subclone_freqs)) {
    stop("n_mut_per_clone must have same length as subclone_freqs")
  }

  # Store clone names
  names(subclone_freqs) <- paste0("Clone", seq_along(subclone_freqs))

  # Fill in default noise parameters if not provided
  if (!"enabled" %in% names(biological_noise)) {
    biological_noise$enabled <- TRUE
  }
  if (!"concentration" %in% names(biological_noise)) {
    biological_noise$concentration <- 50
  }

  if (!"enabled" %in% names(sequencing_noise)) {
    sequencing_noise$enabled <- TRUE
  }
  if (!"mean_depth" %in% names(sequencing_noise)) {
    sequencing_noise$mean_depth <- 100
  }
  if (!"depth_variation" %in% names(sequencing_noise)) {
    sequencing_noise$depth_variation <- "negative_binomial"
  }
  if (!"depth_dispersion" %in% names(sequencing_noise)) {
    sequencing_noise$depth_dispersion <- 20
  }
  if (!"error_rate" %in% names(sequencing_noise)) {
    sequencing_noise$error_rate <- 0.001
  }
  if (!"binomial_sampling" %in% names(sequencing_noise)) {
    sequencing_noise$binomial_sampling <- TRUE
  }

  # Generate mutations
  mutation_list <- list()
  idx <- 1

  # 1. Founder mutations (present in all clones)
  founder_freq <- sum(subclone_freqs)
  mutation_list[[idx]] <- .generate_mutations(
    n_mut_founder,
    founder_freq,
    clone_ids = seq_along(subclone_freqs),
    type = "founder",
    bio_noise = biological_noise,
    seq_noise = sequencing_noise
  )
  idx <- idx + 1

  # 2. Shared mutations between subgroups
  for (i in seq_along(n_mut_shared)) {
    shared_clones <- as.numeric(strsplit(names(n_mut_shared)[i], " ")[[1]])

    # Skip if any clone index is out of bounds
    if (any(shared_clones > length(subclone_freqs))) {
      warning("Shared mutation group '", names(n_mut_shared)[i],
              "' references non-existent clones. Skipping.")
      next
    }

    shared_freq <- sum(subclone_freqs[shared_clones])

    mutation_list[[idx]] <- .generate_mutations(
      n_mut_shared[[i]],
      shared_freq,
      clone_ids = shared_clones,
      type = "shared",
      bio_noise = biological_noise,
      seq_noise = sequencing_noise
    )
    idx <- idx + 1
  }

  # 3. Private mutations for each subclone
  for (i in seq_along(subclone_freqs)) {
    mutation_list[[idx]] <- .generate_mutations(
      n_mut_per_clone[i],
      subclone_freqs[i],
      clone_ids = i,
      type = "private",
      bio_noise = biological_noise,
      seq_noise = sequencing_noise
    )
    idx <- idx + 1
  }

  # 4. Germline variants (if enabled)
  if (!is.null(germline_variants) && germline_variants$enabled) {
    # Fill in defaults
    if (!"n_variants" %in% names(germline_variants)) {
      germline_variants$n_variants <- 50
    }
    if (!"vaf_expected" %in% names(germline_variants)) {
      germline_variants$vaf_expected <- 0.5
    }

    # Calculate tumor purity
    tumor_purity <- sum(subclone_freqs)

    # Germline variants in tumor sample:
    # - In pure tumor cells: VAF = 0.5 (heterozygous)
    # - In normal cells: VAF = 0.5 (heterozygous)
    # - Mixed sample: VAF = 0.5 * tumor_purity + 0.5 * (1 - tumor_purity) = 0.5
    # So germline VAF is always ~0.5 regardless of purity!
    # (This is the key difference from somatic variants)

    germline_vaf <- germline_variants$vaf_expected

    mutation_list[[idx]] <- .generate_mutations(
      germline_variants$n_variants,
      germline_vaf,
      clone_ids = "germline",
      type = "germline",
      bio_noise = biological_noise,
      seq_noise = sequencing_noise
    )
    idx <- idx + 1
  }

  # Combine all mutations
  mutation_list <- mutation_list[!vapply(mutation_list, is.null, logical(1))]
  complete_data <- do.call(rbind, mutation_list)

  # Add genomic coordinates (random, for demonstration)
  complete_data$Chromosome <- sample(paste0("chr", 1:22),
                                      nrow(complete_data), replace = TRUE)
  complete_data$Position <- sample(1e6:2e8, nrow(complete_data), replace = TRUE)
  complete_data$Ref <- sample(c("A", "T", "C", "G"),
                               nrow(complete_data), replace = TRUE)
  complete_data$Alt <- sample(c("A", "T", "C", "G"),
                               nrow(complete_data), replace = TRUE)

  # Reorder columns
  complete_data <- complete_data[, c("Mutation", "Chromosome", "Position",
                                      "Ref", "Alt", "True_VAF", "VAF",
                                      "Depth", "Alt_reads",
                                      "Clone", "Type", "Clone_IDs")]

  # Create clonal structure data frame
  clonal_structure <- data.frame(
    Clone = names(subclone_freqs),
    Frequency = subclone_freqs,
    N_private_mutations = n_mut_per_clone,
    row.names = NULL
  )

  # Store parameters
  params <- list(
    subclone_freqs = subclone_freqs,
    n_mut_per_clone = n_mut_per_clone,
    n_mut_founder = n_mut_founder,
    n_mut_shared = n_mut_shared,
    biological_noise = biological_noise,
    sequencing_noise = sequencing_noise,
    germline_variants = germline_variants
  )

  # Create metadata
  metadata <- list(
    date = as.character(Sys.time()),
    version = as.character(packageVersion("ClonalSim")),
    seed = seed
  )

  # Create and return ClonalSimData object
  result <- ClonalSimData(
    mutations = complete_data,
    params = params,
    clonal_structure = clonal_structure,
    metadata = metadata
  )

  return(result)
}

#' Internal function to generate mutations (not exported)
#'
#' @keywords internal
#' @noRd
.generate_mutations <- function(n_mut, base_freq, clone_ids, type = "private",
                                bio_noise, seq_noise) {

  if (n_mut == 0) return(NULL)

  # Generate mutation names
  if (type == "founder") {
    mut_names <- paste0("Founder_", seq_len(n_mut))
    clone_label <- "Founder"
  } else if (type == "shared") {
    clone_str <- paste(clone_ids, collapse = "_")
    mut_names <- paste0("Shared_C", clone_str, "_mut", seq_len(n_mut))
    clone_label <- paste0("Clone", paste(clone_ids, collapse = "+"))
  } else if (type == "germline") {
    mut_names <- paste0("Germline_", seq_len(n_mut))
    clone_label <- "Germline"
    clone_ids <- "germline"  # Ensure it's a string
  } else {
    mut_names <- paste0("Clone", clone_ids, "_mut", seq_len(n_mut))
    clone_label <- paste0("Clone", clone_ids)
  }

  # Step 1: Apply biological noise
  if (bio_noise$enabled) {
    true_vaf <- applyBiologicalNoise(base_freq, n_mut, bio_noise$concentration)
  } else {
    true_vaf <- rep(base_freq, n_mut)
  }

  # Step 2: Simulate depth
  if (seq_noise$enabled) {
    depth <- simulateDepth(
      n_mut,
      mean_depth = seq_noise$mean_depth,
      distribution = seq_noise$depth_variation,
      dispersion = seq_noise$depth_dispersion
    )
  } else {
    depth <- rep(100L, n_mut)
  }

  # Step 3: Apply technical sequencing noise
  if (seq_noise$enabled && seq_noise$binomial_sampling) {
    seq_result <- simulateSequencingReads(true_vaf, depth, seq_noise$error_rate)
    observed_vaf <- seq_result$vaf
    alt_reads <- seq_result$alt_reads
  } else {
    observed_vaf <- true_vaf
    alt_reads <- as.integer(round(true_vaf * depth))
  }

  # Create data frame
  data.frame(
    Mutation = mut_names,
    True_VAF = true_vaf,
    VAF = observed_vaf,
    Depth = depth,
    Alt_reads = alt_reads,
    Clone = clone_label,
    Type = type,
    Clone_IDs = paste(clone_ids, collapse = ","),
    stringsAsFactors = FALSE
  )
}
