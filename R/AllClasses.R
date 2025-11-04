#' ClonalSimData Class
#'
#' @description S4 class to store results from tumor clonal evolution simulation
#'
#' @slot mutations data.frame containing mutation information with columns:
#'   Mutation, Chromosome, Position, Ref, Alt, True_VAF, VAF, Depth, Alt_reads,
#'   Clone, Type, Clone_IDs
#' @slot params list containing simulation parameters:
#'   subclone_freqs, n_mut_per_clone, n_mut_founder, n_mut_shared,
#'   biological_noise, sequencing_noise
#' @slot clonal_structure data.frame defining hierarchical clone relationships
#' @slot metadata list containing additional metadata (date, version, seed)
#'
#' @name ClonalSimData-class
#' @rdname ClonalSimData-class
#' @exportClass ClonalSimData
#'
#' @examples
#' # Create a simple simulation
#' sim <- simulateTumor(
#'   subclone_freqs = c(0.3, 0.4, 0.3),
#'   n_mut_per_clone = c(30, 40, 30)
#' )
#'
#' # Access mutations
#' head(getMutations(sim))
#'
#' # Access parameters
#' getSimParams(sim)
#'
setClass("ClonalSimData",
         slots = c(
           mutations = "data.frame",
           params = "list",
           clonal_structure = "data.frame",
           metadata = "list"
         ),
         prototype = list(
           mutations = data.frame(),
           params = list(),
           clonal_structure = data.frame(),
           metadata = list()
         )
)

#' Constructor for ClonalSimData
#'
#' @param mutations data.frame with mutation data
#' @param params list with simulation parameters
#' @param clonal_structure data.frame with clone hierarchy
#' @param metadata list with additional metadata
#'
#' @return ClonalSimData object
#' @export
#'
#' @examples
#' # Typically created by simulateTumor(), but can be constructed manually
#' mutations <- data.frame(
#'   Mutation = "mut1",
#'   Chromosome = "chr1",
#'   Position = 1000000,
#'   Ref = "A",
#'   Alt = "T",
#'   True_VAF = 0.5,
#'   VAF = 0.48,
#'   Depth = 100,
#'   Alt_reads = 48,
#'   Clone = "Clone1",
#'   Type = "private"
#' )
#' params <- list(subclone_freqs = c(0.5))
#' sim_data <- ClonalSimData(mutations = mutations, params = params)
#'
ClonalSimData <- function(mutations = data.frame(),
                          params = list(),
                          clonal_structure = data.frame(),
                          metadata = list()) {
  new("ClonalSimData",
      mutations = mutations,
      params = params,
      clonal_structure = clonal_structure,
      metadata = metadata)
}

#' Validity check for ClonalSimData
#'
#' @param object ClonalSimData object
#' @return TRUE if valid, otherwise error message
#' @name ClonalSimData-validity
#' @rdname ClonalSimData-class
#'
setValidity("ClonalSimData", function(object) {
  errors <- character()

  # Check mutations data.frame
  if (nrow(object@mutations) > 0) {
    required_cols <- c("Mutation", "Chromosome", "Position", "Ref", "Alt",
                       "True_VAF", "VAF", "Depth", "Alt_reads", "Clone", "Type")
    missing_cols <- setdiff(required_cols, colnames(object@mutations))
    if (length(missing_cols) > 0) {
      errors <- c(errors, paste("Missing required columns in mutations:",
                                paste(missing_cols, collapse = ", ")))
    }

    # Check VAF bounds
    if (any(object@mutations$VAF < 0 | object@mutations$VAF > 1, na.rm = TRUE)) {
      errors <- c(errors, "VAF values must be between 0 and 1")
    }

    # Check depth is positive
    if (any(object@mutations$Depth <= 0, na.rm = TRUE)) {
      errors <- c(errors, "Depth values must be positive")
    }

    # Check alt_reads <= depth
    if (any(object@mutations$Alt_reads > object@mutations$Depth, na.rm = TRUE)) {
      errors <- c(errors, "Alt_reads cannot exceed Depth")
    }
  }

  # Check params
  if (length(object@params) > 0) {
    if ("subclone_freqs" %in% names(object@params)) {
      freqs <- object@params$subclone_freqs
      if (any(freqs < 0 | freqs > 1)) {
        errors <- c(errors, "Subclone frequencies must be between 0 and 1")
      }
      if (sum(freqs) > 1) {
        errors <- c(errors, "Sum of subclone frequencies cannot exceed 1")
      }
    }
  }

  if (length(errors) == 0) TRUE else errors
})
