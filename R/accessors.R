#' Accessor functions for ClonalSimData objects
#'
#' @name accessors
#' @rdname accessors
NULL

#' Get mutations data frame
#'
#' @param object ClonalSimData object
#' @return data.frame with mutation information
#' @export
#'
#' @examples
#' sim <- simulateTumor()
#' mutations <- getMutations(sim)
#' head(mutations)
#'
getMutations <- function(object) {
  stopifnot(inherits(object, "ClonalSimData"))
  object@mutations
}

#' Get simulation parameters
#'
#' @param object ClonalSimData object
#' @return list with simulation parameters
#' @export
#'
#' @examples
#' sim <- simulateTumor()
#' params <- getSimParams(sim)
#' params$subclone_freqs
#'
getSimParams <- function(object) {
  stopifnot(inherits(object, "ClonalSimData"))
  object@params
}

#' Get true VAF values (biological, without sequencing noise)
#'
#' @param object ClonalSimData object
#' @return numeric vector of true VAF values
#' @export
#'
#' @examples
#' sim <- simulateTumor()
#' true_vaf <- getTrueVAF(sim)
#' summary(true_vaf)
#'
getTrueVAF <- function(object) {
  stopifnot(inherits(object, "ClonalSimData"))
  object@mutations$True_VAF
}

#' Get observed VAF values (with sequencing noise)
#'
#' @param object ClonalSimData object
#' @return numeric vector of observed VAF values
#' @export
#'
#' @examples
#' sim <- simulateTumor()
#' obs_vaf <- getObservedVAF(sim)
#' summary(obs_vaf)
#'
getObservedVAF <- function(object) {
  stopifnot(inherits(object, "ClonalSimData"))
  object@mutations$VAF
}

#' Get clonal structure
#'
#' @param object ClonalSimData object
#' @return data.frame defining clone hierarchy
#' @export
#'
#' @examples
#' sim <- simulateTumor()
#' structure <- getClonalStructure(sim)
#' print(structure)
#'
getClonalStructure <- function(object) {
  stopifnot(inherits(object, "ClonalSimData"))
  object@clonal_structure
}

#' Get metadata
#'
#' @param object ClonalSimData object
#' @return list with metadata
#' @export
#'
#' @examples
#' sim <- simulateTumor()
#' metadata <- getMetadata(sim)
#' metadata$date
#'
getMetadata <- function(object) {
  stopifnot(inherits(object, "ClonalSimData"))
  object@metadata
}
